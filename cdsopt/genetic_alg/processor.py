# -*- coding: utf-8 -*-
"""GA Processor: main loop, pool reuse, checkpointing, convergence."""
from __future__ import annotations

import datetime
import logging
import pickle
import random
from dataclasses import dataclass, field
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Dict, List, Optional

from tqdm import tqdm

from cdsopt.fitness.evaluator import FitnessConfig, FitnessEvaluator
from cdsopt.genetic_alg.individual import Individual, ProteinSpec
from cdsopt.genetic_alg.nsga2 import Objective, environmental_selection, fast_non_dominated_sort, _count_satisfied
from cdsopt.genetic_alg.operators import GeneticOperators

logger = logging.getLogger(__name__)
CHECKPOINT_VERSION = "1.0"


@dataclass
class GAConfig:
    population_size: int = 100
    generations: int = 1000
    mute_rate: float = 0.01
    n_elite: int = 2
    amplification: int = 3
    processes: int = 1
    output_interval: int = 10
    fitness_config: FitnessConfig = field(default_factory=FitnessConfig)
    not_mutate_idx: set[int] = field(default_factory=set)
    rng_seed: int | None = None
    weighted_init: bool = True


class GeneticAlgorithmProcessor:
    def __init__(self, protein_sequence: str, config: GAConfig | None = None, init_cds_list: List[str] | None = None):
        self.cfg = config or GAConfig()
        self.rng = random.Random(self.cfg.rng_seed)
        # Load codon weights for weighted initialization
        codon_weights = None
        if self.cfg.weighted_init:
            from cdsopt.tables.codon_frequency_table import get_table_weights
            from cdsopt.tables.genetic_code import get_code_map_by_genetic_code
            try:
                raw_weights = get_table_weights(self.cfg.fitness_config.species)
                code_map = get_code_map_by_genetic_code(self.cfg.fitness_config.genetic_code)
                cai_weights = {}
                for aa, codons in code_map.items():
                    freqs = [raw_weights.get(c, 0.0) for c in codons]
                    max_freq = max(freqs) if any(freqs) else 1.0
                    for c, f in zip(codons, freqs):
                        rna_codon = c.replace("T", "U")
                        cai_weight = f / max_freq if max_freq > 0 else 0.0
                        # Sharpen distribution to push initial CAI > 0.95
                        cai_weights[rna_codon] = cai_weight ** 10
                codon_weights = cai_weights
            except Exception as e:
                logger.warning("Failed to load codon weights for weighted init: %s", e)
        self.spec = ProteinSpec.from_protein(protein_sequence, codon_weights=codon_weights)
        from cdsopt.fitness.evaluator import build_objectives
        self.objectives = build_objectives(self.cfg.fitness_config)
        self.evaluator = FitnessEvaluator(config=self.cfg.fitness_config)
        self.operators = GeneticOperators(self.spec, base_mute_rate=self.cfg.mute_rate)
        self.population: List[Individual] = []
        self.fitness_list: List[dict] = []
        self.generation = 0
        self._pool: Optional[Pool] = None
        if self.cfg.processes > 1:
            self._pool = Pool(processes=self.cfg.processes)
        self._fixed_ref: Individual | None = None
        self._initialize_population(init_cds_list)

    def _initialize_population(self, init_cds_list: List[str] | None) -> None:
        n_init = 0
        if init_cds_list:
            for cds in init_cds_list[: self.cfg.population_size]:
                try:
                    cds_rna = cds.upper().replace("T", "U")
                    ind = Individual.from_codon_list([cds_rna[i : i + 3] for i in range(0, len(cds_rna), 3)], self.spec)
                    self.population.append(ind)
                    n_init += 1
                except ValueError as e:
                    logger.warning("Skipping invalid init CDS: %s", e)
            logger.info("Loaded %d initial CDS sequences from input", n_init)
        if self.population:
            self._fixed_ref = self.population[0].copy()
        else:
            self._fixed_ref = self.spec.random_individual(self.rng, weighted=True)
        n_random = self.cfg.population_size - len(self.population)
        if n_random > 0:
            logger.info("Generating %d random individuals to fill population (weighted=%s)", n_random, self.cfg.weighted_init)
            for _ in range(n_random):
                self.population.append(self.spec.random_individual(
                    self.rng, weighted=self.cfg.weighted_init,
                    reference=self._fixed_ref, fixed_idx=self.cfg.not_mutate_idx
                ))
        logger.info("Evaluating initial population of %d individuals...", len(self.population))
        self._evaluate()
        init_unique = len(set(self.population))
        logger.info("Initial population evaluation complete. unique=%d, duplicates=%d", init_unique, len(self.population) - init_unique)

    def run(self, output_dir: str | Path = "./outputs") -> None:
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        log_file = out_path / f"{datetime.datetime.now():%Y%m%d_%H%M%S}.log"
        prev_seq_0 = None
        pbar = tqdm(range(self.generation, self.cfg.generations), desc="GA", disable=not __import__("sys").stdout.isatty())
        for gen in pbar:
            self.generation = gen
            self._step()
            curr_seq_0 = self.spec.to_rna(self.population[0])
            if curr_seq_0 != prev_seq_0 or gen % self.cfg.output_interval == 0:
                self._write_generation_fitness(out_path, gen)
                prev_seq_0 = curr_seq_0
            if gen % self.cfg.output_interval == 0:
                self._write_logs(log_file)
                self._save_checkpoint(out_path / "checkpoint.pkl")
            pbar.set_description(f"GA Gen {gen}")
            pbar.set_postfix(self._summary_dict())
        self._write_logs(log_file)
        self._save_checkpoint(out_path / "checkpoint.pkl")
        self._shutdown_pool()

    def _step(self) -> None:
        elites = self.operators.select_elite(self.population, self.fitness_list, self.objectives, self.cfg.n_elite)
        children = self._reproduce()
        combined_pop = elites + children + self.population
        combined_fitness = self._evaluate_population(combined_pop)
        selected = environmental_selection(combined_fitness, self.cfg.population_size, self.objectives)
        self.population = [combined_pop[i].copy() for i in selected]
        self.fitness_list = [combined_fitness[i] for i in selected]
        pre_dedup_unique = len(set(self.population))
        logger.debug("Gen %d post-selection: pop=%d, unique=%d", self.generation, len(self.population), pre_dedup_unique)

        # Deduplicate population and backfill with random individuals to maintain size
        unique_pop: List[Individual] = []
        unique_fitness: List[dict] = []
        seen: set[Individual] = set()
        for ind, fit in zip(self.population, self.fitness_list):
            if ind not in seen:
                unique_pop.append(ind)
                unique_fitness.append(fit)
                seen.add(ind)
        n_backfill = self.cfg.population_size - len(unique_pop)
        if n_backfill > 0:
            logger.debug("Gen %d dedup removed %d duplicates, backfilling %d random individuals", self.generation, len(self.population) - len(unique_pop), n_backfill)
        while len(unique_pop) < self.cfg.population_size:
            new_ind = self.spec.random_individual(
                rng=self.rng, weighted=self.cfg.weighted_init,
                reference=self._fixed_ref, fixed_idx=self.cfg.not_mutate_idx
            )
            unique_pop.append(new_ind)
            unique_fitness.append(self.evaluator.evaluate(self.spec.to_rna(new_ind)))
        self.population = unique_pop
        self.fitness_list = unique_fitness
        logger.debug("Gen %d step complete: pop=%d, unique=%d, elite=%d, children=%d", self.generation, len(self.population), len(set(self.population)), len(elites), len(children))

    def _reproduce(self) -> List[Individual]:
        n_children = self.cfg.population_size * self.cfg.amplification
        mute_rate = self.operators.adaptive_mute_rate(self.population)
        children: List[Individual] = []
        existing = set(self.population)
        diversity_limit_reached = len(existing) * 2 >= self.spec.total_variants
        attempts = 0
        skipped = 0
        while len(children) < n_children:
            attempts += 1
            p1 = self.rng.choice(self.population)
            p2 = self.rng.choice(self.population)
            child = self.operators.single_point_crossover(p1, p2, not_mutate_idx=self.cfg.not_mutate_idx, rng=self.rng)
            child = self.operators.mutate(child, mute_rate=mute_rate, not_mutate_idx=self.cfg.not_mutate_idx, rng=self.rng)
            if not diversity_limit_reached and child in existing:
                skipped += 1
                continue
            children.append(child)
            existing.add(child)
        children_unique = len(set(children))
        logger.debug("Gen %d reproduce: attempts=%d, skipped=%d, children=%d, children_unique=%d, mute_rate=%.4f, diversity_limit=%s",
                     self.generation, attempts, skipped, len(children), children_unique, mute_rate, diversity_limit_reached)
        return children

    def _evaluate(self) -> None:
        self.fitness_list = self._evaluate_population(self.population)

    def _evaluate_population(self, population: List[Individual]) -> List[dict]:
        seqs = [self.spec.to_rna(ind) for ind in population]
        if self.cfg.processes <= 1 or len(seqs) < 4:
            results = []
            for _, seq in enumerate(seqs):
                results.append(self.evaluator.evaluate(seq))
            return results
        chunk_size = max(1, len(seqs) // self.cfg.processes)
        chunks = [seqs[i : i + chunk_size] for i in range(0, len(seqs), chunk_size)]
        results = self._pool.map(_eval_chunk, [(chunk, self.cfg.fitness_config) for chunk in chunks])
        flat: Dict[str, dict] = {}
        for r in results:
            flat.update(r)
        return [flat[seq] for seq in seqs]

    def _summary_dict(self) -> dict:
        if not self.fitness_list:
            return {}
        best = self.fitness_list[0]
        result = {obj.key: f"{best.get(obj.key, 0):.3f}" for obj in self.objectives}
        result["Satisfied"] = f"{_count_satisfied(best, self.objectives)}/{len(self.objectives)}"
        return result

    def _write_logs(self, log_file: Path) -> None:
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"Gen {self.generation}: {self._summary_dict()}\n")

    def _write_generation_fitness(self, out_path: Path, gen: int) -> None:
        import csv
        csv_path = out_path / f"fitness_gen_{gen:03d}.csv"
        fieldnames = ["id", "sequence", "length", "full_length"] + [obj.key for obj in self.objectives]
        has_context = self.cfg.fitness_config.prefix or self.cfg.fitness_config.suffix
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for i, (ind, fit) in enumerate(zip(self.population, self.fitness_list)):
                cds = self.spec.to_rna(ind)
                row = {"id": f"seq_{i}", "sequence": cds, "length": len(cds)}
                if has_context:
                    row["full_length"] = fit.get("full_length", len(cds))
                else:
                    row["full_length"] = len(cds)
                row.update({k: fit.get(k, "") for k in fieldnames[4:]})
                writer.writerow(row)
        logger.info("Generation %d fitness written to %s", gen, csv_path)

    def _save_checkpoint(self, path: Path) -> None:
        data = {
            "version": CHECKPOINT_VERSION,
            "protein_sequence": self.spec.protein_sequence,
            "generation": self.generation,
            "population": [ind.indices for ind in self.population],
            "fitness_list": self.fitness_list,
            "config": self.cfg,
            "fixed_ref": self._fixed_ref.indices if self._fixed_ref else None,
            "rng_state": self.rng.getstate(),

        }
        with open(path, "wb") as f:
            pickle.dump(data, f)
        logger.debug("Checkpoint saved to %s", path)

    @classmethod
    def load_checkpoint(cls, path: str | Path, processes: int | None = None) -> "GeneticAlgorithmProcessor":
        with open(path, "rb") as f:
            data = pickle.load(f)
        if data.get("version") != CHECKPOINT_VERSION:
            raise ValueError(f"Checkpoint version mismatch: expected {CHECKPOINT_VERSION}, got {data.get('version')}")
        cfg = data["config"]
        if processes is not None:
            cfg.processes = processes
        instance = cls.__new__(cls)
        instance.cfg = cfg
        instance.rng = random.Random()
        instance.rng.setstate(data["rng_state"])
        instance.spec = ProteinSpec.from_protein(data["protein_sequence"])
        from cdsopt.fitness.evaluator import build_objectives
        instance.objectives = build_objectives(cfg.fitness_config)
        instance.evaluator = FitnessEvaluator(config=cfg.fitness_config)
        instance.operators = GeneticOperators(instance.spec, base_mute_rate=cfg.mute_rate)
        instance.population = [Individual(idx) for idx in data["population"]]
        instance.fitness_list = data["fitness_list"]
        instance.generation = data["generation"]
        fixed_ref_indices = data.get("fixed_ref")
        instance._fixed_ref = Individual(fixed_ref_indices) if fixed_ref_indices is not None else None

        instance._pool = None
        if cfg.processes > 1:
            instance._pool = Pool(processes=cfg.processes)
        return instance

    def _shutdown_pool(self) -> None:
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
            self._pool = None

    def __del__(self):
        if hasattr(self, "_pool"):
            self._shutdown_pool()


def _eval_chunk(args: tuple) -> Dict[str, dict]:
    seqs, cfg = args
    from cdsopt.fitness.evaluator import FitnessEvaluator
    ev = FitnessEvaluator(config=cfg)
    return {seq: ev.evaluate(seq) for seq in seqs}
