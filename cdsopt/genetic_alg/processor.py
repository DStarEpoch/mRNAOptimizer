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
from cdsopt.genetic_alg.nsga2 import Objective, environmental_selection, fast_non_dominated_sort
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
    output_interval: int = 1
    early_stop_patience: int = 50
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
        self._initialize_population(init_cds_list)
        self._best_front_history: List[dict] = []
        self._stagnation_count = 0

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
        n_random = self.cfg.population_size - len(self.population)
        if n_random > 0:
            logger.info("Generating %d random individuals to fill population (weighted=%s)", n_random, self.cfg.weighted_init)
            for _ in range(n_random):
                self.population.append(self.spec.random_individual(self.rng, weighted=self.cfg.weighted_init))
        logger.info("Evaluating initial population of %d individuals...", len(self.population))
        self._evaluate()
        logger.info("Initial population evaluation complete")

    def run(self, output_dir: str | Path = "./outputs") -> None:
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        log_file = out_path / f"{datetime.datetime.now():%Y%m%d_%H%M%S}.log"
        pbar = tqdm(range(self.generation, self.cfg.generations), desc="GA", disable=not __import__("sys").stdout.isatty())
        for gen in pbar:
            self.generation = gen
            self._step()
            if gen % self.cfg.output_interval == 0:
                self._write_logs(log_file)
                self._save_checkpoint(out_path / "checkpoint.pkl")
            if self._check_convergence():
                logger.info("Early stopping at generation %d", gen)
                break
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

    def _reproduce(self) -> List[Individual]:
        n_children = self.cfg.population_size * self.cfg.amplification
        mute_rate = self.operators.adaptive_mute_rate(self.population)
        children: List[Individual] = []
        while len(children) < n_children:
            p1 = self.rng.choice(self.population)
            p2 = self.rng.choice(self.population)
            child = self.operators.single_point_crossover(p1, p2, rng=self.rng)
            child = self.operators.mutate(child, mute_rate=mute_rate, not_mutate_idx=self.cfg.not_mutate_idx, rng=self.rng)
            children.append(child)
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
        results = self._pool.map(_eval_chunk, chunks)
        flat: Dict[str, dict] = {}
        for r in results:
            flat.update(r)
        return [flat[seq] for seq in seqs]

    def _check_convergence(self) -> bool:
        if self.cfg.early_stop_patience <= 0:
            return False
        fronts = fast_non_dominated_sort(self.fitness_list, self.objectives)
        if not fronts or not fronts[0]:
            return False
        front_fitness = [self.fitness_list[i] for i in fronts[0]]
        summary = {obj.key: sum(f.get(obj.key, 0.0) for f in front_fitness) / len(front_fitness) for obj in self.objectives}
        self._best_front_history.append(summary)
        if len(self._best_front_history) < 2:
            return False
        prev, curr = self._best_front_history[-2], self._best_front_history[-1]
        improved = any(abs(curr[k] - prev[k]) > self.cfg.fitness_config.cai_tolerance for k in curr)
        self._stagnation_count = 0 if improved else self._stagnation_count + 1
        return self._stagnation_count >= self.cfg.early_stop_patience

    def _summary_dict(self) -> dict:
        if not self.fitness_list:
            return {}
        best = self.fitness_list[0]
        return {k: f"{best.get(k, 0):.3f}" for k in ("CAI", "avg_MFE", "CG_content", "AUP")}

    def _write_logs(self, log_file: Path) -> None:
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"Gen {self.generation}: {self._summary_dict()}\n")

    def _save_checkpoint(self, path: Path) -> None:
        data = {
            "version": CHECKPOINT_VERSION,
            "protein_sequence": self.spec.protein_sequence,
            "generation": self.generation,
            "population": [ind.indices for ind in self.population],
            "fitness_list": self.fitness_list,
            "config": self.cfg,
            "rng_state": self.rng.getstate(),
            "best_front_history": self._best_front_history,
            "stagnation_count": self._stagnation_count,
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
        instance._best_front_history = data.get("best_front_history", [])
        instance._stagnation_count = data.get("stagnation_count", 0)
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


def _eval_chunk(seqs: List[str]) -> Dict[str, dict]:
    from cdsopt.fitness.evaluator import FitnessEvaluator, FitnessConfig
    ev = FitnessEvaluator(config=FitnessConfig())
    return {seq: ev.evaluate(seq) for seq in seqs}
