# -*- coding: utf-8 -*-
"""CLI entry point for cdsopt."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

import click

from cdsopt.fitness.evaluator import FitnessConfig
from cdsopt.genetic_alg.processor import GAConfig, GeneticAlgorithmProcessor
from cdsopt.io_utils import read_protein_sequence, read_cds_sequences, write_results

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def _load_yaml(path: Path) -> Dict[str, Any]:
    import yaml
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


_PARAM_DEFAULTS = {
    "species": "human", "pop_size": 100, "generations": 1000, "processes": 1,
    "seed": None, "output_dir": "./outputs", "enable_tai": False, "enable_cpb": False,
    "fold_engine": "auto", "amplification": 3, "mute_rate": 0.01, "n_elite": 2,
    "early_stop_patience": 50, "not_mutate_idx": "",
    "target_cai": 0.9, "tolerance_cai": 0.001,
    "target_tai": 0.9, "tolerance_tai": 0.001,
    "target_cg": 0.6, "tolerance_cg": 0.005,
    "target_avg_mfe": -0.4, "tolerance_avg_mfe": 0.05,
    "target_aup": 0.4, "tolerance_aup": 0.01,
    "target_cpb": 0.5, "tolerance_cpb": 0.01,
    "weighted_init": True,
}


def _validate_init_cds(path: Path | None, protein: str, genetic_code: int) -> list[str] | None:
    if path is None:
        return None
    from cdsopt.genetic_alg.individual import ProteinSpec, Individual
    spec = ProteinSpec.from_protein(protein, genetic_code)
    sequences = read_cds_sequences(path)
    valid, errors = [], []
    for seq in sequences:
        codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
        try:
            Individual.from_codon_list(codons, spec)
            valid.append(seq)
        except ValueError as e:
            errors.append(f"  - {seq[:30]}...: {e}")
    if errors:
        logger.warning("Skipped %d invalid init CDS sequences:\n%s", len(errors), "\n".join(errors))
    if not valid:
        raise ValueError(f"None of the {len(sequences)} init CDS sequences match the protein '{protein}'")
    logger.info("Using %d validated init CDS sequences", len(valid))
    return valid


def _build_fitness_config(yaml_cfg: dict, **p) -> FitnessConfig:
    return FitnessConfig(
        species=p["species"],
        fold_engine=p["fold_engine"],
        enable_tai=p["enable_tai"],
        enable_cpb=p["enable_cpb"],
        target_cai=p["target_cai"],
        cai_tolerance=p["tolerance_cai"],
        target_tai=p["target_tai"],
        tai_tolerance=p["tolerance_tai"],
        target_avg_mfe=p["target_avg_mfe"],
        avg_mfe_tolerance=p["tolerance_avg_mfe"],
        target_aup=p["target_aup"],
        aup_tolerance=p["tolerance_aup"],
        target_cpb=p["target_cpb"],
        cpb_tolerance=p["tolerance_cpb"],
        target_cg_content=p["target_cg"],
        cg_content_tolerance=p["tolerance_cg"],
    )


@click.group()
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging")
@click.pass_context
def app(ctx: click.Context, verbose: bool) -> None:
    """mRNA sequence optimizer based on multi-objective genetic algorithms."""
    _setup_logging(verbose)
    ctx.ensure_object(dict)


@click.command()
@click.argument("input_protein_seq", type=click.Path(exists=True, path_type=Path))
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path), default=None, help="YAML configuration file")
@click.option("-s", "--species", default="human", help="Host species for codon usage")
@click.option("-n", "--pop-size", default=100, show_default=True, help="Population size")
@click.option("-g", "--generations", default=1000, show_default=True, help="Max generations")
@click.option("-p", "--processes", default=1, show_default=True, help="Parallel processes")
@click.option("--seed", type=int, default=None, help="Random seed")
@click.option("-o", "--output-dir", default="./outputs", show_default=True, help="Output directory")
@click.option("--enable-tai", is_flag=True, help="Enable tAI objective")
@click.option("--enable-cpb", is_flag=True, help="Enable codon pair bias objective")
@click.option("--fold-engine", default="auto", show_default=True, help="RNA folding engine: auto, vienna, or linearfold")
@click.option("--amplification", default=3, show_default=True, help="Offspring amplification factor (offspring = pop_size * amplification)")
@click.option("--mute-rate", default=0.05, show_default=True)
@click.option("--n-elite", default=10, show_default=True, help="Number of elite individuals preserved each generation")
@click.option("--early-stop-patience", default=50, show_default=True)
@click.option("--not-mutate-idx", default="", help="Comma-separated 0-based indices that should not mutate")
@click.option("--target-cai", type=float, default=0.9, show_default=True)
@click.option("--target-cg", type=float, default=0.6, show_default=True)
@click.option("--target-avg-mfe", type=float, default=-0.4, show_default=True)
@click.option("--target-aup", type=float, default=0.3, show_default=True)
@click.option("--target-tai", type=float, default=0.9, show_default=True)
@click.option("--target-cpb", type=float, default=0.5, show_default=True)
@click.option("--tolerance-cai", type=float, default=0.01, show_default=True)
@click.option("--tolerance-cg", type=float, default=0.005, show_default=True)
@click.option("--tolerance-avg-mfe", type=float, default=0.01, show_default=True)
@click.option("--tolerance-aup", type=float, default=0.05, show_default=True)
@click.option("--tolerance-tai", type=float, default=0.01, show_default=True)
@click.option("--tolerance-cpb", type=float, default=0.01, show_default=True)
@click.option("--init-cds", type=click.Path(exists=True, path_type=Path), default=None, help="FASTA file with initial CDS sequences")
@click.option("--weighted-init", is_flag=True, help="Generate random initial population weighted by species CAI")
def optimize(**kwargs) -> None:
    """Optimize mRNA coding sequence for a given protein.

    INPUT_PROTEIN_SEQ: Protein sequence file (FASTA or plain text).
    """
    yaml_cfg = _load_yaml(kwargs["config"]) if kwargs.get("config") else {}

    # Resolve parameters: CLI > YAML > default
    def _resolve(key: str, default: Any) -> Any:
        val = kwargs[key]
        return val if val != default else yaml_cfg.get(key, default)

    protein = read_protein_sequence(kwargs["input_protein_seq"])
    logger.info("Loaded protein sequence: %d residues", len(protein))

    params = {k: _resolve(k, d) for k, d in _PARAM_DEFAULTS.items()}

    # not_mutate_idx can be a comma-separated string or a list in YAML
    nmi = kwargs["not_mutate_idx"]
    _yaml_nmi = yaml_cfg.get("not_mutate_idx")
    if _yaml_nmi is not None and nmi == "":
        nmi = ",".join(str(x) for x in _yaml_nmi) if isinstance(_yaml_nmi, list) else str(_yaml_nmi)
    not_mutate_set = {int(x.strip()) for x in nmi.split(",") if x.strip()}

    fitness_config = _build_fitness_config(yaml_cfg, **params)
    weighted_init = kwargs.get("weighted_init") or yaml_cfg.get("weighted_init", False)
    cfg = GAConfig(
        population_size=params["pop_size"],
        generations=params["generations"],
        processes=params["processes"],
        rng_seed=params["seed"],
        amplification=params["amplification"],
        mute_rate=params["mute_rate"],
        n_elite=params["n_elite"],
        early_stop_patience=params["early_stop_patience"],
        fitness_config=fitness_config,
        not_mutate_idx=not_mutate_set,
        weighted_init=weighted_init,
    )

    genetic_code = yaml_cfg.get("genetic_code", 1)
    init_cds_path = kwargs.get("init_cds") or yaml_cfg.get("init_cds")
    init_cds_list = _validate_init_cds(init_cds_path, protein, genetic_code)

    proc = GeneticAlgorithmProcessor(protein, config=cfg, init_cds_list=init_cds_list)
    proc.run(output_dir=params["output_dir"])
    write_results(proc, Path(params["output_dir"]))
    logger.info("Optimization complete. Results written to %s", params["output_dir"])


@click.command()
@click.argument("checkpoint", type=click.Path(exists=True, path_type=Path))
@click.option("-j", "--processes", type=int, default=None, help="Override number of processes")
@click.option("-o", "--output-dir", type=click.Path(path_type=Path), default=None, help="Override output directory")
def resume(checkpoint: Path, processes: int | None, output_dir: Path | None) -> None:
    """Resume optimization from a checkpoint file."""
    proc = GeneticAlgorithmProcessor.load_checkpoint(checkpoint, processes=processes)
    out = output_dir or Path("./outputs")
    logger.info("Resumed from generation %d", proc.generation)
    proc.run(output_dir=out)
    write_results(proc, out)
    logger.info("Optimization complete. Results written to %s", out)


app.add_command(optimize)
app.add_command(resume)

if __name__ == "__main__":
    app()
