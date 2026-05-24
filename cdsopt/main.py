# -*- coding: utf-8 -*-
"""CLI entry point for cdsopt."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

import click

from cdsopt.fitness.evaluator import FitnessConfig, FitnessEvaluator
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
    "seed": None, "output_dir": "./outputs",
    "fold_engine": "auto", "amplification": 3, "mute_rate": 0.05, "n_elite": 10,
    "early_stop_patience": 50, "not_mutate_idx": "",
    "target_cai": 0.9, "tolerance_cai": 0.02,
    "target_tai": None, "tolerance_tai": 0.03,
    "target_cg": None, "tolerance_cg": 0.005,
    "target_avg_mfe": -0.4, "tolerance_avg_mfe": 0.05,
    "target_aup": None, "tolerance_aup": 0.05,
    "target_cpb": None, "tolerance_cpb": 0.01,
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
    # Resolve targets: CLI > YAML > default (None means disabled)
    def _resolve_target(key, default):
        val = p.get(key, default)
        if val is not None:
            return val
        return yaml_cfg.get(key, default)

    return FitnessConfig(
        species=p["species"],
        fold_engine=p["fold_engine"],
        target_cai=p["target_cai"],
        cai_tolerance=p["tolerance_cai"],
        target_tai=_resolve_target("target_tai", None),
        tai_tolerance=p["tolerance_tai"],
        target_avg_mfe=p["target_avg_mfe"],
        avg_mfe_tolerance=p["tolerance_avg_mfe"],
        target_aup=_resolve_target("target_aup", None),
        aup_tolerance=p["tolerance_aup"],
        target_cpb=_resolve_target("target_cpb", None),
        cpb_tolerance=p["tolerance_cpb"],
        target_cg_content=_resolve_target("target_cg", None),
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
@click.option("--fold-engine", default="auto", show_default=True, help="RNA folding engine: auto, vienna, or linearfold")
@click.option("--amplification", default=3, show_default=True, help="Offspring amplification factor (offspring = pop_size * amplification)")
@click.option("--mute-rate", default=0.05, show_default=True)
@click.option("--n-elite", default=10, show_default=True, help="Number of elite individuals preserved each generation")
@click.option("--early-stop-patience", default=50, show_default=True)
@click.option("--not-mutate-idx", default="", help="Comma-separated 0-based indices that should not mutate")
@click.option("--target-cai", type=float, default=0.9, show_default=True)
@click.option("--target-cg", type=float, default=None, show_default=True)
@click.option("--target-avg-mfe", type=float, default=-0.4, show_default=True)
@click.option("--target-aup", type=float, default=None, show_default=True)
@click.option("--target-tai", type=float, default=None, show_default=True)
@click.option("--target-cpb", type=float, default=None, show_default=True)
@click.option("--tolerance-cai", type=float, default=0.02, show_default=True)
@click.option("--tolerance-cg", type=float, default=0.005, show_default=True)
@click.option("--tolerance-avg-mfe", type=float, default=0.05, show_default=True)
@click.option("--tolerance-aup", type=float, default=0.05, show_default=True)
@click.option("--tolerance-tai", type=float, default=0.03, show_default=True)
@click.option("--tolerance-cpb", type=float, default=0.01, show_default=True)
@click.option("--init-cds", type=click.Path(exists=True, path_type=Path), default=None, help="FASTA file with initial CDS sequences")
@click.option("--weighted-init", is_flag=True, help="Generate random initial population weighted by species CAI")
def optimize(**kwargs) -> None:
    """Optimize mRNA coding sequence for a given protein.

    INPUT_PROTEIN_SEQ: Protein sequence file (FASTA or plain text).
    """
    yaml_cfg = _load_yaml(kwargs["config"]) if kwargs.get("config") else {}
    ctx = click.get_current_context()

    # Resolve parameters: CLI explicit > YAML > default
    def _resolve(key: str, default: Any) -> Any:
        source = ctx.get_parameter_source(key)
        if source == click.core.ParameterSource.COMMANDLINE:
            return kwargs[key]
        return yaml_cfg.get(key, default)

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


@click.command()
@click.argument("cds", type=click.Path(exists=True, path_type=Path))
@click.option("-s", "--species", default="human", help="Host species for codon usage", show_default=True)
@click.option("--fold-engine", default="auto", show_default=True)
@click.option("-o", "--output", type=click.Path(path_type=Path), default=None, help="Output file (CSV); default: stdout")
def report(cds: Path, species: str, fold_engine: str, output: Path | None) -> None:
    """Evaluate all fitness parameters for given CDS sequences fasta file."""
    fitness_config = FitnessConfig(
        species=species,
        fold_engine=fold_engine,
        target_cai=0.0,
        target_avg_mfe=0.0,
        target_tai=0.0,
        target_cg_content=0.0,
        target_aup=0.0,
        target_cpb=0.0,
    )
    evaluator = FitnessEvaluator(config=fitness_config)

    from Bio import SeqIO
    records = list(SeqIO.parse(cds, "fasta"))
    logger.info("Loaded %d CDS sequences from %s", len(records), cds)

    rows = []
    for record in records:
        seq = str(record.seq).upper().replace("T", "U")
        fit = evaluator.evaluate(seq)
        rows.append({
            "id": record.id,
            "length": len(seq),
            "CAI": fit.get("CAI", ""),
            "tAI": fit.get("tAI", ""),
            "CG_content": fit.get("CG_content", ""),
            "MFE": fit.get("MFE", ""),
            "avg_MFE": fit.get("avg_MFE", ""),
            "AUP": fit.get("AUP", ""),
            "CPB": fit.get("CPB", ""),
        })

    import csv, sys
    fieldnames = ["id", "length", "CAI", "tAI", "CG_content", "MFE", "avg_MFE", "AUP", "CPB"]
    fobj = open(output, "w", newline="", encoding="utf-8") if output else sys.stdout
    writer = csv.DictWriter(fobj, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    if output:
        fobj.close()
        logger.info("Report written to %s", output)


app.add_command(optimize)
app.add_command(resume)
app.add_command(report)

if __name__ == "__main__":
    app()
