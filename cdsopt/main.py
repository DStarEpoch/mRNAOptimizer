# -*- coding: utf-8 -*-
"""CLI entry point for cdsopt."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

import click

from cdsopt.fitness.evaluator import FitnessConfig, FitnessEvaluator
from cdsopt.genetic_alg.processor import GAConfig, GeneticAlgorithmProcessor
from cdsopt.io_utils import read_protein_sequence, read_cds_sequences, read_single_sequence, write_results
from cdsopt.utils.fold_tools import estimate_fold
from cdsopt.utils.motif_filter import ForbiddenMotifConfig, RESTRICTION_ENZYMES
from cdsopt.utils.scoring import count_cg

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
    "fold_engine": "auto", "amplification": 3, "mute_rate": 0.05, 
    "n_elite": 10,"not_mutate_idx": "",
    "target_cai": 0.9, "tolerance_cai": 0.02,
    "target_tai": None, "tolerance_tai": 0.03,
    "target_cg": None, "tolerance_cg": 0.005,
    "target_avg_mfe": -0.4, "tolerance_avg_mfe": 0.05,
    "target_aup": None, "tolerance_aup": 0.05,
    "target_cpb": None, "tolerance_cpb": 0.01,
    "prefix": "", "suffix": "",
    "forbidden_motifs": None,
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
        prefix=p.get("prefix", ""),
        suffix=p.get("suffix", ""),
        forbidden_motifs=p.get("forbidden_motifs") or ForbiddenMotifConfig(),
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
@click.option("--prefix", type=click.Path(exists=True, path_type=Path), default=None, help="FASTA/plain text with fixed 5' flanking sequence")
@click.option("--suffix", type=click.Path(exists=True, path_type=Path), default=None, help="FASTA/plain text with fixed 3' flanking sequence")
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

    # Load forbidden motif config from YAML if present.
    fm_cfg = yaml_cfg.get("forbidden_motifs")
    if isinstance(fm_cfg, dict):
        params["forbidden_motifs"] = ForbiddenMotifConfig(**fm_cfg)
    elif isinstance(fm_cfg, ForbiddenMotifConfig):
        params["forbidden_motifs"] = fm_cfg

    # Load fixed flanking sequences if provided.
    if kwargs.get("prefix"):
        params["prefix"] = read_single_sequence(kwargs["prefix"])
    if kwargs.get("suffix"):
        params["suffix"] = read_single_sequence(kwargs["suffix"])

    # not_mutate_idx can be a comma-separated string or a list in YAML.
    # Supports:
    #   - single int: 5  (fixes codon index 5)
    #   - tuple/list (start, end): (10, 20)  (fixes codons 10..20 inclusive)
    def _parse_nmi_item(item):
        if isinstance(item, int):
            return {item}
        if isinstance(item, (list, tuple)) and len(item) == 2:
            return set(range(int(item[0]), int(item[1]) + 1))
        if isinstance(item, str):
            item = item.strip()
            if not item:
                return set()
            # String tuple like "(10, 20)"
            if item.startswith("(") and item.endswith(")"):
                parts = [p.strip() for p in item[1:-1].split(",")]
                if len(parts) == 2:
                    return set(range(int(parts[0]), int(parts[1]) + 1))
            # Range like "10-20"
            if "-" in item:
                parts = item.split("-")
                if len(parts) == 2 and all(p.strip().isdigit() for p in parts):
                    return set(range(int(parts[0]), int(parts[1]) + 1))
            return {int(item)}
        return set()

    nmi = kwargs["not_mutate_idx"]
    _yaml_nmi = yaml_cfg.get("not_mutate_idx")
    not_mutate_set: set[int] = set()
    if _yaml_nmi is not None and nmi == "":
        for item in _yaml_nmi if isinstance(_yaml_nmi, list) else [_yaml_nmi]:
            not_mutate_set.update(_parse_nmi_item(item))
    else:
        for item in nmi.split(","):
            not_mutate_set.update(_parse_nmi_item(item))

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


def _looks_like_cds(rna_seq: str, genetic_code: int = 1) -> bool:
    """Heuristic check whether an RNA sequence is a CDS (or CDS fragment).

    A sequence is treated as CDS-like when:
      - its length is a multiple of 3 and >= 3 nt;
      - it only contains A/C/G/U;
      - every codon exists in the selected genetic code table;
      - no internal stop codon is present.

    Start and stop codons are intentionally *not* required so that partial
    CDS regions can also be analyzed.
    """
    seq = rna_seq.upper().replace("T", "U")
    if len(seq) % 3 != 0 or len(seq) < 3:
        return False
    if not set(seq).issubset({"A", "C", "G", "U"}):
        return False

    from cdsopt.tables.genetic_code import get_code_map_by_genetic_code
    code_map = get_code_map_by_genetic_code(genetic_code)
    valid_codons: set[str] = set()
    stop_codons: set[str] = set()
    for aa, codons in code_map.items():
        for codon in codons:
            codon = codon.upper().replace("T", "U")
            valid_codons.add(codon)
            if aa == "*":
                stop_codons.add(codon)

    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        if codon not in valid_codons:
            return False
        if codon in stop_codons:
            return False
    return True


@click.command()
@click.argument("seq", type=click.Path(exists=True, path_type=Path))
@click.option("-s", "--species", default="human", help="Host species for codon usage", show_default=True)
@click.option("--fold-engine", default="auto", show_default=True, help="RNA folding engine: auto, vienna, or linearfold")
@click.option("-o", "--output", type=click.Path(path_type=Path), default=None, help="Output file (CSV); default: stdout")
@click.option("--prefix", type=click.Path(exists=True, path_type=Path), default=None, help="Fixed 5' flanking sequence for full-length folding")
@click.option("--suffix", type=click.Path(exists=True, path_type=Path), default=None, help="Fixed 3' flanking sequence for full-length folding")
def report(seq: Path, species: str, fold_engine: str, output: Path | None,
           prefix: Path | None, suffix: Path | None) -> None:
    """Evaluate CDS sequences, optionally with fixed flanking sequences.

    SEQ: FASTA file containing RNA sequences.  Each record is treated as a CDS.
    CAI/tAI/CG/CPB are computed on the CDS, while MFE/avg_MFE/AUP are computed
    on prefix + CDS + suffix when --prefix/--suffix are provided.
    """
    prefix_seq = read_single_sequence(prefix) if prefix else ""
    suffix_seq = read_single_sequence(suffix) if suffix else ""

    # Report command scans all known motifs by default.
    report_motif_config = ForbiddenMotifConfig(
        enzymes=list(RESTRICTION_ENZYMES.keys()),
        polyt_min_len=6,
        polya_signals=True,
        homopolymer_min_len=6,
        include_other_motifs=True,
    )
    fitness_config = FitnessConfig(
        species=species,
        fold_engine=fold_engine,
        target_cai=0.0,
        target_avg_mfe=0.0,
        target_tai=0.0,
        target_cg_content=0.0,
        target_aup=0.0,
        target_cpb=0.0,
        prefix=prefix_seq,
        suffix=suffix_seq,
        forbidden_motifs=report_motif_config,
    )
    evaluator = FitnessEvaluator(config=fitness_config)

    from Bio import SeqIO
    records = list(SeqIO.parse(seq, "fasta"))
    logger.info("Loaded %d sequences from %s", len(records), seq)

    rows: list[dict] = []
    for record in records:
        rna_seq = str(record.seq).upper().replace("T", "U")
        fit = evaluator.evaluate(rna_seq)
        rows.append({
            "id": record.id,
            "full_length": fit.get("full_length", ""),
            "cds_length": fit.get("cds_length", ""),
            "prefix_length": fit.get("prefix_length", ""),
            "suffix_length": fit.get("suffix_length", ""),
            "CAI": fit.get("CAI", ""),
            "tAI": fit.get("tAI", ""),
            "CG_content": fit.get("CG_content", ""),
            "MFE": fit.get("MFE", ""),
            "avg_MFE": fit.get("avg_MFE", ""),
            "AUP": fit.get("AUP", ""),
            "CPB": fit.get("CPB", ""),
            "motif": fit.get("motif", ""),
        })

    import csv, sys
    fieldnames = [
        "id", "full_length", "cds_length", "prefix_length", "suffix_length",
        "CAI", "tAI", "CG_content", "MFE", "avg_MFE", "AUP", "CPB", "motif",
    ]
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