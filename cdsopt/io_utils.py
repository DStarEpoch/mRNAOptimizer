# -*- coding: utf-8 -*-
"""Input / output helpers."""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Dict, List

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    _HAS_BIOPYTHON = True
except Exception:
    _HAS_BIOPYTHON = False

logger = logging.getLogger(__name__)


def _read_fasta_records(path: Path) -> List:
    if not _HAS_BIOPYTHON:
        raise RuntimeError("biopython is required to parse FASTA files")
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def read_protein_sequence(path: Path) -> str:
    return read_single_sequence(path, alphabet="protein")


def read_single_sequence(path: Path, alphabet: str = "rna") -> str:
    """Read a single sequence from a FASTA or plain-text file.

    :param alphabet: 'rna' keeps A/C/G/T/U; 'protein' keeps amino-acid letters.
    """
    text = path.read_text(encoding="utf-8").strip()
    if text.startswith(">"):
        seq = str(_read_fasta_records(path)[0].seq)
    else:
        seq = "".join(text.split())
    if not seq:
        raise ValueError(f"Empty sequence in {path}")
    seq = seq.upper()
    if alphabet == "rna":
        # Accept T or U; downstream code normalizes to U.
        if any(c not in "ACGTU" for c in seq):
            raise ValueError(f"Sequence in {path} contains non-RNA characters")
    return seq


def read_cds_sequences(path: Path) -> List[str]:
    """Read DNA/RNA CDS sequences from a FASTA file."""
    return [str(r.seq).upper() for r in _read_fasta_records(path)]


def write_fasta(sequences: List[str], descriptions: List[str], path: Path) -> None:
    if not _HAS_BIOPYTHON:
        with open(path, "w", encoding="utf-8") as f:
            for desc, seq in zip(descriptions, sequences):
                f.write(f">{desc}\n{seq}\n")
        return
    records = [SeqRecord(Seq(seq), id=f"seq_{i}", description=desc) for i, (seq, desc) in enumerate(zip(sequences, descriptions))]
    with open(path, "w", encoding="utf-8") as f:
        SeqIO.write(records, f, "fasta")


def write_csv(rows: List[Dict], path: Path, fieldnames: List[str] | None = None) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    headers = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(data: Dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)


def write_results(processor, output_dir: Path, prefix: str = "") -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pop, fitness, spec = processor.population, processor.fitness_list, processor.spec
    if not pop:
        logger.warning("Empty population, nothing to write")
        return

    fc = processor.cfg.fitness_config
    has_context = bool(fc.prefix or fc.suffix)

    sequences = [spec.to_rna(ind) for ind in pop]
    descriptions, rows = [], []
    full_sequences = [] if has_context else None
    for i, (seq, fit) in enumerate(zip(sequences, fitness)):
        desc = " ".join([f"gen{processor.generation}"] + [f"{k}={fit[k]:.4f}" for k in ("CAI", "avg_MFE", "CG_content", "AUP", "CPB") if k in fit])
        descriptions.append(desc)
        row = {"id": f"seq_{i}", "sequence": seq, "length": len(seq)}
        if has_context:
            full_seq = fc.prefix + seq + fc.suffix
            full_sequences.append(full_seq)
            row["full_sequence"] = full_seq
            row["full_length"] = len(full_seq)
        row.update({k: fit.get(k) for k in fit if k not in ("rna_seq", "sequence", "motif_counts")})
        rows.append(row)

    write_fasta(sequences, descriptions, output_dir / f"{prefix}pareto_front.fasta")
    if has_context and full_sequences:
        write_fasta(full_sequences, descriptions, output_dir / f"{prefix}pareto_front_full.fasta")
    write_csv(rows, output_dir / f"{prefix}fitness.csv")

    cfg = processor.cfg
    fc = cfg.fitness_config
    summary = {
        "protein_sequence": spec.protein_sequence,
        "protein_length": len(spec.protein_sequence),
        "final_generation": processor.generation,
        "population_size": len(pop),
        "config": {k: getattr(cfg, k) for k in ("population_size", "generations", "mute_rate", "n_elite", "amplification", "processes")},
        "fitness_config": {k: getattr(fc, k) for k in ("species", "genetic_code", "target_cai", "cai_tolerance", "target_avg_mfe", "avg_mfe_tolerance", "target_tai", "tai_tolerance", "target_cg_content", "cg_content_tolerance", "target_aup", "aup_tolerance", "target_cpb", "cpb_tolerance", "fold_engine", "cache_maxsize", "prefix", "suffix")},
    }
    write_summary(summary, output_dir / f"{prefix}summary.json")
