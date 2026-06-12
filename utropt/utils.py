# -*- coding: utf-8 -*-
"""General RNA sequence utilities for utropt."""
from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def to_rna(seq: str) -> str:
    """Normalize DNA/RNA sequence to uppercase RNA alphabet."""
    return seq.upper().replace("T", "U")


def reverse_complement_rna(seq: str) -> str:
    """Return reverse complement of an RNA sequence."""
    comp = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}
    return "".join(comp.get(b, b) for b in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    """Fraction of G/C in sequence."""
    seq = to_rna(seq)
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def cpg_count(seq: str) -> int:
    """Count CpG dinucleotides (DNA notation: CG; RNA notation: CG)."""
    seq = to_rna(seq)
    return seq.count("CG")


def find_motifs(seq: str, patterns: list[str]) -> list[dict]:
    """Find all occurrences of given motifs in a sequence."""
    seq = to_rna(seq)
    hits = []
    for pattern in patterns:
        pattern = to_rna(pattern)
        start = 0
        while True:
            idx = seq.find(pattern, start)
            if idx == -1:
                break
            hits.append({"position": idx, "pattern": pattern})
            start = idx + 1
    return hits
