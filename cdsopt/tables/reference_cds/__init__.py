# -*- coding: utf-8 -*-
"""Reference CDS sequences for CPB calculation."""
from __future__ import annotations

import os
from typing import List


def get_available_species() -> List[str]:
    """Return list of species with reference CDS files."""
    fa_dir = os.path.dirname(__file__)
    species = []
    for filename in os.listdir(fa_dir):
        if filename.endswith(".fa"):
            species.append(filename[:-3].replace("_", " "))
    return species


def get_reference_sequences(species: str = "human") -> List[str]:
    """Load reference DNA sequences for a given species.

    :param species: Species name (default "human").
    :return: List of DNA sequences (T-based, no headers).
    :raises ValueError: If species is not available.
    :raises FileNotFoundError: If the reference file is missing.
    """
    species = species.lower().replace(" ", "_")
    fa_path = os.path.join(os.path.dirname(__file__), f"{species}.fa")

    available = get_available_species()
    if species.replace("_", " ") not in available:
        raise ValueError(
            f"Reference CDS for species '{species}' is not available. "
            f"Available species: {available}"
        )

    if not os.path.exists(fa_path):
        raise FileNotFoundError(f"Reference CDS file not found: {fa_path}")

    sequences: List[str] = []
    with open(fa_path, "r", encoding="utf-8") as f:
        lines = f.read().strip().splitlines()

    current_seq: List[str] = []
    for line in lines:
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq).upper().replace("U", "T"))
                current_seq = []
        else:
            current_seq.append(line.strip())
    if current_seq:
        sequences.append("".join(current_seq).upper().replace("U", "T"))

    return sequences
