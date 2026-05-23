# -*- coding: utf-8 -*-
"""
Codon Pair Bias (CPB) calculation wrapper around codonbias.

Uses codonbias.scores.CodonPairBias (Coleman et al., Science, 2008).
"""
from __future__ import annotations

from typing import List

from codonbias.scores import CodonPairBias as _CodonPairBias

from cdsopt.tables.codon_frequency_table import get_table_weights
from cdsopt.tables.genetic_code import get_code_map_by_genetic_code


def _build_reference_sequences(species: str = "human", genetic_code: int = 1, n_seqs: int = 100) -> List[str]:
    """Build synthetic reference sequences from codon frequency table."""
    import random

    weights = get_table_weights(species)
    code_map = get_code_map_by_genetic_code(genetic_code)

    # Build weighted choice lists per amino acid
    aa_codons = {}
    aa_weights = {}
    for aa, codons in code_map.items():
        aa_codons[aa] = codons
        aa_weights[aa] = [weights.get(c, 0.0) for c in codons]

    sequences = []
    rng = random.Random(42)
    for _ in range(n_seqs):
        # Random protein of 50 residues (excluding stop)
        residues = [c for c in code_map.keys() if c != "*"]
        prot = [rng.choice(residues) for _ in range(50)]
        seq = "".join(
            rng.choices(aa_codons[r], weights=aa_weights[r], k=1)[0]
            for r in prot
        )
        sequences.append(seq)
    return sequences


def calc_cpb(
    rna_seq: str,
    species: str = "human",
    genetic_code: int = 1,
    ref_sequences: List[str] | None = None,
) -> float:
    """
    Calculate Codon Pair Bias (CPB) for an RNA sequence.

    Wraps codonbias.scores.CodonPairBias (Coleman et al., Science 2008).
    The score is the mean log-ratio of observed vs expected codon-pair
    probabilities.  Positive = over-represented pairs, negative = under-represented.

    :param rna_seq: RNA sequence
    :param species: species for generating default reference sequences
    :param genetic_code: NCBI genetic code id
    :param ref_sequences: optional list of reference DNA sequences for training.
                          If None, synthetic references are built from the
                          species codon frequency table.
    :return: CPB score (float)
    """
    dna_seq = rna_seq.upper().replace("U", "T")

    if ref_sequences is None:
        ref_sequences = _build_reference_sequences(species, genetic_code)

    cpb = _CodonPairBias(ref_sequences, k_mer=2, genetic_code=genetic_code)
    return float(cpb.get_score(dna_seq))
