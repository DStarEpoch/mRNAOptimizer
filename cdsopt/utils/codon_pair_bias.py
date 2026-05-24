# -*- coding: utf-8 -*-
"""
Codon Pair Bias (CPB) calculation wrapper around codonbias.

Uses codonbias.scores.CodonPairBias (Coleman et al., Science, 2008).
"""
from __future__ import annotations

from typing import List

from codonbias.scores import CodonPairBias as _CodonPairBias

from cdsopt.tables.reference_cds import get_reference_sequences


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
                          If None, real reference sequences are loaded from
                          cdsopt.tables.reference_cds.
    :return: CPB score (float)
    """
    dna_seq = rna_seq.upper().replace("U", "T")

    if ref_sequences is None:
        ref_sequences = get_reference_sequences(species)

    cpb = _CodonPairBias(ref_sequences, k_mer=2, genetic_code=genetic_code)
    return float(cpb.get_score(dna_seq))
