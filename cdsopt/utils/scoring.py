# -*- coding: utf-8 -*-
"""Scoring utilities: CAI, tAI, CG content."""
from __future__ import annotations

from copy import deepcopy

from cai2 import CAI
from codonbias.scores import TrnaAdaptationIndex

from cdsopt.tables.codon_frequency_table import get_table_weights
from cdsopt.tables.genetic_code import get_code_map_by_genetic_code
from cdsopt.tables.tgcn import get_tgcn


def calc_cai(rna_seq: str, species: str = "human", genetic_code: int = 1, weights: dict | None = None) -> float:
    if weights is None:
        weights = deepcopy(get_table_weights(species))
        code_map = get_code_map_by_genetic_code(genetic_code)
        for res in code_map:
            max_codon_freq = max(weights[codon] for codon in code_map[res] if codon in weights)
            for codon in code_map[res]:
                if codon in weights:
                    weights[codon] /= max_codon_freq
    return CAI(rna_seq.upper().replace("U", "T"), weights)


def count_cg(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("C") + sequence.count("G")) / len(sequence)


def calc_tai(rna_seq: str, genetic_code: int = 1, species: str = "human", tGCN=None, genome_id: str = None, domain: str = None):
    if tGCN is None:
        if genome_id is not None and domain is not None:
            tai = TrnaAdaptationIndex(genome_id=genome_id, domain=domain, genetic_code=genetic_code)
        else:
            tai = TrnaAdaptationIndex(tGCN=get_tgcn(species), genetic_code=genetic_code)
    else:
        tai = TrnaAdaptationIndex(tGCN=tGCN, genetic_code=genetic_code)
    return tai.get_score(rna_seq.upper().replace("U", "T"))
