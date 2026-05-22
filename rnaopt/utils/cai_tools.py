# -*- coding: utf-8 -*-
from cai2 import CAI
from rnaopt.tables.codon_frequency_table import get_table_weights
from rnaopt.tables.genetic_code import get_code_map_by_genetic_code

def calc_cai(rna_seq: str, species: str = "human", genetic_code: int = 1) -> float:
    """
    Calculate the Codon Adaptation Index (CAI) for a given RNA sequence.

    :param rna_seq: The RNA sequence to evaluate.
    :param species: The species for which the CAI is calculated (default is "human").
    :param genetic_code: The genetic code to use (default is 1).
    :return: The CAI value as a float.
    """
    weight = get_table_weights(species)
    code_map = get_code_map_by_genetic_code(genetic_code)
    # calculate relative frequency of each codon to max frequency of codon in the table
    for res in code_map:
        max_codon_freq = max(weight[codon] for codon in code_map[res] if codon in weight)
        for codon in code_map[res]:
            if codon in weight:
                weight[codon] /= max_codon_freq
    dna_seq = rna_seq.upper().replace("U", "T")
    return CAI(dna_seq, weight)
