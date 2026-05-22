# -*- coding: utf-8 -*-

def count_cg(sequence: str) -> float:
    """
    Calculate the CG content of an RNA sequence.

    :param sequence: The RNA sequence to analyze.
    :return: The CG content as a float.
    """
    if not sequence:
        return 0.0
    cg_count = sequence.count('C') + sequence.count('G')
    total_count = len(sequence)
    return cg_count / total_count if total_count > 0 else 0.0
