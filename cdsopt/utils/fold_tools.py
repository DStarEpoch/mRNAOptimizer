# -*- coding: utf-8 -*-
from ViennaRNA import RNA

def estimate_fold(sequence: str) -> dict:
    """
    Estimate the secondary structure of an RNA sequence using ViennaRNA.

    :param sequence: The RNA sequence to fold.
    :return: A dictionary containing the MFE and the folded structure.
    """
    fc = RNA.fold_compound(sequence)
    ss, mfe = fc.mfe()
    # calculate Average Unpaired Probability (AUP) from secondary structure
    aup = ss.count('.') / len(ss) if len(ss) > 0 else 0.0
    return {
        "mfe": mfe,
        "structure": ss,
        "aup": aup
    }
