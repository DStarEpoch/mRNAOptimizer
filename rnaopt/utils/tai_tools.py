# -*- coding: utf-8 -*-
from codonbias.scores import TrnaAdaptationIndex
try:
    from tables.tgcn import get_tgcn
except ImportError:
    from rnaopt.tables.tgcn import get_tgcn


def calc_tai(rna_seq: str, genetic_code: int = 1, species: str = "human", tGCN=None, genome_id: str = None, domain: str = None):
    """
    Calculate the tRNA Adaptation Index (tAI) for a given RNA sequence.

    :param rna_seq: The RNA sequence to evaluate.
    :param genetic_code: The NCBI genetic code ID to use (default is 1).
    :param species: The species for which the tAI is calculated (default is "human").
                    Used to load local tGCN data when `tGCN`, `genome_id`, and `domain` are not provided.
    :param tGCN: A pandas DataFrame with columns 'anti_codon' and 'GCN'.
                 If provided, it takes precedence over `species`, `genome_id`, and `domain`.
    :param genome_id: Genome ID of the organism for fetching tGCN from GtRNAdb.
    :param domain: Taxonomic domain of the organism for fetching tGCN from GtRNAdb.
    :return: The tAI value as a float.
    :raises ValueError: If insufficient data is provided to compute tAI.
    """
    if tGCN is None:
        if genome_id is not None and domain is not None:
            try:
                tai = TrnaAdaptationIndex(genome_id=genome_id, domain=domain, genetic_code=genetic_code)
            except ImportError as e:
                raise ImportError(
                    "Fetching tGCN from GtRNAdb requires 'lxml'. "
                    "Install it with: pip install lxml"
                ) from e
        else:
            tGCN = get_tgcn(species)
            tai = TrnaAdaptationIndex(tGCN=tGCN, genetic_code=genetic_code)
    else:
        tai = TrnaAdaptationIndex(tGCN=tGCN, genetic_code=genetic_code)

    dna_seq = rna_seq.upper().replace("U", "T")
    return tai.get_score(dna_seq)
