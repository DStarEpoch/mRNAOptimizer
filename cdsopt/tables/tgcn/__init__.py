# -*- coding: utf-8 -*-
import os
import pandas as pd


def get_tgcn(species: str = "human") -> pd.DataFrame:
    """
    Get the tRNA gene copy number (tGCN) table for a given species.

    :param species: The species name (default is "human").
    :return: A pandas DataFrame with columns 'anti_codon' and 'GCN'.
    :raises FileNotFoundError: If the tGCN table for the species is not found.
    """
    species = species.lower().replace(" ", "_")
    table_path = os.path.join(os.path.dirname(__file__), f"{species}.csv")
    if not os.path.exists(table_path):
        raise FileNotFoundError(
            f"tGCN table for species '{species}' not found at {table_path}"
        )
    return pd.read_csv(table_path)
