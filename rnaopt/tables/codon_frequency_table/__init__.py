# -*- coding: utf-8 -*-
# @Time    : 2025/8/5 9:53
# @Author  : yuyeqing
# @File    : __init__.py
import os
import re
import csv


def get_available_species() -> list:
    """
    Get the list of available species for which codon frequency tables are provided.

    :return: A list of species names.
    """
    csv_dir = os.path.dirname(__file__)
    species = []
    for filename in os.listdir(csv_dir):
        if filename.endswith('.csv'):
            name = filename[:-4].replace('_', ' ')
            species.append(name)
    return species

def get_table_weights(species: str = "human") -> dict:
    """
    Get the codon frequency table weights for a given species.

    :param species: The species for which the codon frequency table is requested (default is "human").
    :return: A dictionary containing codon frequencies.
    """
    species = species.lower()
    table_file_name = "_".join(re.split(r'[ _]+', species)) + ".csv"
    table_file_path = os.path.join(os.path.dirname(__file__),  table_file_name)

    available_species = get_available_species()
    species_friendly = species.replace('_', ' ')
    if species_friendly not in available_species:
        raise ValueError(f"Species '{species}' is not available. Available species are: {available_species}")

    if not os.path.exists(table_file_path):
        raise FileNotFoundError(f"Codon frequency table for species '{species}' not found at {table_file_path}")

    weight = dict()
    with open(table_file_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        for row in reader:
            if len(row) == 3:
                codon, _, freq = row
                weight[codon] = float(freq)
    return weight
