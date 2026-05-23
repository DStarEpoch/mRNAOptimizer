from .standard_code import RESIDUE_DNA_MAP as STANDARD_CODE_RESIDUE_DNA_MAP
from .standard_code import DNA_TO_RESIDUE as STANDARD_CODE_DNA_TO_RESIDUE
from .vertebrate_mitochondrial_code import RESIDUE_DNA_MAP as MITOCHONDRIAL_CODE_RESIDUE_DNA_MAP
from .vertebrate_mitochondrial_code import DNA_TO_RESIDUE as MITOCHONDRIAL_CODE_DNA_TO_RESIDUE


GENETIC_CODE_MAP = {
    1: STANDARD_CODE_RESIDUE_DNA_MAP,
    2: MITOCHONDRIAL_CODE_RESIDUE_DNA_MAP
}

DNA_TO_RESIDUE_MAP = {
    1: STANDARD_CODE_DNA_TO_RESIDUE,
    2: MITOCHONDRIAL_CODE_DNA_TO_RESIDUE
}


def get_code_map_by_genetic_code(genetic_code: int = 1) -> dict:
    code_map = GENETIC_CODE_MAP[genetic_code]
    return code_map

def get_dna_to_residue_map_by_genetic_code(genetic_code: int = 1) -> dict:
    dna_to_residue_map = DNA_TO_RESIDUE_MAP[genetic_code]
    return dna_to_residue_map
    
