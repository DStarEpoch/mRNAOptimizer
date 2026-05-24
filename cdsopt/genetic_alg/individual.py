# -*- coding: utf-8 -*-
"""
Individual representation and protein sequence specification.

Individuals are encoded as integer lists where each position i holds an index
into the list of synonymous codons for residue i of the protein sequence.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import List, Optional

from cdsopt.tables.genetic_code import get_code_map_by_genetic_code


@dataclass(frozen=True)
class ProteinSpec:
    """Immutable specification of a protein and its synonymous codon choices."""

    protein_sequence: str
    codon_choices: List[List[str]]  # RNA codons (U-based)
    genetic_code: int = 1
    cai_weights: Optional[List[List[float]]] = None  # normalized weights per position
    total_variants: int = field(init=False)

    def __post_init__(self):
        assert len(self.protein_sequence) == len(self.codon_choices), \
            "Protein length must match codon choices length"
        object.__setattr__(self, 'total_variants', math.prod(len(c) for c in self.codon_choices))

    @classmethod
    def from_protein(cls, protein_sequence: str, genetic_code: int = 1, codon_weights: Optional[dict] = None) -> "ProteinSpec":
        protein_sequence = protein_sequence.upper()
        dna_map = get_code_map_by_genetic_code(genetic_code)
        codon_choices = []
        weight_lists = [] if codon_weights else None
        for residue in protein_sequence:
            if residue not in dna_map:
                raise ValueError(f"Residue '{residue}' not found in genetic code {genetic_code}")
            # Store as RNA codons (U-based) for consistency with folding tools
            rna_codons = [c.replace("T", "U") for c in dna_map[residue]]
            codon_choices.append(rna_codons)
            if codon_weights is not None:
                weights = [codon_weights.get(c, 1.0) for c in rna_codons]
                total = sum(weights)
                if total > 0:
                    weights = [w / total for w in weights]
                else:
                    weights = [1.0 / len(rna_codons)] * len(rna_codons)
                weight_lists.append(weights)
        return cls(protein_sequence, codon_choices, genetic_code, weight_lists)

    @property
    def length(self) -> int:
        return len(self.protein_sequence)

    def to_rna(self, individual: Individual | List[int]) -> str:
        """Convert an individual to its RNA sequence."""
        return "".join(
            self.codon_choices[i][individual[i]] for i in range(self.length)
        )

    def random_individual(self, rng: random.Random | None = None, weighted: bool = False) -> "Individual":
        """Generate a random individual. If weighted=True, sample by CAI weights."""
        gen = rng or random
        indices = []
        for i, choices in enumerate(self.codon_choices):
            if weighted and self.cai_weights and self.cai_weights[i]:
                idx = gen.choices(range(len(choices)), weights=self.cai_weights[i], k=1)[0]
            else:
                idx = gen.randrange(len(choices))
            indices.append(idx)
        return Individual(indices)


class Individual:
    """Lightweight wrapper around an integer list representing codon choices."""

    __slots__ = ("indices",)

    def __init__(self, indices: List[int]):
        self.indices = indices

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, idx: int) -> int:
        return self.indices[idx]

    def __setitem__(self, idx: int, value: int) -> None:
        self.indices[idx] = value

    def __hash__(self) -> int:
        return hash(tuple(self.indices))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Individual):
            return NotImplemented
        return self.indices == other.indices

    def copy(self) -> "Individual":
        return Individual(self.indices.copy())

    def to_tuple(self) -> tuple:
        return tuple(self.indices)

    @staticmethod
    def from_codon_list(codon_list: List[str], spec: ProteinSpec) -> "Individual":
        """Convert a list of RNA codons to an Individual."""
        indices = []
        for i, codon in enumerate(codon_list):
            codon = codon.upper().replace("T", "U")
            choices = spec.codon_choices[i]
            if codon not in choices:
                raise ValueError(
                    f"Codon '{codon}' at position {i} is not a valid synonymous codon "
                    f"for residue '{spec.protein_sequence[i]}'. Valid choices: {choices}"
                )
            indices.append(choices.index(codon))
        return Individual(indices)
