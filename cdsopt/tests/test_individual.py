# -*- coding: utf-8 -*-
import pytest
from cdsopt.genetic_alg.individual import ProteinSpec, Individual


class TestProteinSpec:
    def test_from_protein_standard_code(self):
        spec = ProteinSpec.from_protein("MFK")
        assert spec.protein_sequence == "MFK"
        assert spec.length == 3
        # M -> AUG
        assert spec.codon_choices[0] == ["AUG"]
        # F -> UUU, UUC
        assert set(spec.codon_choices[1]) == {"UUU", "UUC"}
        # K -> AAA, AAG
        assert set(spec.codon_choices[2]) == {"AAA", "AAG"}

    def test_to_rna(self):
        spec = ProteinSpec.from_protein("MFK")
        ind = Individual([0, 0, 0])
        rna = spec.to_rna(ind)
        assert rna == "AUGUUUAAA"

    def test_random_individual(self):
        spec = ProteinSpec.from_protein("MFK")
        rng = __import__("random").Random(42)
        ind = spec.random_individual(rng)
        assert len(ind) == 3
        assert ind[0] == 0  # M only has 1 codon
        assert 0 <= ind[1] <= 1
        assert 0 <= ind[2] <= 1

    def test_invalid_residue_raises(self):
        with pytest.raises(ValueError):
            ProteinSpec.from_protein("MFKX")


class TestIndividual:
    def test_copy(self):
        ind = Individual([0, 1, 2])
        copy = ind.copy()
        assert copy == ind
        assert copy.indices is not ind.indices

    def test_hash_and_eq(self):
        a = Individual([0, 1, 2])
        b = Individual([0, 1, 2])
        c = Individual([0, 1, 3])
        assert a == b
        assert hash(a) == hash(b)
        assert a != c

    def test_from_codon_list(self):
        spec = ProteinSpec.from_protein("MFK")
        ind = Individual.from_codon_list(["AUG", "UUC", "AAG"], spec)
        assert ind.indices == [0, 1, 1]
        assert spec.to_rna(ind) == "AUGUUCAAG"

    def test_from_codon_list_invalid_codon_raises(self):
        spec = ProteinSpec.from_protein("MFK")
        with pytest.raises(ValueError):
            Individual.from_codon_list(["AUG", "UUU", "XXX"], spec)

    def test_from_codon_list_t_to_u_conversion(self):
        spec = ProteinSpec.from_protein("MFK")
        ind = Individual.from_codon_list(["ATG", "TTC", "AAG"], spec)
        assert ind.indices == [0, 1, 1]
