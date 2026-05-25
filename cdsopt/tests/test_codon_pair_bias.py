# -*- coding: utf-8 -*-
import pytest
from cdsopt.utils.codon_pair_bias import calc_cpb
from cdsopt.tables.reference_cds import get_reference_sequences


class TestGetReferenceSequences:
    def test_returns_real_human_sequences(self):
        seqs = get_reference_sequences(species="human")
        assert len(seqs) >= 1
        for s in seqs:
            assert len(s) > 0
            assert set(s).issubset({"A", "T", "G", "C"})

    def test_reproducible(self):
        s1 = get_reference_sequences(species="human")
        s2 = get_reference_sequences(species="human")
        assert s1 == s2


class TestCalcCpb:
    def test_returns_float(self):
        score = calc_cpb("AUGUUUAAAGGG")
        assert isinstance(score, float)

    def test_score_range(self):
        score = calc_cpb("AUGUUUAAAGGG")
        assert -10.0 < score < 10.0

    def test_rna_input_accepted(self):
        score = calc_cpb("AUGUUUAAAGGG")
        assert isinstance(score, float)

    def test_dna_input_accepted(self):
        score = calc_cpb("ATGTTTAAAGGG")
        assert isinstance(score, float)

    def test_custom_ref_sequences(self):
        ref = ["ATGTTTAAAGGG", "ATGGGGAAATTT"]
        score = calc_cpb("ATGTTTAAAGGG", ref_sequences=ref)
        assert isinstance(score, float)

    def test_short_sequence(self):
        score = calc_cpb("AUGUUU")
        assert isinstance(score, float)
