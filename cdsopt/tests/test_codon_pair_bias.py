# -*- coding: utf-8 -*-
import pytest
from cdsopt.utils.codon_pair_bias import calc_cpb, _build_reference_sequences


class TestBuildReferenceSequences:
    def test_builds_dna_sequences(self):
        seqs = _build_reference_sequences(species="human", n_seqs=10)
        assert len(seqs) == 10
        for s in seqs:
            assert len(s) == 150  # 50 codons * 3
            assert set(s).issubset({"A", "T", "G", "C"})

    def test_reproducible_with_same_seed(self):
        s1 = _build_reference_sequences(species="human", n_seqs=5)
        s2 = _build_reference_sequences(species="human", n_seqs=5)
        assert s1 == s2


class TestCalcCpb:
    def test_returns_float(self):
        score = calc_cpb("AUGUUUAAAGGG")
        assert isinstance(score, float)

    def test_score_range(self):
        # CPB can be positive or negative depending on pair representation
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
        # Only 2 codons -> 1 pair, should still work
        score = calc_cpb("AUGUUU")
        assert isinstance(score, float)
