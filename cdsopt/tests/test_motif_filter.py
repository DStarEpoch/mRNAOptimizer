# -*- coding: utf-8 -*-
"""Tests for cdsopt.utils.motif_filter."""

import pytest

from cdsopt.utils.motif_filter import (
    ForbiddenMotifConfig,
    RESTRICTION_ENZYMES,
    find_motifs,
    find_polyt,
    find_homopolymers,
    scan_forbidden_motifs,
    summarize_hits,
    format_motif_counts,
    reverse_complement,
)


def test_reverse_complement():
    assert reverse_complement("GCTCTTC") == "GAAGAGC"
    assert reverse_complement("GGTCTC") == "GAGACC"
    assert reverse_complement("AATAAA") == "TTTATT"


def test_find_motifs_forward_and_reverse():
    seq = "AAA" + RESTRICTION_ENZYMES["BspQ1"] + "CCC"
    hits = find_motifs(seq, {"BspQ1": RESTRICTION_ENZYMES["BspQ1"]})
    assert len(hits) == 1
    assert hits[0].name == "BspQ1"
    assert hits[0].strand == "+"

    rev_seq = "AAA" + reverse_complement(RESTRICTION_ENZYMES["BspQ1"]) + "CCC"
    hits = find_motifs(rev_seq, {"BspQ1": RESTRICTION_ENZYMES["BspQ1"]})
    assert len(hits) == 1
    assert hits[0].strand == "-"


def test_find_polyt():
    seq = "AAATTTTTTAAA"
    hits = find_polyt(seq, 6)
    assert len(hits) == 1
    assert hits[0].name == "polyT"
    assert hits[0].motif == "TTTTTT"


def test_find_homopolymers():
    seq = "AAACCCCCGTTTTT"
    hits = find_homopolymers(seq, 5)
    names = {h.name for h in hits}
    assert "poly-C" in names
    assert "poly-T" in names


def test_forbidden_motif_config_builds_motif_dict():
    cfg = ForbiddenMotifConfig(enzymes=["BspQ1", "BsaI"], motifs=["AATAAA"], polya_signals=False)
    d = cfg.build_motif_dict()
    assert "BspQ1" in d
    assert "BsaI" in d
    assert "motif_0" in d
    assert "AAUAAA" not in d  # polyA disabled


def test_scan_forbidden_motifs():
    cfg = ForbiddenMotifConfig(enzymes=["BspQ1"], polyt_min_len=6, homopolymer_min_len=6, polya_signals=False)
    seq = RESTRICTION_ENZYMES["BspQ1"] + "A" + "T" * 7
    hits = scan_forbidden_motifs(seq, cfg)
    counts = summarize_hits(hits)
    assert counts.get("BspQ1", 0) == 1
    assert counts.get("polyT", 0) == 1


def test_format_motif_counts():
    assert format_motif_counts({"BspQ1": 2, "polyT": 0}) == "BspQ1:2"
    assert format_motif_counts({"BspQ1": 1, "BsaI": 2}) == "BsaI:2|BspQ1:1"
    assert format_motif_counts({}) == ""


def test_unknown_enzyme_raises():
    with pytest.raises(ValueError):
        ForbiddenMotifConfig(enzymes=["NotARealEnzyme"])
