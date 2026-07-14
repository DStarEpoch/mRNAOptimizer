# -*- coding: utf-8 -*-
"""Forbidden motif detection and repair for CDS optimization.

The module supports hard constraints: sequences containing forbidden DNA/RNA
motifs can be detected and optionally repaired by synonymous codon changes.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


#: Common restriction enzyme recognition sequences (DNA alphabet).
RESTRICTION_ENZYMES: Dict[str, str] = {
    "BspQ1": "GCTCTTC",
    "SapI": "GCTCTTC",
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "SmaI": "CCCGGG",
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "XbaI": "TCTAGA",
    "NotI": "GCGGCCGC",
    "PacI": "TTAATTAA",
    "AscI": "GGCGCGCC",
    "FseI": "GGCCGGCC",
    "PmeI": "GTTTAAAC",
    "Swal": "ATTTAAAT",
    "NdeI": "CATATG",
    "NheI": "GCTAGC",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SalI": "GTCGAC",
}

#: Eukaryotic polyadenylation signal variants (DNA alphabet).
POLYADENYLATION_SIGNALS: Dict[str, str] = {
    "AAUAAA": "AATAAA",
    "AUUAAA": "ATTAAA",
    "AAGAAA": "AAGAAA",
    "AATATA": "AATATA",
    "AATACA": "AATACA",
    "AATAAT": "AATAAT",
    "AATACT": "AATACT",
    "AACAAA": "AACAAA",
}

#: Other motifs often avoided in mRNA templates.
OTHER_MOTIFS: Dict[str, str] = {
    "ATTTA_AU_rich": "ATTTA",
    "T7_promoter_core": "TAATACGACTCACTATA",
}


def _norm_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _norm_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTU", "TGCAA")
    return seq.translate(comp)[::-1]


@dataclass
class ForbiddenMotifConfig:
    """Configuration for motifs that should be avoided in the optimized CDS.

    All motif sequences are interpreted as DNA (T, not U).  The optimizer works
    with RNA internally, but the configuration uses the more familiar DNA
    notation.
    """

    enzymes: List[str] = field(default_factory=list)
    motifs: List[str] = field(default_factory=list)
    polyt_min_len: Optional[int] = 6
    polya_signals: bool = True
    homopolymer_min_len: Optional[int] = 6
    include_other_motifs: bool = False

    def __post_init__(self):
        # Validate enzymes
        bad = [e for e in self.enzymes if e not in RESTRICTION_ENZYMES]
        if bad:
            raise ValueError(f"Unknown enzymes: {bad}. Known: {list(RESTRICTION_ENZYMES.keys())}")

    def build_motif_dict(self) -> Dict[str, str]:
        """Return a flat dictionary of named motifs to scan for."""
        motif_dict: Dict[str, str] = {}
        for enzyme in self.enzymes:
            seq = RESTRICTION_ENZYMES[enzyme]
            motif_dict[enzyme] = seq
        for i, motif in enumerate(self.motifs):
            motif_dict[f"motif_{i}"] = _norm_dna(motif)
        if self.polya_signals:
            motif_dict.update(POLYADENYLATION_SIGNALS)
        if self.include_other_motifs:
            motif_dict.update(OTHER_MOTIFS)
        return motif_dict


@dataclass
class MotifHit:
    name: str
    motif: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    strand: str  # "+" or "-"

    def involved_codon_indices(self) -> List[int]:
        """Return the 0-based codon indices that overlap with this hit."""
        first_codon = self.start // 3
        last_codon = (self.end - 1) // 3
        return list(range(first_codon, last_codon + 1))


def find_motifs(seq: str, motif_dict: Dict[str, str]) -> List[MotifHit]:
    """Find all occurrences of named motifs (forward and reverse complement)."""
    seq_dna = _norm_dna(seq)
    hits: List[MotifHit] = []
    for name, motif in motif_dict.items():
        motif_dna = _norm_dna(motif)
        rev = reverse_complement(motif_dna)
        k = len(motif_dna)
        for i in range(len(seq_dna) - k + 1):
            sub = seq_dna[i:i + k]
            if sub == motif_dna:
                hits.append(MotifHit(name, motif_dna, i, i + k, "+"))
            elif sub == rev:
                hits.append(MotifHit(name, rev, i, i + k, "-"))
    return hits


def find_polyt(seq: str, min_len: int) -> List[MotifHit]:
    """Find runs of consecutive T (RNA/DNA both accepted)."""
    if min_len is None or min_len <= 0:
        return []
    seq_dna = _norm_dna(seq)
    hits: List[MotifHit] = []
    i = 0
    while i < len(seq_dna):
        if seq_dna[i] == "T":
            j = i + 1
            while j < len(seq_dna) and seq_dna[j] == "T":
                j += 1
            length = j - i
            if length >= min_len:
                hits.append(MotifHit("polyT", "T" * length, i, j, "+"))
            i = j
        else:
            i += 1
    return hits


def find_homopolymers(seq: str, min_len: Optional[int]) -> List[MotifHit]:
    """Find runs of any single nucleotide longer than allowed."""
    if min_len is None or min_len <= 0:
        return []
    seq_dna = _norm_dna(seq)
    hits: List[MotifHit] = []
    i = 0
    while i < len(seq_dna):
        base = seq_dna[i]
        j = i + 1
        while j < len(seq_dna) and seq_dna[j] == base:
            j += 1
        length = j - i
        if length >= min_len:
            hits.append(MotifHit(f"poly-{base}", base * length, i, j, "+"))
        i = j
    return hits


def scan_forbidden_motifs(seq: str, config: ForbiddenMotifConfig) -> List[MotifHit]:
    """Scan a sequence for all forbidden motifs defined by *config*."""
    hits: List[MotifHit] = []
    motif_dict = config.build_motif_dict()
    if motif_dict:
        hits.extend(find_motifs(seq, motif_dict))
    if config.polyt_min_len:
        hits.extend(find_polyt(seq, config.polyt_min_len))
    if config.homopolymer_min_len:
        hits.extend(find_homopolymers(seq, config.homopolymer_min_len))
    hits.sort(key=lambda h: h.start)
    return hits


def summarize_hits(hits: List[MotifHit]) -> Dict[str, int]:
    """Return a dict mapping motif name to occurrence count."""
    counts: Dict[str, int] = {}
    for h in hits:
        counts[h.name] = counts.get(h.name, 0) + 1
    return counts


def format_motif_counts(counts: Dict[str, int], sep: str = "|") -> str:
    """Format non-zero motif counts as a compact string."""
    if not counts:
        return ""
    items = sorted(counts.items())
    return sep.join(f"{name}:{cnt}" for name, cnt in items if cnt > 0)


def involved_codon_indices(hits: List[MotifHit]) -> set[int]:
    """Return the union of codon indices involved in any hit."""
    idx: set[int] = set()
    for h in hits:
        idx.update(h.involved_codon_indices())
    return idx


def repair_individual(
    indices: List[int],
    spec,
    config: ForbiddenMotifConfig,
    rng: random.Random,
    max_attempts: int = 50,
) -> Optional[List[int]]:
    """Try to repair a codon index list by synonymous changes that remove forbidden motifs.

    Returns a new list of indices if successful, or None if repair failed.
    """
    from cdsopt.genetic_alg.individual import Individual

    indices = list(indices)
    for _ in range(max_attempts):
        seq = spec.to_rna(indices)
        hits = scan_forbidden_motifs(seq, config)
        if not hits:
            return indices

        # Pick a random hit and try to change one of its involved codons.
        hit = rng.choice(hits)
        codon_positions = hit.involved_codon_indices()
        rng.shuffle(codon_positions)
        repaired = False
        for pos in codon_positions:
            choices = spec.codon_choices[pos]
            if len(choices) <= 1:
                continue
            current = indices[pos]
            # Try a few alternative synonymous codons
            alts = [i for i in range(len(choices)) if i != current]
            rng.shuffle(alts)
            for alt in alts[:3]:
                new_indices = list(indices)
                new_indices[pos] = alt
                new_seq = spec.to_rna(new_indices)
                new_hits = scan_forbidden_motifs(new_seq, config)
                if len(new_hits) < len(hits):
                    indices = new_indices
                    repaired = True
                    break
            if repaired:
                break
        if not repaired:
            # No progress could be made on this hit; give up.
            return None
    return None
