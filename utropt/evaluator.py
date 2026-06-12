# -*- coding: utf-8 -*-
"""
UTR and full-mRNA evaluation module.

Combines rule-based metrics with GEMORNA deep-learning predictions:
  - 5' UTR: Kozak, uORFs, MFE, GC, CpG, GEMORNA TIE proxy
  - 3' UTR: PAS, AREs, miRNA seeds, RBP motifs, MFE, GC, GEMORNA stability proxy
  - Full mRNA: global MFE/avg_MFE/AUP, translational-start-window MFE, GC
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cdsopt.tables.genetic_code import get_code_map_by_genetic_code
from cdsopt.utils.fold_tools import estimate_fold

from utropt.gemorna_adapter import predict_3utr, predict_5utr
from utropt.mirna_seeds import miRNASeedDB
from utropt.rbp_motifs import ATtRACTParser
from utropt.utils import cpg_count, find_motifs, gc_content, to_rna

logger = logging.getLogger(__name__)


@dataclass
class UTRConfig:
    """Configuration for UTR evaluation."""

    species: str = "human"
    genetic_code: int = 1
    fold_engine: str = "auto"
    kozak_consensus: str = "GCCGCCACCATGG"
    translational_window: int = 50  # 5'UTR end + CDS first N nt
    attract_path: Path | str | None = None
    mirna_seed_db: dict[str, str] | None = None
    use_gemorna: bool = True
    mirna_allow_wobble: bool = False


class UTREvaluator:
    """Evaluate 5'/3' UTRs and full-length mRNA sequences."""

    # PAS consensus and common variants (strength-weighted)
    PAS_PATTERNS = {
        "AAUAAA": 1.0,
        "AUUAAA": 0.8,
        "AGUAAA": 0.5,
        "UAUAAA": 0.4,
        "CAUAAA": 0.4,
        "GAUAAA": 0.3,
    }

    # ARE / HuR-like motifs
    ARE_PATTERNS = ["AUUUA", "AUUUUA", "AUUUUUA", "UUUUAUU"]

    def __init__(self, config: UTRConfig | None = None):
        self.config = config or UTRConfig()
        self._code_map = get_code_map_by_genetic_code(self.config.genetic_code)
        self._stop_codons = self._get_stop_codons()
        self._valid_codons = self._get_valid_codons()

        # Optional ATtRACT parser
        if self.config.attract_path:
            self.rbp_parser = ATtRACTParser(self.config.attract_path)
        else:
            default_attract = (
                Path(__file__).resolve().parent.parent / "data" / "ATtRACT" / "ATtRACT_db.txt"
            )
            self.rbp_parser = ATtRACTParser(default_attract) if default_attract.exists() else None
            if self.rbp_parser is None:
                logger.warning("ATtRACT database not found; RBP motif scanning disabled.")

        # miRNA seed DB
        self.mirna_db = miRNASeedDB(self.config.mirna_seed_db)

    # --------------------------------------------------------------------- #
    # Helper methods
    # --------------------------------------------------------------------- #

    def _get_stop_codons(self) -> set[str]:
        stops: set[str] = set()
        for aa, codons in self._code_map.items():
            if aa == "*":
                for codon in codons:
                    stops.add(to_rna(codon))
        return stops

    def _get_valid_codons(self) -> set[str]:
        valid: set[str] = set()
        for codons in self._code_map.values():
            for codon in codons:
                valid.add(to_rna(codon))
        return valid

    def _fold(self, seq: str, need_aup: bool = True) -> dict[str, Any]:
        """Run RNA folding; return empty dict on failure."""
        try:
            return estimate_fold(seq, engine=self.config.fold_engine, need_aup=need_aup)
        except Exception as e:
            logger.warning("Folding failed for %s...: %s", seq[:30], e)
            return {"mfe": 0.0, "structure": "." * len(seq), "aup": 0.0}

    def _gemorna_5utr(self, seq: str) -> float | None:
        if not self.config.use_gemorna:
            return None
        try:
            return float(predict_5utr(seq))
        except Exception as e:
            logger.warning("GEMORNA 5'UTR prediction failed: %s", e)
            return None

    def _gemorna_3utr(self, seq: str) -> float | None:
        if not self.config.use_gemorna:
            return None
        try:
            return float(predict_3utr(seq))
        except Exception as e:
            logger.warning("GEMORNA 3'UTR prediction failed: %s", e)
            return None

    # --------------------------------------------------------------------- #
    # 5' UTR evaluation
    # --------------------------------------------------------------------- #

    def kozak_score(self, five_utr: str, cds: str) -> float:
        """Compare -6..+4 window around AUG to Kozak consensus.

        Classic consensus: GCC(A/G)CCAUGG (10 nt).
        We score against both variants and return the best match.
        """
        window = to_rna(five_utr[-6:]) + to_rna(cds[:4])
        if len(window) < 10:
            window += "N" * (10 - len(window))

        consensuses = ["GCCGCCATGG", "GCCACCATGG"]
        best = 0.0
        for consensus in consensuses:
            matches = sum(1 for a, b in zip(window, consensus) if a == b)
            best = max(best, matches / len(consensus))
        return best

    def find_uorfs(self, seq_5utr: str) -> list[dict]:
        """Find AUG-initiated upstream ORFs in 5'UTR."""
        seq = to_rna(seq_5utr)
        uorfs = []
        for i in range(len(seq) - 2):
            if seq[i : i + 3] == "AUG":
                for j in range(i + 3, len(seq) - 2, 3):
                    codon = seq[j : j + 3]
                    if codon in self._stop_codons:
                        uorfs.append({"start": i, "end": j + 3, "length": j + 3 - i})
                        break
        return uorfs

    def evaluate_5utr(self, seq_5utr: str, cds_start: str = "") -> dict[str, Any]:
        """Evaluate a 5'UTR sequence."""
        seq = to_rna(seq_5utr)
        fold = self._fold(seq, need_aup=True)
        uorfs = self.find_uorfs(seq)

        result = {
            "length": len(seq),
            "gc_content": gc_content(seq),
            "cpg_count": cpg_count(seq),
            "kozak_score": self.kozak_score(seq, cds_start),
            "uorf_count": len(uorfs),
            "uorfs": uorfs,
            "mfe": float(fold.get("mfe", 0.0)),
            "avg_mfe": float(fold.get("mfe", 0.0)) / len(seq) if seq else 0.0,
            "aup": float(fold.get("aup", 0.0)),
        }

        gemorna = self._gemorna_5utr(seq)
        if gemorna is not None:
            result["gemorna_tie"] = gemorna

        return result

    # --------------------------------------------------------------------- #
    # 3' UTR evaluation
    # --------------------------------------------------------------------- #

    def evaluate_pas(self, seq_3utr: str) -> dict[str, Any]:
        """Evaluate polyadenylation signals."""
        seq = to_rna(seq_3utr)
        hits = []
        for variant, strength in self.PAS_PATTERNS.items():
            start = 0
            while True:
                idx = seq.find(variant, start)
                if idx == -1:
                    break
                hits.append(
                    {
                        "position": idx,
                        "signal": variant,
                        "strength": strength,
                        "distance_to_end": len(seq) - idx,
                    }
                )
                start = idx + 1
        return {
            "pas_count": len(hits),
            "pas_hits": hits,
            "has_classical_pas": any(h["signal"] == "AAUAAA" for h in hits),
            "max_pas_strength": max((h["strength"] for h in hits), default=0.0),
        }

    def evaluate_ares(self, seq_3utr: str) -> dict[str, Any]:
        """Evaluate AU-rich elements."""
        seq = to_rna(seq_3utr)
        hits = find_motifs(seq, self.ARE_PATTERNS)
        return {
            "are_count": len(hits),
            "are_hits": hits,
        }

    def evaluate_3utr(self, seq_3utr: str) -> dict[str, Any]:
        """Evaluate a 3'UTR sequence."""
        seq = to_rna(seq_3utr)
        fold = self._fold(seq, need_aup=True)
        pas = self.evaluate_pas(seq)
        ares = self.evaluate_ares(seq)
        mirna = self.mirna_db.summary(seq, allow_wobble=self.config.mirna_allow_wobble)

        result = {
            "length": len(seq),
            "gc_content": gc_content(seq),
            "mfe": float(fold.get("mfe", 0.0)),
            "avg_mfe": float(fold.get("mfe", 0.0)) / len(seq) if seq else 0.0,
            "aup": float(fold.get("aup", 0.0)),
            **pas,
            **ares,
            "mirna_seed_hits": mirna["total_hits"],
            "mirna_seed_mirnas": mirna["unique_mirnas"],
        }

        if self.rbp_parser is not None:
            try:
                rbp_summary = self.rbp_parser.summary(seq, organism="Homo_sapiens")
                result["rbp_total_hits"] = rbp_summary["total_hits"]
                result["rbp_stabilizing_hits"] = rbp_summary["stabilizing_hits"]
                result["rbp_destabilizing_hits"] = rbp_summary["destabilizing_hits"]
                result["rbp_unique"] = rbp_summary["unique_rbps"]
            except Exception as e:
                logger.warning("RBP motif scanning failed: %s", e)

        gemorna = self._gemorna_3utr(seq)
        if gemorna is not None:
            result["gemorna_stability"] = gemorna

        return result

    # --------------------------------------------------------------------- #
    # Full mRNA evaluation
    # --------------------------------------------------------------------- #

    def evaluate_mrna(
        self,
        five_utr: str,
        cds: str,
        three_utr: str,
        poly_a: str = "",
    ) -> dict[str, Any]:
        """Evaluate a full-length mRNA assembled from UTR + CDS + polyA."""
        five = to_rna(five_utr)
        cds_rna = to_rna(cds)
        three = to_rna(three_utr)
        poly = to_rna(poly_a)

        full_seq = five + cds_rna + three + poly
        full_fold = self._fold(full_seq, need_aup=True)

        # Translational start window: last part of 5'UTR + first N nt of CDS
        win_len = self.config.translational_window
        five_part = five[-win_len:] if len(five) >= win_len else five
        cds_part = cds_rna[: win_len - len(five_part)]
        start_window = five_part + cds_part
        start_fold = self._fold(start_window, need_aup=True)

        result = {
            "total_length": len(full_seq),
            "five_utr_length": len(five),
            "cds_length": len(cds_rna),
            "three_utr_length": len(three),
            "polya_length": len(poly),
            "gc_content": gc_content(full_seq),
            "mfe_full": float(full_fold.get("mfe", 0.0)),
            "avg_mfe_full": float(full_fold.get("mfe", 0.0)) / len(full_seq) if full_seq else 0.0,
            "aup_full": float(full_fold.get("aup", 0.0)),
            "mfe_start_window": float(start_fold.get("mfe", 0.0)),
            "avg_mfe_start_window": float(start_fold.get("mfe", 0.0)) / len(start_window)
            if start_window
            else 0.0,
            "aup_start_window": float(start_fold.get("aup", 0.0)),
            "five_utr": self.evaluate_5utr(five, cds_rna),
            "three_utr": self.evaluate_3utr(three),
        }

        return result
