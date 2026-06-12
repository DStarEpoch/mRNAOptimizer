# -*- coding: utf-8 -*-
"""
RBP motif scanning using ATtRACT database.

ATtRACT file (tab-delimited) columns:
  Gene_name Gene_id Mutated Organism Motif Len Experiment_description Database
  Pubmed Experiment_description Family Matrix_id Score

Download: https://attract.cnic.es/download
Expected location: mRNAOptimizer/data/ATtRACT/ATtRACT_db.txt
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__name__)

# Simplified functional annotation for well-known RBPs.
# Effect labels are approximate; used for quick filtering only.
_RBP_EFFECT = {
    # Stabilizing / translational activation
    "ELAVL1": "stabilizing",   # HuR
    "PABPC1": "stabilizing",
    "PABPC4": "stabilizing",
    "HNRNPA2B1": "stabilizing",
    "IGF2BP1": "stabilizing",
    "IGF2BP2": "stabilizing",
    "IGF2BP3": "stabilizing",
    # Destabilizing / decay-promoting
    "ZFP36": "destabilizing",  # TTP
    "ZFP36L1": "destabilizing",
    "ZFP36L2": "destabilizing",
    "HNRNPC": "destabilizing",
    "AUF1": "destabilizing",   # HNRNPD
    "HNRNPD": "destabilizing",
    # Translational control
    "CPEB1": "translational_control",
    "CPEB2": "translational_control",
    "CPEB3": "translational_control",
    "CPEB4": "translational_control",
}


@dataclass(frozen=True)
class RBPMotif:
    gene_name: str
    gene_id: str
    organism: str
    motif: str
    experiment: str
    database: str
    family: str
    score: float

    @property
    def effect(self) -> str | None:
        return _RBP_EFFECT.get(self.gene_name.upper())


class ATtRACTParser:
    """Parse and scan ATtRACT RBP motif database."""

    DEFAULT_PATH = Path(__file__).resolve().parent.parent / "data" / "ATtRACT" / "ATtRACT_db.txt"

    def __init__(self, db_path: Path | str | None = None):
        self.db_path = Path(db_path) if db_path else self.DEFAULT_PATH
        self._motifs: list[RBPMotif] | None = None

    def load(self, organism: str | None = "Homo_sapiens") -> list[RBPMotif]:
        """Load motifs, optionally filtered by organism."""
        if self._motifs is not None:
            motifs = self._motifs
        else:
            if not self.db_path.exists():
                raise FileNotFoundError(f"ATtRACT database not found: {self.db_path}")

            motifs = []
            with self.db_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    try:
                        score = float(str(row["Score"]).replace("*", "").strip())
                    except ValueError:
                        score = 0.0
                    motifs.append(
                        RBPMotif(
                            gene_name=row["Gene_name"].strip(),
                            gene_id=row["Gene_id"].strip(),
                            organism=row["Organism"].strip(),
                            motif=row["Motif"].upper().replace("T", "U").strip(),
                            experiment=row.get("Experiment_description", "").strip(),
                            database=row.get("Database", "").strip(),
                            family=row.get("Family", "").strip(),
                            score=score,
                        )
                    )
            self._motifs = motifs

        if organism:
            return [m for m in motifs if m.organism == organism]
        return list(motifs)

    def scan(
        self,
        seq: str,
        organism: str | None = "Homo_sapiens",
        score_threshold: float = 0.0,
    ) -> list[dict]:
        """Scan a sequence for ATtRACT RBP motif hits.

        :param seq: RNA sequence (A/C/G/U/T).
        :param organism: Filter motifs by organism; None for all.
        :param score_threshold: Minimum motif score to report.
        :return: List of hit dictionaries.
        """
        seq = seq.upper().replace("T", "U")
        motifs = self.load(organism=organism)
        hits = []

        for motif in motifs:
            if motif.score < score_threshold:
                continue
            pattern = motif.motif
            start = 0
            while True:
                idx = seq.find(pattern, start)
                if idx == -1:
                    break
                hits.append(
                    {
                        "position": idx,
                        "rbp": motif.gene_name,
                        "gene_id": motif.gene_id,
                        "motif": pattern,
                        "family": motif.family,
                        "score": motif.score,
                        "effect": motif.effect,
                    }
                )
                start = idx + 1

        return hits

    def summary(
        self,
        seq: str,
        organism: str | None = "Homo_sapiens",
        score_threshold: float = 0.0,
    ) -> dict:
        """Return a summary of RBP motif hits."""
        hits = self.scan(seq, organism=organism, score_threshold=score_threshold)
        return {
            "total_hits": len(hits),
            "stabilizing_hits": len([h for h in hits if h["effect"] == "stabilizing"]),
            "destabilizing_hits": len([h for h in hits if h["effect"] == "destabilizing"]),
            "translational_control_hits": len(
                [h for h in hits if h["effect"] == "translational_control"]
            ),
            "unique_rbps": sorted({h["rbp"] for h in hits}),
        }
