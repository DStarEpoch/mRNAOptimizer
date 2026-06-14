# -*- coding: utf-8 -*-
"""
Build UTR combinations from GEMORNA-generated parts and recommend the best ones.

Scores each UTR independently with :class:`utropt.evaluator.UTREvaluator`
(no full-length folding, so it runs fast).  CDS optimization is intentionally
**not** handled here; use ``cdsopt`` for that.
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Callable

from utropt.evaluator import UTRConfig, UTREvaluator
from utropt.gemorna_generator import generate_3utr_candidates, generate_5utr_candidates
from utropt.utils import to_rna

logger = logging.getLogger(__name__)


@dataclass
class Candidate:
    """One assembled full-length mRNA candidate."""

    name: str
    five_utr: str
    cds: str
    three_utr: str
    poly_a: str
    scores: dict

    @property
    def full_seq(self) -> str:
        return to_rna(self.five_utr + self.cds + self.three_utr + self.poly_a)


class UTRRecommender:
    """Generate and rank mRNA candidates from UTR parts and a fixed CDS."""

    def __init__(self, config: UTRConfig | None = None):
        self.config = config or UTRConfig()
        self.evaluator = UTREvaluator(self.config)

    def build_combinations(
        self,
        tag: str,
        five_utrs: list[str],
        three_utrs: list[str],
        cds: str,
        poly_a: str = "",
    ) -> list[Candidate]:
        """Cartesian-product all UTR parts and score full-length mRNA via cdsopt fold_tools."""
        candidates: list[Candidate] = []
        for idx, (five, three) in enumerate(product(five_utrs, three_utrs), 1):
            name = f"{tag}_candidate_{idx}"
            # evaluate_mrna internally uses cdsopt.utils.fold_tools.estimate_fold.
            scores = self.evaluator.evaluate_mrna(five, cds, three, poly_a)
            candidates.append(
                Candidate(
                    name=name,
                    five_utr=five,
                    cds=cds,
                    three_utr=three,
                    poly_a=poly_a,
                    scores=scores,
                )
            )
        return candidates

    @staticmethod
    def default_score(c: Candidate) -> float:
        """Default composite objective: higher is better.

        Combines:
          - 5'UTR GEMORNA TIE proxy
          - 3'UTR GEMORNA stability proxy
          - 5'UTR MFE (less negative is better)
          - 3'UTR MFE (less negative is better)
          - penalties for uORFs / weak Kozak / missing PAS
        """
        five = c.scores.get("five_utr", {})
        three = c.scores.get("three_utr", {})

        tie = five.get("gemorna_tie", 0.0)
        stability = three.get("gemorna_stability", 0.0)
        five_mfe = five.get("mfe", 0.0)
        three_mfe = three.get("mfe", 0.0)

        uorf_penalty = 2.0 * five.get("uorf_count", 0)
        kozak = five.get("kozak_score", 0.0)
        pas_penalty = 0.0 if three.get("has_classical_pas") else 1.0

        score = 0.0
        score += 2.0 * float(tie)
        score += 2.0 * float(stability)
        score += -0.1 * abs(float(five_mfe))
        score += -0.05 * abs(float(three_mfe))
        score += 1.0 * float(kozak)
        score -= uorf_penalty
        score -= pas_penalty
        return score

    def recommend(
        self,
        candidates: list[Candidate],
        top_n: int | None = 3,
        scorer: Callable[[Candidate], float] | None = None,
    ) -> list[Candidate]:
        """Return the top-N candidates by ``scorer`` (default: :meth:`default_score`).

        If ``top_n`` is None, return all candidates sorted by score.
        """
        scorer = scorer or self.default_score
        ranked = sorted(candidates, key=scorer, reverse=True)
        return ranked if top_n is None else ranked[:top_n]

    def generate_and_recommend(
        self,
        tag: str,
        original_five_utr: str,
        original_three_utr: str,
        cds: str,
        poly_a: str = "",
        n_generated: int = 5,
        top_n: int = 3,
        utr5_lengths: str | list[str] = "short",
        utr3_lengths: str | list[str] = "long",
    ) -> dict[str, list[Candidate]]:
        """Run the GEMORNA UTR design workflow for one fixed CDS.

        Returns two recommendation groups:
          - ``baseline``: original UTRs + provided CDS
          - ``optimized_utr``: GEMORNA UTRs + provided CDS
        """
        if isinstance(utr5_lengths, str):
            utr5_lengths = [utr5_lengths]
        if isinstance(utr3_lengths, str):
            utr3_lengths = [utr3_lengths]

        logger.info("[%s] Generating GEMORNA UTR parts...", tag)
        gen_5utrs: list[str] = []
        for length in utr5_lengths:
            gen_5utrs.extend(generate_5utr_candidates(length=length, n=n_generated))
        gen_3utrs: list[str] = []
        for length in utr3_lengths:
            gen_3utrs.extend(generate_3utr_candidates(length=length, n=n_generated))

        results: dict[str, list[Candidate]] = {}

        # Baseline: original UTR + fixed CDS
        baseline = self.build_combinations(
            tag=f"{tag}_baseline",
            five_utrs=[original_five_utr],
            three_utrs=[original_three_utr],
            cds=cds,
            poly_a=poly_a,
        )
        results["baseline"] = baseline[:top_n] if (top_n is not None and len(baseline) <= top_n) else self.recommend(baseline, top_n=top_n)

        # Optimized UTR + fixed CDS
        results["optimized_utr"] = self.recommend(
            self.build_combinations(
                tag=f"{tag}_optUTR",
                five_utrs=gen_5utrs,
                three_utrs=gen_3utrs,
                cds=cds,
                poly_a=poly_a,
            ),
            top_n=top_n,
        )

        return results


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def write_candidate_fasta(candidate: Candidate, path: Path) -> None:
    """Write a candidate as a multi-record FASTA with each part separated."""
    parts = [
        ("5UTR", candidate.five_utr),
        ("CDS", candidate.cds),
        ("3UTR", candidate.three_utr),
    ]
    if candidate.poly_a:
        parts.append(("polyA", candidate.poly_a))

    lines = []
    for part_name, seq in parts:
        header = f">{candidate.name}|{part_name}|len={len(seq)}"
        lines.append(header)
        # Wrap sequence to 80 chars per line for readability.
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_recommendation_report(
    results: dict[str, list[Candidate]],
    out_dir: Path,
) -> None:
    """Write all candidates as FASTA + a summary CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for combo, candidates in results.items():
        combo_dir = out_dir / combo
        combo_dir.mkdir(parents=True, exist_ok=True)

        # Write a combined FASTA with every candidate in this combo.
        all_fasta_path = combo_dir / "all.fasta"
        all_fasta_lines: list[str] = []
        for rank, cand in enumerate(candidates, 1):
            all_fasta_lines.append(f">{cand.name}|rank={rank}")
            all_fasta_lines.append(cand.full_seq)
        all_fasta_path.write_text("\n".join(all_fasta_lines) + "\n", encoding="utf-8")

        # Write individual multi-record FASTAs for top 3, and collect CSV rows for all.
        for rank, cand in enumerate(candidates, 1):
            if rank <= 3:
                write_candidate_fasta(cand, combo_dir / f"top{rank}.fasta")
            rows.append(
                {
                    "combo": combo,
                    "rank": rank,
                    "name": cand.name,
                    "total_length": cand.scores.get("total_length"),
                    "five_utr_length": cand.scores.get("five_utr_length"),
                    "cds_length": cand.scores.get("cds_length"),
                    "three_utr_length": cand.scores.get("three_utr_length"),
                    "gc_content": cand.scores.get("gc_content"),
                    "mfe_full": cand.scores.get("mfe_full"),
                    "avg_mfe_full": cand.scores.get("avg_mfe_full"),
                    "aup_full": cand.scores.get("aup_full"),
                    "mfe_start_window": cand.scores.get("mfe_start_window"),
                    "avg_mfe_start_window": cand.scores.get("avg_mfe_start_window"),
                    "five_utr_mfe": cand.scores.get("five_utr", {}).get("mfe"),
                    "three_utr_mfe": cand.scores.get("three_utr", {}).get("mfe"),
                    "gemorna_tie": cand.scores.get("five_utr", {}).get("gemorna_tie"),
                    "gemorna_stability": cand.scores.get("three_utr", {}).get("gemorna_stability"),
                    "kozak_score": cand.scores.get("five_utr", {}).get("kozak_score"),
                    "uorf_count": cand.scores.get("five_utr", {}).get("uorf_count"),
                    "has_classical_pas": cand.scores.get("three_utr", {}).get("has_classical_pas"),
                }
            )

    fieldnames = [
        "combo",
        "rank",
        "name",
        "total_length",
        "five_utr_length",
        "cds_length",
        "three_utr_length",
        "gc_content",
        "mfe_full",
        "avg_mfe_full",
        "aup_full",
        "mfe_start_window",
        "avg_mfe_start_window",
        "five_utr_mfe",
        "three_utr_mfe",
        "gemorna_tie",
        "gemorna_stability",
        "kozak_score",
        "uorf_count",
        "has_classical_pas",
    ]
    with open(out_dir / "recommendation_report.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
