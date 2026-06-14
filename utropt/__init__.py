# -*- coding: utf-8 -*-
"""utropt: UTR evaluation and design tools for mRNAOptimizer."""
from __future__ import annotations

__version__ = "0.2.0"

from utropt.evaluator import UTRConfig, UTREvaluator
from utropt.gemorna_adapter import predict_3utr, predict_5utr
from utropt.gemorna_generator import generate_3utr_candidates, generate_5utr_candidates
from utropt.recommender import UTRRecommender
from utropt.utils import cpg_count, gc_content, to_rna

__all__ = [
    "UTRConfig",
    "UTREvaluator",
    "predict_5utr",
    "predict_3utr",
    "generate_5utr_candidates",
    "generate_3utr_candidates",
    "UTRRecommender",
    "to_rna",
    "gc_content",
    "cpg_count",
]
