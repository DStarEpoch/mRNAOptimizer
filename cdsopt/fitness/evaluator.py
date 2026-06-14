# -*- coding: utf-8 -*-
"""Multi-objective fitness evaluator with optional fixed flanking sequences."""
from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from multiprocessing import Pool
from typing import Dict, List, Optional

from cdsopt.fitness.cache import FitnessCache
from cdsopt.genetic_alg.nsga2 import Objective
from cdsopt.utils.codon_pair_bias import calc_cpb as _calc_cpb
from cdsopt.utils.fold_tools import estimate_fold as _estimate_fold
from cdsopt.utils.scoring import calc_cai as _calc_cai, calc_tai as _calc_tai, count_cg as _count_cg

logger = logging.getLogger(__name__)


@dataclass
class FitnessConfig:
    species: str = "human"
    genetic_code: int = 1
    # CAI and avg_MFE are always enabled (target must be set)
    target_cai: float = 0.9
    cai_tolerance: float = 0.001
    target_avg_mfe: float = -0.4
    avg_mfe_tolerance: float = 0.05
    # Others are enabled only when target is not None
    target_tai: Optional[float] = None
    tai_tolerance: float = 0.001
    target_cg_content: Optional[float] = None
    cg_content_tolerance: float = 0.005
    target_aup: Optional[float] = None
    aup_tolerance: float = 0.01
    target_cpb: Optional[float] = None
    cpb_tolerance: float = 0.01
    fold_engine: str = "auto"
    cache_maxsize: int = 100_000
    # Fixed flanking sequences. When both are empty (default), the evaluator
    # folds and evaluates only the CDS. Otherwise the full molecule is
    # prefix + CDS + suffix, while codon-level metrics are still computed on
    # the CDS alone.
    prefix: str = ""
    suffix: str = ""


def _norm_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def build_objectives(cfg: FitnessConfig) -> List[Objective]:
    objs: List[Objective] = []
    objs.append(Objective("CAI", target=cfg.target_cai, tolerance=cfg.cai_tolerance))
    objs.append(Objective("avg_MFE", target=cfg.target_avg_mfe, tolerance=cfg.avg_mfe_tolerance))
    if cfg.target_tai is not None:
        objs.append(Objective("tAI", target=cfg.target_tai, tolerance=cfg.tai_tolerance))
    if cfg.target_cg_content is not None:
        objs.append(Objective("CG_content", target=cfg.target_cg_content, tolerance=cfg.cg_content_tolerance))
    if cfg.target_aup is not None:
        objs.append(Objective("AUP", target=cfg.target_aup, tolerance=cfg.aup_tolerance))
    if cfg.target_cpb is not None:
        objs.append(Objective("CPB", target=cfg.target_cpb, tolerance=cfg.cpb_tolerance))
    return objs


class FitnessEvaluator:
    def __init__(self, config: FitnessConfig | None = None):
        self.cfg = config or FitnessConfig()
        self.cache = FitnessCache(maxsize=self.cfg.cache_maxsize)
        self._resolved_fold_engine = (
            "linearfold" if shutil.which("linearfold") else "vienna"
        ) if self.cfg.fold_engine == "auto" else self.cfg.fold_engine

        from cdsopt.tables.codon_frequency_table import get_table_weights
        from cdsopt.tables.genetic_code import get_code_map_by_genetic_code
        raw_weights = get_table_weights(self.cfg.species)
        code_map = get_code_map_by_genetic_code(self.cfg.genetic_code)
        self._cai_weights: Dict[str, float] = {}
        for aa, codons in code_map.items():
            freqs = [raw_weights.get(c, 0.0) for c in codons]
            max_freq = max(freqs) if any(freqs) else 1.0
            for c, f in zip(codons, freqs):
                self._cai_weights[c] = f / max_freq if max_freq > 0 else 0.0

        self._prefix = _norm_rna(self.cfg.prefix)
        self._suffix = _norm_rna(self.cfg.suffix)
        self._has_context = bool(self._prefix or self._suffix)

        logger.debug("FitnessEvaluator ready (species=%s, engine=%s, context=%s, objectives=%s)", self.cfg.species, self._resolved_fold_engine, self._has_context, self._active_objectives())

    def _active_objectives(self) -> List[str]:
        objs = ["CAI", "avg_MFE"]
        if self.cfg.target_tai is not None:
            objs.append("tAI")
        if self.cfg.target_cg_content is not None:
            objs.append("CG_content")
        if self.cfg.target_aup is not None:
            objs.append("AUP")
        if self.cfg.target_cpb is not None:
            objs.append("CPB")
        return objs

    def _cache_key(self, cds_seq: str) -> str:
        if not self._has_context:
            return cds_seq
        return f"{cds_seq}|{self._prefix}|{self._suffix}"

    def evaluate(self, rna_seq: str) -> dict:
        """Evaluate a CDS sequence.

        If fixed flanking sequences were configured, MFE/avg_MFE/AUP/structure
        are computed on prefix + CDS + suffix, while CAI/tAI/CG/CPB remain
        CDS-only metrics.
        """
        cds_seq = _norm_rna(rna_seq)
        cache_key = self._cache_key(cds_seq)
        cached = self.cache.get(cache_key)
        if cached is not None:
            return cached

        result: dict = {"rna_seq": cds_seq}
        _run = lambda fn, key, default: self._try_evaluate(fn, key, default, result, cds_seq)

        # CDS-only metrics.
        _run(lambda: float(_calc_cai(cds_seq, weights=self._cai_weights)), "CAI", 0.0)
        if self.cfg.target_tai is not None:
            _run(lambda: float(_calc_tai(cds_seq, genetic_code=self.cfg.genetic_code, species=self.cfg.species)), "tAI", 0.0)
        if self.cfg.target_cg_content is not None:
            _run(lambda: float(_count_cg(cds_seq)), "CG_content", 0.0)
        if self.cfg.target_cpb is not None:
            _run(lambda: float(_calc_cpb(cds_seq, species=self.cfg.species, genetic_code=self.cfg.genetic_code)), "CPB", 0.0)

        # Folding: full molecule if context is set, otherwise CDS alone.
        full_seq = (self._prefix + cds_seq + self._suffix) if self._has_context else cds_seq
        need_aup = self.cfg.target_aup is not None

        def _fold():
            fd = _estimate_fold(full_seq, engine=self._resolved_fold_engine, need_aup=need_aup)
            result["MFE"] = float(fd["mfe"])
            result["avg_MFE"] = result["MFE"] / len(full_seq) if len(full_seq) > 0 else 0.0
            result["structure"] = fd["structure"]
            if need_aup:
                result["AUP"] = float(fd["aup"])

        fail_vals = {"MFE": 0.0, "avg_MFE": 0.0, "structure": "." * len(full_seq)}
        if need_aup:
            fail_vals["AUP"] = 0.0
        self._try_evaluate(_fold, "fold", None, result, cds_seq, on_fail=lambda: result.update(fail_vals))

        # Length annotations for reporting.
        result["full_length"] = len(full_seq)
        result["cds_length"] = len(cds_seq)
        result["prefix_length"] = len(self._prefix)
        result["suffix_length"] = len(self._suffix)

        self.cache.set(cache_key, result)
        return result

    def _try_evaluate(self, fn, key, default, result, rna_seq, on_fail=None):
        try:
            val = fn()
            if val is not None and key != "fold":
                result[key] = val
        except Exception as e:
            logger.warning("%s evaluation failed for %s: %s", key, rna_seq[:20], e)
            if on_fail:
                on_fail()
            elif default is not None:
                result[key] = default

    def evaluate_batch(self, rna_seqs: List[str], processes: int = 1) -> Dict[str, dict]:
        if processes <= 1:
            return {seq: self.evaluate(seq) for seq in rna_seqs}
        to_eval, results = [], {}
        for seq in rna_seqs:
            cached = self.cache.get(self._cache_key(_norm_rna(seq)))
            if cached is not None:
                results[seq] = cached
            else:
                to_eval.append(seq)
        if not to_eval:
            return results
        cfg_dict = {k: getattr(self.cfg, k) for k in (
            "species", "genetic_code",
            "target_cai", "cai_tolerance", "target_avg_mfe", "avg_mfe_tolerance",
            "target_tai", "tai_tolerance", "target_cg_content", "cg_content_tolerance",
            "target_aup", "aup_tolerance", "target_cpb", "cpb_tolerance",
            "fold_engine", "prefix", "suffix",
        )}
        with Pool(processes=processes) as pool:
            mapped = pool.map(_eval_worker, [(seq, cfg_dict) for seq in to_eval])
        for seq, fit in mapped:
            results[seq] = fit
            self.cache.set(self._cache_key(_norm_rna(seq)), fit)
        return results


def _eval_worker(args: tuple) -> tuple:
    seq, cfg_dict = args
    cfg = FitnessConfig(**cfg_dict)
    return seq, FitnessEvaluator(config=cfg).evaluate(seq)
