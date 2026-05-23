# -*- coding: utf-8 -*-
"""Multi-objective fitness evaluator."""
from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List

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
    enable_cai: bool = True
    enable_tai: bool = False
    enable_cg: bool = True
    enable_fold: bool = True
    enable_cpb: bool = False
    target_cai: float = 0.9
    cai_tolerance: float = 0.001
    target_tai: float = 0.9
    tai_tolerance: float = 0.001
    target_avg_mfe: float = -0.4
    avg_mfe_tolerance: float = 0.05
    target_cg_content: float = 0.6
    cg_content_tolerance: float = 0.005
    target_aup: float = 0.4
    aup_tolerance: float = 0.01
    target_cpb: float = 0.5
    cpb_tolerance: float = 0.01
    fold_engine: str = "auto"
    cache_maxsize: int = 100_000


def build_objectives(cfg: FitnessConfig) -> List[Objective]:
    objs: List[Objective] = []
    if cfg.enable_cai:
        objs.append(Objective("CAI", target=cfg.target_cai, tolerance=cfg.cai_tolerance, direction="maximize"))
    if cfg.enable_tai:
        objs.append(Objective("tAI", target=cfg.target_tai, tolerance=cfg.tai_tolerance, direction="maximize"))
    if cfg.enable_cg:
        objs.append(Objective("CG_content", target=cfg.target_cg_content, tolerance=cfg.cg_content_tolerance))
    if cfg.enable_fold:
        objs.append(Objective("avg_MFE", target=cfg.target_avg_mfe, tolerance=cfg.avg_mfe_tolerance))
        objs.append(Objective("AUP", target=cfg.target_aup, tolerance=cfg.aup_tolerance))
    if cfg.enable_cpb:
        objs.append(Objective("CPB", target=cfg.target_cpb, tolerance=cfg.cpb_tolerance, direction="maximize"))
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

        logger.debug("FitnessEvaluator ready (species=%s, engine=%s, objectives=%s)", self.cfg.species, self._resolved_fold_engine, self._active_objectives())

    def _active_objectives(self) -> List[str]:
        objs = []
        if self.cfg.enable_cai: objs.append("CAI")
        if self.cfg.enable_tai: objs.append("tAI")
        if self.cfg.enable_cg: objs.append("CG")
        if self.cfg.enable_fold: objs.extend(["avg_MFE", "AUP"])
        if self.cfg.enable_cpb: objs.append("CPB")
        return objs

    def evaluate(self, rna_seq: str) -> dict:
        cached = self.cache.get(rna_seq)
        if cached is not None:
            return cached

        result: dict = {"rna_seq": rna_seq}
        _run = lambda fn, key, default: self._try_evaluate(fn, key, default, result, rna_seq)

        if self.cfg.enable_cai:
            _run(lambda: float(_calc_cai(rna_seq, weights=self._cai_weights)), "CAI", 0.0)
        if self.cfg.enable_tai:
            _run(lambda: float(_calc_tai(rna_seq, genetic_code=self.cfg.genetic_code, species=self.cfg.species)), "tAI", 0.0)
        if self.cfg.enable_cg:
            _run(lambda: float(_count_cg(rna_seq)), "CG_content", 0.0)
        if self.cfg.enable_fold:
            def _fold():
                fd = _estimate_fold(rna_seq, engine=self._resolved_fold_engine)
                result["MFE"] = float(fd["mfe"])
                result["avg_MFE"] = result["MFE"] / len(rna_seq) if len(rna_seq) > 0 else 0.0
                result["AUP"] = float(fd["aup"])
                result["structure"] = fd["structure"]
            self._try_evaluate(_fold, "fold", None, result, rna_seq, on_fail=lambda: result.update({"MFE": 0.0, "avg_MFE": 0.0, "AUP": 0.0, "structure": "." * len(rna_seq)}))
        if self.cfg.enable_cpb:
            _run(lambda: float(_calc_cpb(rna_seq, species=self.cfg.species, genetic_code=self.cfg.genetic_code)), "CPB", 0.0)

        self.cache.set(rna_seq, result)
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
            cached = self.cache.get(seq)
            if cached is not None:
                results[seq] = cached
            else:
                to_eval.append(seq)
        if not to_eval:
            return results
        cfg_dict = {k: getattr(self.cfg, k) for k in ("species", "genetic_code", "enable_cai", "enable_tai", "enable_cg", "enable_fold", "enable_cpb", "fold_engine")}
        with Pool(processes=processes) as pool:
            mapped = pool.map(_eval_worker, [(seq, cfg_dict) for seq in to_eval])
        for seq, fit in mapped:
            results[seq] = fit
            self.cache.set(seq, fit)
        return results


def _eval_worker(args: tuple) -> tuple:
    seq, cfg_dict = args
    cfg = FitnessConfig(**cfg_dict)
    return seq, FitnessEvaluator(config=cfg).evaluate(seq)
