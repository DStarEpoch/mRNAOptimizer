# -*- coding: utf-8 -*-
import pytest
from cdsopt.fitness.cache import FitnessCache
from cdsopt.fitness.evaluator import FitnessEvaluator, FitnessConfig, build_objectives


class TestFitnessCache:
    def test_cache_hit_and_miss(self):
        cache = FitnessCache(maxsize=10)
        assert cache.get("AUG") is None
        assert cache.misses == 1
        cache.set("AUG", {"CAI": 0.9})
        assert cache.get("AUG") == {"CAI": 0.9}
        assert cache.hits == 1

    def test_cache_returns_copy(self):
        cache = FitnessCache()
        original = {"CAI": 0.9}
        cache.set("AUG", original)
        retrieved = cache.get("AUG")
        retrieved["CAI"] = 0.5
        assert cache.get("AUG") == {"CAI": 0.9}

    def test_cache_lru_eviction(self):
        cache = FitnessCache(maxsize=2)
        cache.set("A", {"v": 1})
        cache.set("B", {"v": 2})
        cache.set("C", {"v": 3})
        assert cache.get("A") is None
        assert cache.get("B") == {"v": 2}
        assert cache.get("C") == {"v": 3}

    def test_cache_clear(self):
        cache = FitnessCache()
        cache.set("AUG", {"CAI": 0.9})
        cache.clear()
        assert len(cache) == 0
        assert cache.hits == 0
        assert cache.misses == 0
        assert cache.hit_rate == 0.0

    def test_hit_rate(self):
        cache = FitnessCache()
        cache.set("A", {"v": 1})
        cache.get("A")
        cache.get("B")
        assert cache.hit_rate == 0.5


class TestBuildObjectives:
    def test_default_objectives(self):
        cfg = FitnessConfig()
        objs = build_objectives(cfg)
        keys = [o.key for o in objs]
        assert "CAI" in keys
        assert "avg_MFE" in keys
        # CG_content and AUP are disabled by default (target=None)
        assert "CG_content" not in keys
        assert "AUP" not in keys
        assert "tAI" not in keys
        assert "CPB" not in keys

    def test_tai_enabled(self):
        cfg = FitnessConfig(target_tai=0.9)
        objs = build_objectives(cfg)
        keys = [o.key for o in objs]
        assert "tAI" in keys

    def test_tolerance_passed_through(self):
        cfg = FitnessConfig(cai_tolerance=0.05)
        objs = build_objectives(cfg)
        cai_obj = next(o for o in objs if o.key == "CAI")
        assert cai_obj.tolerance == 0.05
        assert cai_obj.target == 0.9

    def test_avg_mfe_not_mfe(self):
        cfg = FitnessConfig()
        objs = build_objectives(cfg)
        keys = [o.key for o in objs]
        assert "avg_MFE" in keys
        assert "MFE" not in keys


class TestFitnessEvaluator:
    def test_default_objectives(self):
        ev = FitnessEvaluator()
        result = ev.evaluate("AUGUUUAAAGGG")
        assert "CAI" in result
        assert "avg_MFE" in result
        assert "MFE" in result
        # CG_content and AUP are disabled by default
        assert "CG_content" not in result
        assert "AUP" not in result
        assert 0.0 <= result["CAI"] <= 1.0
        assert result["avg_MFE"] == pytest.approx(result["MFE"] / len("AUGUUUAAAGGG"), abs=1e-9)

    def test_cache_used_on_second_eval(self):
        ev = FitnessEvaluator()
        r1 = ev.evaluate("AUGUUUAAAGGG")
        r2 = ev.evaluate("AUGUUUAAAGGG")
        assert r1 == r2
        assert len(ev.cache) == 1
        assert ev.cache.hits >= 1

    def test_batch_evaluation(self):
        ev = FitnessEvaluator()
        seqs = ["AUGUUUAAAGGG", "AUGGGGUUUAAA"]
        results = ev.evaluate_batch(seqs, processes=1)
        assert len(results) == 2
        for seq in seqs:
            assert seq in results
            assert "CAI" in results[seq]

    def test_tai_disabled_by_default(self):
        ev = FitnessEvaluator()
        result = ev.evaluate("AUGUUUAAAGGG")
        assert "tAI" not in result

    def test_tai_when_enabled(self):
        cfg = FitnessConfig(target_tai=0.9)
        ev = FitnessEvaluator(config=cfg)
        result = ev.evaluate("AUGUUUAAAGGG")
        assert "tAI" in result
        assert isinstance(result["tAI"], float)

    def test_cpb_when_enabled(self):
        cfg = FitnessConfig(target_cpb=0.0)
        ev = FitnessEvaluator(config=cfg)
        result = ev.evaluate("AUGUUUAAAGGG")
        assert "CPB" in result
        assert isinstance(result["CPB"], float)

    def test_invalid_sequence_handled_gracefully(self):
        ev = FitnessEvaluator()
        result = ev.evaluate("AUG")
        assert isinstance(result, dict)
        assert "CAI" in result
