# -*- coding: utf-8 -*-
import tempfile
import pytest
from cdsopt.genetic_alg.processor import GeneticAlgorithmProcessor, GAConfig
from cdsopt.genetic_alg.individual import ProteinSpec
from cdsopt.fitness.evaluator import FitnessConfig


@pytest.fixture
def short_protein():
    return "MFK"


class TestInitialization:
    def test_random_initialization(self, short_protein):
        cfg = GAConfig(population_size=10, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg)
        assert len(proc.population) == 10
        assert len(proc.fitness_list) == 10
        assert all("CAI" in f for f in proc.fitness_list)

    def test_init_cds_list(self, short_protein):
        # Provide valid RNA CDS sequences
        init = ["AUGUUUAAA", "AUGUUCAAG"]
        cfg = GAConfig(population_size=5, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg, init_cds_list=init)
        assert len(proc.population) == 5
        # First 2 should come from init_cds_list
        assert proc.spec.to_rna(proc.population[0]) in init

    def test_seed_reproducibility(self, short_protein):
        cfg = GAConfig(population_size=10, rng_seed=42, processes=1)
        p1 = GeneticAlgorithmProcessor(short_protein, config=cfg)
        p2 = GeneticAlgorithmProcessor(short_protein, config=cfg)
        assert [ind.indices for ind in p1.population] == [
            ind.indices for ind in p2.population
        ]


class TestCheckpoint:
    def test_save_and_load(self, short_protein):
        cfg = GAConfig(population_size=10, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg)
        proc.generation = 5

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            path = f.name

        proc._save_checkpoint(path)
        loaded = GeneticAlgorithmProcessor.load_checkpoint(path)

        assert loaded.generation == 5
        assert len(loaded.population) == 10
        assert loaded.cfg.population_size == 10

    def test_version_mismatch_raises(self, short_protein):
        import pickle

        cfg = GAConfig(population_size=5, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg)
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            path = f.name
        proc._save_checkpoint(path)

        # Tamper with version
        with open(path, "rb") as f:
            data = pickle.load(f)
        data["version"] = "0.0"
        with open(path, "wb") as f:
            pickle.dump(data, f)

        with pytest.raises(ValueError, match="version mismatch"):
            GeneticAlgorithmProcessor.load_checkpoint(path)


class TestStep:
    def test_step_reduces_to_population_size(self, short_protein):
        cfg = GAConfig(population_size=10, amplification=2, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg)
        proc._step()
        assert len(proc.population) == 10
        assert len(proc.fitness_list) == 10

    def test_elitism_preserves_best(self, short_protein):
        cfg = GAConfig(population_size=10, n_elite=2, amplification=1, processes=1)
        proc = GeneticAlgorithmProcessor(short_protein, config=cfg)
        best_before = [f["CAI"] for f in proc.fitness_list]
        proc._step()
        best_after = [f["CAI"] for f in proc.fitness_list]
        # Elite should keep at least some good individuals
        assert max(best_after) >= max(best_before) * 0.9
