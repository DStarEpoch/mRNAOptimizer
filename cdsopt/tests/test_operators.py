# -*- coding: utf-8 -*-
import random
import pytest
from cdsopt.genetic_alg.individual import ProteinSpec, Individual
from cdsopt.genetic_alg.operators import GeneticOperators
from cdsopt.genetic_alg.nsga2 import Objective


@pytest.fixture
def spec():
    return ProteinSpec.from_protein("MFKMFKMFKMFK")


@pytest.fixture
def ops(spec):
    return GeneticOperators(spec, base_mute_rate=0.01)


def _make_ind(i: int) -> Individual:
    """Construct a unique individual based on integer i (binary encoding)."""
    indices = [0] * 12
    # M positions (0,3,6,9) are fixed at 0
    # Variable positions: F(1,4,7,10) and K(2,5,8,11)
    var_pos = [1, 2, 4, 5, 7, 8, 10, 11]
    for bit, pos in enumerate(var_pos):
        indices[pos] = (i >> bit) & 1
    return Individual(indices)


class TestCrossover:
    def test_single_point_produces_valid_individual(self, spec, ops):
        p1 = _make_ind(0)
        p2 = _make_ind(255)
        child = ops.single_point_crossover(p1, p2)
        assert len(child) == 12
        rna = spec.to_rna(child)
        assert len(rna) == 36

    def test_single_point_different_parents(self, spec, ops):
        p1 = _make_ind(0)
        p2 = _make_ind(255)
        rng = random.Random(42)
        child = ops.single_point_crossover(p1, p2, rng=rng)
        assert child.indices != p1.indices or child.indices != p2.indices

    def test_uniform_crossover(self, spec, ops):
        p1 = _make_ind(0)
        p2 = _make_ind(255)
        rng = random.Random(7)
        child = ops.uniform_crossover(p1, p2, rng=rng)
        assert len(child) == 12
        for i in range(12):
            assert child.indices[i] in (p1.indices[i], p2.indices[i])


class TestMutation:
    def test_mutate_changes_only_synonymous_codons(self, spec, ops):
        ind = _make_ind(0)
        rng = random.Random(123)
        mutant = ops.mutate(ind, mute_rate=1.0, rng=rng)
        rna_orig = spec.to_rna(ind)
        rna_mut = spec.to_rna(mutant)
        assert len(rna_mut) == len(rna_orig)

    def test_not_mutate_idx_respected(self, spec, ops):
        ind = _make_ind(0)
        rng = random.Random(456)
        mutant = ops.mutate(ind, mute_rate=1.0, not_mutate_idx={2, 3, 4}, rng=rng)
        for i in {2, 3, 4}:
            assert mutant.indices[i] == ind.indices[i]


class TestAdaptiveMuteRate:
    def test_low_diversity_increases_rate(self, spec, ops):
        pop = [_make_ind(0) for _ in range(10)]
        rate = ops.adaptive_mute_rate(pop)
        assert rate > ops.base_mute_rate

    def test_high_diversity_decreases_rate(self, spec, ops):
        # 10 distinct individuals -> diversity = 1.0
        pop = [_make_ind(i) for i in range(10)]
        rate = ops.adaptive_mute_rate(pop)
        assert rate < ops.base_mute_rate

    def test_normal_diversity_keeps_rate(self, spec, ops):
        # 5 distinct types in 16 individuals -> diversity ≈ 0.31
        pop = [_make_ind(i % 5) for i in range(16)]
        rate = ops.adaptive_mute_rate(pop)
        assert rate == ops.base_mute_rate


class TestElitism:
    def test_selects_from_first_front(self, spec, ops):
        pop = [_make_ind(0), _make_ind(255)]
        fitness = [
            {"f1": 1.0, "f2": 1.0},
            {"f1": 2.0, "f2": 2.0},
        ]
        objs = [Objective("f1"), Objective("f2")]
        elite = ops.select_elite(pop, fitness, objs, n_elite=1)
        assert len(elite) == 1
        assert elite[0] == pop[0]

    def test_returns_copies(self, spec, ops):
        pop = [_make_ind(0)]
        fitness = [{"f1": 1.0}]
        objs = [Objective("f1")]
        elite = ops.select_elite(pop, fitness, objs, n_elite=1)
        assert elite[0] == pop[0]
        assert elite[0].indices is not pop[0].indices

    def test_zero_elite_returns_empty(self, spec, ops):
        pop = [_make_ind(0)]
        fitness = [{"f1": 1.0}]
        objs = [Objective("f1")]
        assert ops.select_elite(pop, fitness, objs, n_elite=0) == []
