# -*- coding: utf-8 -*-
import pytest
from cdsopt.genetic_alg.nsga2 import (
    Objective,
    dominates,
    fast_non_dominated_sort,
    crowding_distance,
    environmental_selection,
)


class TestDominates:
    def test_a_dominates_b(self):
        objs = [Objective("f1"), Objective("f2")]
        a = {"f1": 1.0, "f2": 2.0}
        b = {"f1": 2.0, "f2": 3.0}
        assert dominates(a, b, objs) is True
        assert dominates(b, a, objs) is False

    def test_no_domination_tradeoff(self):
        objs = [Objective("f1"), Objective("f2")]
        a = {"f1": 1.0, "f2": 3.0}
        b = {"f1": 2.0, "f2": 2.0}
        assert dominates(a, b, objs) is False
        assert dominates(b, a, objs) is False

    def test_target_distance(self):
        objs = [Objective("f1", target=0.0)]
        a = {"f1": 1.0}
        b = {"f1": 2.0}
        assert dominates(a, b, objs) is True
        assert dominates(b, a, objs) is False

    def test_equal_values_no_domination(self):
        objs = [Objective("f1")]
        a = {"f1": 1.0}
        b = {"f1": 1.0}
        assert dominates(a, b, objs) is False
        assert dominates(b, a, objs) is False

    def test_tolerance_considers_equal(self):
        objs = [Objective("f1", target=1.0, tolerance=0.1)]
        a = {"f1": 1.05}
        b = {"f1": 1.08}
        # Both within 1.0 ± 0.1, so _obj_value returns 0.0 for both
        assert dominates(a, b, objs) is False
        assert dominates(b, a, objs) is False

    def test_tolerance_one_outside(self):
        objs = [Objective("f1", target=1.0, tolerance=0.1)]
        a = {"f1": 1.05}   # inside
        b = {"f1": 1.50}   # outside, distance = 0.5
        assert dominates(a, b, objs) is True
        assert dominates(b, a, objs) is False


class TestFastNonDominatedSort:
    def test_three_fronts(self):
        objs = [Objective("f1"), Objective("f2")]
        fitness = [
            {"f1": 1.0, "f2": 1.0},  # rank 0 (non-dominated)
            {"f1": 2.0, "f2": 2.0},  # rank 2 (dominated by 0 and 2)
            {"f1": 1.5, "f2": 1.5},  # rank 1 (dominated by 0, dominates 1)
        ]
        fronts = fast_non_dominated_sort(fitness, objs)
        assert len(fronts) == 3
        assert 0 in fronts[0]
        assert 2 in fronts[1]
        assert 1 in fronts[2]

    def test_all_non_dominated(self):
        objs = [Objective("f1"), Objective("f2")]
        fitness = [
            {"f1": 1.0, "f2": 3.0},
            {"f1": 3.0, "f2": 1.0},
        ]
        fronts = fast_non_dominated_sort(fitness, objs)
        assert len(fronts) == 1
        assert sorted(fronts[0]) == [0, 1]


class TestCrowdingDistance:
    def test_boundary_infinite(self):
        objs = [Objective("f1")]
        fitness = [
            {"f1": 0.0},
            {"f1": 1.0},
            {"f1": 2.0},
        ]
        front = [0, 1, 2]
        cd = crowding_distance(front, fitness, objs)
        assert cd[0] == float("inf")
        assert cd[2] == float("inf")
        assert cd[1] < float("inf")

    def test_single_front_all_inf(self):
        objs = [Objective("f1")]
        fitness = [{"f1": 1.0}, {"f1": 2.0}]
        cd = crowding_distance([0, 1], fitness, objs)
        assert cd[0] == float("inf")
        assert cd[1] == float("inf")

    def test_tolerance_clamped_values_get_zero_distance(self):
        objs = [Objective("f1", target=1.0, tolerance=0.2)]
        fitness = [
            {"f1": 0.9},   # inside, _obj_value = 0.0
            {"f1": 1.1},   # inside, _obj_value = 0.0
            {"f1": 2.0},   # outside, _obj_value = 1.0
        ]
        front = [0, 1, 2]
        cd = crowding_distance(front, fitness, objs)
        # 0 and 1 both have _obj_value = 0.0, so they cluster together
        assert cd[0] == float("inf")
        assert cd[2] == float("inf")


class TestEnvironmentalSelection:
    def test_selects_correct_number(self):
        objs = [Objective("f1"), Objective("f2")]
        fitness = [
            {"f1": 1.0, "f2": 1.0},
            {"f1": 2.0, "f2": 2.0},
            {"f1": 3.0, "f2": 3.0},
        ]
        selected = environmental_selection(fitness, 2, objs)
        assert len(selected) == 2

    def test_returns_all_when_n_exceeds_pop(self):
        objs = [Objective("f1")]
        fitness = [{"f1": 1.0}, {"f1": 2.0}]
        selected = environmental_selection(fitness, 5, objs)
        assert sorted(selected) == [0, 1]
