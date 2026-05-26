# -*- coding: utf-8 -*-
"""Genetic operators: crossover, mutation, adaptive mutation rate, elitism."""
from __future__ import annotations

import random
from typing import List

from cdsopt.genetic_alg.individual import Individual, ProteinSpec
from cdsopt.genetic_alg.nsga2 import Objective, fast_non_dominated_sort


class GeneticOperators:
    def __init__(self, spec: ProteinSpec, base_mute_rate: float = 0.01):
        self.spec = spec
        self.base_mute_rate = base_mute_rate

    def single_point_crossover(self, parent1: Individual, parent2: Individual, not_mutate_idx: set[int] | None = None, rng: random.Random | None = None) -> Individual:
        gen = rng or random
        length = len(parent1)
        if length <= 1:
            return parent1.copy()
        cp = gen.randint(1, length - 1)
        child_indices = parent1.indices[:cp] + parent2.indices[cp:]
        if not_mutate_idx:
            child_indices = list(child_indices)
            for i in not_mutate_idx:
                child_indices[i] = parent1.indices[i]
        return Individual(child_indices)

    def uniform_crossover(self, parent1: Individual, parent2: Individual, rng: random.Random | None = None) -> Individual:
        gen = rng or random
        return Individual([gen.choice([p1, p2]) for p1, p2 in zip(parent1.indices, parent2.indices)])

    def mutate(self, individual: Individual, mute_rate: float | None = None, not_mutate_idx: set[int] | None = None, rng: random.Random | None = None) -> Individual:
        gen = rng or random
        rate = mute_rate if mute_rate is not None else self.base_mute_rate
        indices = individual.indices.copy()
        skip = not_mutate_idx or set()
        for i in range(len(indices)):
            if i in skip or gen.random() >= rate:
                continue
            n = len(self.spec.codon_choices[i])
            if n > 1:
                new_idx = gen.randrange(n)
                while new_idx == indices[i]:
                    new_idx = gen.randrange(n)
                indices[i] = new_idx
        return Individual(indices)

    def adaptive_mute_rate(self, population: List[Individual]) -> float:
        if len(population) < 2:
            return self.base_mute_rate
        diversity = len({ind.to_tuple() for ind in population}) / len(population)
        if diversity < 0.3:
            return min(0.5, self.base_mute_rate * 3.0)
        if diversity > 0.8:
            return max(0.001, self.base_mute_rate / 2.0)
        return self.base_mute_rate

    def select_elite(self, population: List[Individual], fitness_list: List[dict], objectives: List[Objective], n_elite: int) -> List[Individual]:
        if n_elite <= 0:
            return []
        fronts = fast_non_dominated_sort(fitness_list, objectives)
        if not fronts:
            return []
        return [population[idx].copy() for idx in fronts[0][:n_elite]]
