# -*- coding: utf-8 -*-
"""NSGA-II: fast non-dominated sort + crowding distance."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True)
class Objective:
    key: str
    target: float | None = None
    tolerance: float = 0.0


def _obj_value(fitness: dict, obj: Objective) -> float:
    raw = float(fitness.get(obj.key, 0.0))
    if obj.target is None:
        return raw
    dist = abs(raw - obj.target)
    return 0.0 if dist <= obj.tolerance else dist


def dominates(a: dict, b: dict, objectives: List[Objective]) -> bool:
    better = False
    for obj in objectives:
        va, vb = _obj_value(a, obj), _obj_value(b, obj)
        if va > vb:
            return False
        if va < vb:
            better = True
    return better


def fast_non_dominated_sort(fitness_list: List[dict], objectives: List[Objective]) -> List[List[int]]:
    n = len(fitness_list)
    domination_count = [0] * n
    dominated_solutions: List[List[int]] = [[] for _ in range(n)]
    fronts: List[List[int]] = [[]]

    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            if dominates(fitness_list[p], fitness_list[q], objectives):
                dominated_solutions[p].append(q)
            elif dominates(fitness_list[q], fitness_list[p], objectives):
                domination_count[p] += 1
        if domination_count[p] == 0:
            fronts[0].append(p)

    i = 0
    while fronts[i]:
        next_front: List[int] = []
        for p in fronts[i]:
            for q in dominated_solutions[p]:
                domination_count[q] -= 1
                if domination_count[q] == 0:
                    next_front.append(q)
        i += 1
        fronts.append(next_front)

    if fronts and not fronts[-1]:
        fronts.pop()
    return fronts


def crowding_distance(front: List[int], fitness_list: List[dict], objectives: List[Objective]) -> Dict[int, float]:
    if len(front) <= 2:
        return {idx: float("inf") for idx in front}

    distances: Dict[int, float] = {idx: 0.0 for idx in front}
    for obj in objectives:
        sorted_front = sorted(front, key=lambda i: _obj_value(fitness_list[i], obj))
        f_min = _obj_value(fitness_list[sorted_front[0]], obj)
        f_max = _obj_value(fitness_list[sorted_front[-1]], obj)
        distances[sorted_front[0]] = float("inf")
        distances[sorted_front[-1]] = float("inf")
        if f_max == f_min:
            continue
        for j in range(1, len(sorted_front) - 1):
            prev_val = _obj_value(fitness_list[sorted_front[j - 1]], obj)
            next_val = _obj_value(fitness_list[sorted_front[j + 1]], obj)
            distances[sorted_front[j]] += (next_val - prev_val) / (f_max - f_min)
    return distances


def _count_satisfied(fitness: dict, objectives: List[Objective]) -> int:
    """Count how many objectives are within target ± tolerance."""
    count = 0
    for obj in objectives:
        if obj.target is None:
            continue
        raw = float(fitness.get(obj.key, 0.0))
        if abs(raw - obj.target) <= obj.tolerance:
            count += 1
    return count


def environmental_selection(fitness_list: List[dict], n_select: int, objectives: List[Objective]) -> List[int]:
    if n_select >= len(fitness_list):
        return list(range(len(fitness_list)))

    fronts = fast_non_dominated_sort(fitness_list, objectives)
    selected: List[int] = []

    def _sort_key(i):
        fit = fitness_list[i]
        satisfied = _count_satisfied(fit, objectives)
        cai_dist = 0.0
        avg_mfe_dist = 0.0
        for obj in objectives:
            if obj.key == "CAI":
                cai_dist = abs(fit.get(obj.key, 0.0) - obj.target)
            elif obj.key == "avg_MFE":
                avg_mfe_dist = abs(fit.get(obj.key, 0.0) - obj.target)
        return (-satisfied, cai_dist, avg_mfe_dist)

    for front in fronts:
        front_sorted = sorted(front, key=_sort_key)
        if len(selected) + len(front_sorted) <= n_select:
            selected.extend(front_sorted)
        else:
            selected.extend(front_sorted[: n_select - len(selected)])
            break
    return selected
