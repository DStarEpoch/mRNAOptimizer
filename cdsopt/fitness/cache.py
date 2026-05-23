# -*- coding: utf-8 -*-
"""
Fitness cache with LRU eviction and hit/miss statistics.

Returns *copies* of cached dicts so callers cannot accidentally mutate
internal state.
"""
from __future__ import annotations

from collections import OrderedDict
from typing import Dict


class FitnessCache:
    """Simple LRU cache for sequence -> fitness dict."""

    def __init__(self, maxsize: int = 100_000):
        self.maxsize = maxsize
        self._cache: OrderedDict[str, dict] = OrderedDict()
        self.hits: int = 0
        self.misses: int = 0

    def get(self, seq: str) -> dict | None:
        """Return a *copy* of the cached fitness on hit, None on miss."""
        if seq not in self._cache:
            self.misses += 1
            return None
        self._cache.move_to_end(seq)
        self.hits += 1
        return self._cache[seq].copy()

    def set(self, seq: str, fitness: dict) -> None:
        """Store fitness. Evicts oldest entry when at capacity."""
        if seq in self._cache:
            self._cache[seq] = fitness
            self._cache.move_to_end(seq)
            return

        self._cache[seq] = fitness
        if len(self._cache) > self.maxsize:
            self._cache.popitem(last=False)

    def clear(self) -> None:
        """Empty the cache and reset statistics."""
        self._cache.clear()
        self.hits = 0
        self.misses = 0

    @property
    def hit_rate(self) -> float:
        total = self.hits + self.misses
        return self.hits / total if total > 0 else 0.0

    def __len__(self) -> int:
        return len(self._cache)
