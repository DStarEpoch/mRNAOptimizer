# AGENTS.md — cdsopt / mRNAOptimizer

> This file is for coding agents. Human contributors should read `README.md` first.

## Project Overview

`cdsopt` is a Python package for optimizing mRNA coding sequences (CDS) using **multi-objective genetic algorithms** (NSGA-II). Given a protein amino-acid sequence, it searches for synonymous codon combinations that balance multiple objectives such as codon adaptation index (CAI), GC content, RNA stability (MFE / AUP), and codon pair bias (CPB).

- **Entry point**: `python -m cdsopt` (currently a stub CLI)
- **Test runner**: `python run_tests.py -v`
- **Python**: >= 3.13
- **Package layout**: `cdsopt/` is a flat package under `mRNAOptimizer/`

## Architecture

```
cdsopt/
├── __main__.py              # python -m cdsopt
├── main.py                  # CLI entry (click), currently stub
├── run_tests.py             # Unified pytest wrapper
│
├── fitness/                 # Fitness evaluation pipeline
│   ├── cache.py             # LRU cache with hit/miss stats
│   └── evaluator.py         # Multi-objective evaluator (CAI, tAI, CG, MFE, AUP, CPB)
│
├── genetic_alg/             # NSGA-II core
│   ├── individual.py        # ProteinSpec + Individual encoding
│   ├── nsga2.py             # Fast non-dominated sort, crowding distance, env. selection
│   ├── operators.py         # Crossover, synonymous mutation, adaptive mute rate, elitism
│   └── processor.py         # Main GA loop, checkpointing, convergence monitoring
│
├── utils/                   # Supporting utilities
│   ├── cai_tools.py         # CAI weight pre-computation
│   ├── cg_tools.py          # GC content calculation
│   ├── codon_pair_bias.py   # CPB via codonbias.scores.CodonPairBias
│   ├── fold_tools.py        # RNA folding: ViennaRNA (exact AUP) / LinearFold (fast MFE)
│   └── tai_tools.py         # tAI weight pre-computation
│
├── tables/                  # Lookup tables
│   ├── codon_frequency_table/
│   ├── genetic_code/        # Standard + vertebrate mitochondrial
│   └── tgcn/
│
└── tests/                   # 70 tests (pytest)
    ├── test_individual.py
    ├── test_fold_tools.py
    ├── test_codon_pair_bias.py
    ├── test_fitness.py
    ├── test_nsga2.py
    ├── test_operators.py
    └── test_processor.py
```

## Completed Modules

| Module | Status | Notes |
|--------|--------|-------|
| `individual.py` | ✅ | ProteinSpec pre-generates synonymous codon lists per position. Individual is an integer-index list. Immutable hash/eq for dedup. |
| `fold_tools.py` | ✅ | Auto-detects `linearfold` CLI in PATH, falls back to ViennaRNA. ViennaRNA mode uses `RNA.fold_compound` + `fc.pf()` for exact AUP (partition function). |
| `codon_pair_bias.py` | ✅ | Wraps `codonbias.scores.CodonPairBias`. Accepts RNA or DNA. Builds random synthetic reference sequences when none provided. |
| `fitness/cache.py` | ✅ | LRU cache (maxsize), returns deep copies, hit/miss rate tracking. |
| `fitness/evaluator.py` | ✅ | `FitnessConfig` dataclass with enable flags and per-objective target/tolerance. `evaluate()` and `evaluate_batch()` (multiprocess). Pre-computes CAI weights once. |
| `nsga2.py` | ✅ | Standard NSGA-II. Supports target/tolerance: `_obj_value()` returns 0.0 if inside `target ± tolerance`. |
| `operators.py` | ✅ | Single-point / uniform crossover, synonymous mutation, adaptive mute rate (diversity-driven), elite selection. |
| `processor.py` | ✅ | Main loop: elitism → reproduce (amplification) → batch evaluate → NSGA-II select. Reusable `Pool`. Checkpoint with version validation. Convergence early-stop. |
| Unit + integration tests | ✅ | 75 passed, 1 skipped. Includes CLI smoke tests. |

## Key Design Decisions

1. **Tolerance mechanism**: `_obj_value()` returns `0.0` for any objective value within `target ± tolerance`. This makes NSGA-II treat such solutions as equally optimal for that objective.
2. **`avg_MFE` not `MFE`**: The evaluator divides MFE by sequence length so the objective is length-normalized.
3. **Pool reuse**: `GeneticAlgorithmProcessor` creates a `multiprocessing.Pool` on init and reuses it across generations. Shuts down in `__del__` (guarded by `hasattr`).
4. **Circular import guard**: `genetic_alg/__init__.py` does **not** import `processor`. `processor.py` uses a local import for `build_objectives`.
5. **CAI2 mutability**: `cai2.calc_cai` mutates the weights dict passed to it. We always pass a deep copy or pre-computed read-only weights.
6. **Type annotation**: `from multiprocessing.pool import Pool` (not `multiprocessing.Pool`) to avoid IDE "Variable not allowed in type" warnings. Use `Optional[Pool]`.

## Current State (as of latest commit)

- **All core algorithm modules implemented and tested**.
- **CLI (`main.py`) is a stub**: only an empty `click.Command()` skeleton.
- **End-to-end integration tests pass**: CLI smoke tests verify `optimize`, `resume`, and YAML config paths.
- **Result serialization implemented**: `pareto_front.fasta`, `fitness.csv`, `summary.json` written after each run.
- **No visualization**: Pareto front plots, convergence curves not implemented.

## Future Tasks

### High Priority

1. ~~Implement CLI (`main.py`)~~ ✅
2. ~~Result output & serialization~~ ✅
3. ~~Integration / smoke test~~ ✅
4. **YAML config file support** ✅ (`--config` flag, CLI args override YAML values)

### Medium Priority

4. **Convergence visualization**
   - Plot Pareto front evolution (e.g., CAI vs. avg_MFE over generations)
   - Plot hypervolume or best-front summary over time
   - Optional Matplotlib / Plotly dependency

5. **Objective weighting / preference articulation**
   - Allow user to specify which objectives matter most
   - Consider weighted sum or reference-point-based selection

6. **Sequence constraints**
   - Enforce specific codons at given positions (e.g., start codon, restriction sites)
   - Avoid specific motifs (e.g., cryptic splice sites, Shine-Dalgarno-like sequences)
   - `not_mutate_idx` already exists in `GAConfig`; extend to motif avoidance

### Low Priority / Nice to Have

7. **Alternative fold engines**
   - Integrate RNAfold (ViennaRNA CLI) as another engine option
   - Support LinearFold partition function mode if available

8. **Performance**
   - Profile batch evaluation bottleneck
   - Consider caching RNA fold results across generations (many sequences may be similar)
   - Numba / Cython for hot loops if needed

9. **Documentation**
   - API docs (sphinx / mkdocs)
   - Example notebooks
   - Configuration file support (YAML/JSON) for batch runs

## Testing Conventions

- Use `pytest`.
- Each module has a corresponding `test_<module>.py`.
- Tests should avoid real `cai2` / `ViennaRNA` calls on tiny sequences (they produce NaN or unstable values). Override `fitness_list` with clean numeric dicts when testing convergence logic.
- `run_tests.py` is the unified entry point; it wraps `pytest` with logging flags.

## Coding Style

- `from __future__ import annotations` in every file.
- Type hints everywhere (`List[dict]`, `Optional[Pool]`, etc.).
- Use `dataclasses` for config objects.
- Log via `logging.getLogger(__name__)`; INFO level for setup messages.
- Docstrings: Google style or plain concise style (match existing code).
