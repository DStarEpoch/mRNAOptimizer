# AGENTS.md — cdsopt / mRNAOptimizer

> This file is for coding agents. Human contributors should read `README.md` first.

## Project Overview

`cdsopt` is a Python package for optimizing mRNA coding sequences (CDS) using **multi-objective genetic algorithms** (NSGA-II). Given a protein amino-acid sequence, it searches for synonymous codon combinations that balance multiple objectives such as codon adaptation index (CAI), GC content, RNA stability (MFE / AUP), and codon pair bias (CPB).

- **Entry point**: `python -m cdsopt` (click CLI with optimize / resume / report)
- **Test runner**: `python run_tests.py -v`
- **Python**: >= 3.10 (tested mainly on 3.10 for GEMORNA compatibility)
- **Package layout**: `cdsopt/` is a flat package under `mRNAOptimizer/`

## Architecture

```
cdsopt/
├── __main__.py              # python -m cdsopt
├── main.py                  # CLI entry (click): optimize / resume / report
├── run_tests.py             # Unified pytest wrapper
│
├── fitness/                 # Fitness evaluation pipeline
│   ├── cache.py             # LRU cache with hit/miss stats
│   └── evaluator.py         # Multi-objective evaluator (CAI, tAI, CG, MFE, AUP, CPB)
│
├── genetic_alg/             # NSGA-II core
│   ├── individual.py        # ProteinSpec + Individual encoding (supports weighted init)
│   ├── nsga2.py             # Fast non-dominated sort, env. selection with satisfied priority
│   ├── operators.py         # Crossover, synonymous mutation, adaptive mute rate, elitism
│   └── processor.py         # Main GA loop, checkpointing, generation fitness dump, dedup
│
├── io_utils.py              # Read protein/CDS FASTA, write results (FASTA/CSV/JSON)
│
├── utils/                   # Supporting utilities
│   ├── scoring.py           # CAI, tAI, CG content (merged from cai_tools+cg_tools+tai_tools)
│   ├── codon_pair_bias.py   # CPB via codonbias.scores.CodonPairBias
│   └── fold_tools.py        # RNA folding: ViennaRNA (exact AUP) / LinearFold (fast MFE)
│
├── tables/                  # Lookup tables and reference sequences
│   ├── codon_frequency_table/
│   ├── genetic_code/        # Standard + vertebrate mitochondrial
│   ├── tgcn/
│   └── reference_cds/       # Real reference CDS for CPB (e.g. human ACTB NCBI RefSeq)
│
└── tests/                   # pytest suite
    ├── test_cli.py
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
| `individual.py` | ✅ | ProteinSpec pre-generates synonymous codon lists. `cai_weights` + `total_variants`. Supports `weighted_init`: samples by species CAI weights (`power=10`). `random_individual(weighted=True)` for high-CAI initialization. |
| `scoring.py` | ✅ | Merged `cai_tools+cg_tools+tai_tools`. `calc_cai`, `calc_tai`, `count_cg`. |
| `fold_tools.py` | ✅ | Auto-detects `linearfold` CLI, falls back to ViennaRNA. `need_aup` flag skips `pf()/bpp()` when AUP not needed. |
| `codon_pair_bias.py` | ✅ | Wraps `codonbias.scores.CodonPairBias`. Loads **real reference CDS** from `tables/reference_cds/` (e.g. human ACTB NM_001101.5). No synthetic fallback. |
| `fitness/cache.py` | ✅ | LRU cache, hit/miss tracking. |
| `fitness/evaluator.py` | ✅ | `FitnessConfig` with per-objective `target` (Optional[float]) and `tolerance`. No `enable_*` flags: target=None means disabled. CAI and avg_MFE are always enabled. `evaluate_batch()` for multiprocessing. |
| `nsga2.py` | ✅ | `_obj_value()` returns 0.0 inside `target ± tolerance`. `environmental_selection` sorts **every front** by: 1) satisfied count desc, 2) CAI distance asc, 3) avg_MFE distance asc. Fixes bug where un-sorted fronts put low-satisfied individuals at index 0. |
| `operators.py` | ✅ | Single-point / uniform crossover, synonymous mutation, adaptive mute rate, elite selection. |
| `processor.py` | ✅ | Main loop with `prev_seq_0` tracking. Dumps `fitness_gen_{gen:03d}.csv` whenever seq_0 changes. **Reproduce dedup**: `existing.add(child)` prevents sibling duplicates. **Post-selection dedup**: removes duplicates after env. selection and backfills with random individuals. No convergence early-stop. Reusable Pool. Checkpoint with version validation. |
| `io_utils.py` | ✅ | `read_protein_sequence`, `read_cds_sequences`, `write_results`. `fitness.csv` no longer contains `rna_seq` column. |
| `main.py` | ✅ | Three commands: `optimize`, `resume`, `report`. YAML config support. **Critical fix**: uses `ctx.get_parameter_source()` to distinguish CLI defaults from explicit arguments; YAML values are no longer masked by Click defaults. |
| Unit + integration tests | ✅ | pytest suite. |

## Key Design Decisions

1. **Target-driven objective toggling**: `FitnessConfig` no longer has `enable_*` flags. An objective is active iff its `target` is not `None`. Only `target_cai` and `target_avg_mfe` have non-None defaults.
2. **Tolerance mechanism**: `_obj_value()` returns `0.0` for any objective value within `target ± tolerance`. NSGA-II treats such solutions as equally optimal for that objective.
3. **Satisfied-priority sorting**: Within every Pareto front, individuals are ranked by: (1) number of satisfied objectives (desc), (2) CAI distance to target (asc), (3) avg_MFE distance to target (asc). This applies to **all fronts**, not just truncated ones.
4. **`avg_MFE` not `MFE`**: The evaluator divides MFE by sequence length so the objective is length-normalized.
5. **Weighted initialization**: `ProteinSpec.random_individual(weighted=True)` samples codons by species CAI weights (sharpened with `power=10`), pushing initial CAI > 0.95. **Caution**: `power=10` is very aggressive and can cause population homogenization (initial CAI ≈ 0.96+). Consider lowering `power` or `target_cai` if diversity collapses.
6. **Pool reuse**: `GeneticAlgorithmProcessor` creates a `multiprocessing.Pool` on init and reuses it across generations.
7. **Circular import guard**: `genetic_alg/__init__.py` does **not** import `processor`. `processor.py` uses a local import for `build_objectives`.
8. **CAI2 mutability**: `cai2.calc_cai` mutates the weights dict passed to it. We always pass a deep copy or pre-computed read-only weights.
9. **Generation fitness dump**: Whenever `seq_0` (population[0]) changes, a full `fitness_gen_{gen:03d}.csv` is written with all individuals' sequences and fitness values.
10. **Population deduplication (two layers)**:
    - **Reproduce layer**: `existing = set(self.population)` + `existing.add(child)` prevents children from duplicating parents or siblings.
    - **Post-selection layer**: After `environmental_selection`, the final population is deduplicated (`set[Individual]`). Missing slots are backfilled with new random individuals to maintain `population_size`.
11. **CLI parameter resolution**: `_resolve()` uses `click.get_current_context().get_parameter_source(key)` to detect whether the user explicitly passed a CLI argument. Only explicit CLI arguments override YAML; Click defaults do not.
12. **Real reference sequences for CPB**: `tables/reference_cds/` stores genuine NCBI RefSeq CDS (e.g. human ACTB). `get_reference_sequences(species)` scans `.fa` files and raises `ValueError` for unsupported species. No synthetic random-sequence fallback.

## CLI Commands

### `cdsopt optimize`
```bash
cdsopt optimize protein.fa --config config.yaml
```
- `--config`: YAML config file
- `--init-cds`: FASTA with initial CDS sequences
- `--weighted-init`: Use species CAI weights for random initialization
- `--species`, `--pop-size`, `--generations`, `--processes`, `--seed`, etc.
- Per-objective `--target-*` and `--tolerance-*` options

### `cdsopt resume`
```bash
cdsopt resume checkpoint.pkl --processes 4
```

### `cdsopt report`
```bash
cdsopt report cds.fa --species human --fold-engine vienna -o report.csv
```
Evaluates all objective raw values for given CDS sequences (no targets involved).

## UTR Optimization Subproject (`utropt/`)

A new subproject is being developed under `mRNAOptimizer/utropt/` to handle **5' UTR / 3' UTR evaluation and design**. It is intentionally decoupled from `cdsopt/` but designed to work with it.

### Goals

1. Evaluate 5' and 3' UTRs with biologically meaningful metrics (Kozak, uORF, MFE, ARE, PAS, miRNA binding, GC, etc.).
2. Design UTR variants by either:
   - Rule-based replacement from a curated UTR library (Kozak, β-globin, bGH pA, SV40 pA, etc.)
   - Generative design via GEMORNA (`submodules/GEMORNA`) for 5'UTR / 3'UTR
3. Assemble full-length mRNA candidates from 5'UTR + CDS + 3'UTR + polyA (polyA kept unchanged).
4. Provide an independent CLI: `python -m utropt`.

### Planned Module Layout

```
utropt/
├── __init__.py
├── cli.py                  # evaluate / design / score commands
├── library.py              # Curated 5'/3' UTR parts
├── assembler.py            # assemble_mrna(five_utr, cds, three_utr, poly_a)
├── evaluator.py            # 5'/3'/full-mRNA metrics
├── designer.py             # rule-based + GEMORNA-based generation
├── utils.py                # FASTA helpers, coordinate mapping
└── tests/
    ├── test_library.py
    ├── test_assembler.py
    └── test_evaluator.py
```

### Implementation Status

| Module | Status | Notes |
|--------|--------|-------|
| `utropt/__init__.py` | ✅ | Exposes public API: `UTREvaluator`, `UTRConfig`, `UTRRecommender`, `predict_5utr`, `predict_3utr`, generators |
| `utropt/utils.py` | ✅ | RNA sequence utilities, motif scanning helpers |
| `utropt/gemorna_adapter.py` | ✅ | Direct import wrapper for GEMORNA 5'/3' UTR predictors; verified on Python 3.13 |
| `utropt/gemorna_generator.py` | ✅ | Thin wrapper around GEMORNA's official CLI (`generate_5utr_candidates`, `generate_3utr_candidates`). Linux/WSL only. CDS generation excluded; use `cdsopt`.
| `utropt/mirna_seeds.py` | ✅ | Built-in human miRNA seed DB + TargetScan/miRBase parsers |
| `utropt/rbp_motifs.py` | ✅ | ATtRACT parser and RBP motif scanner |
| `utropt/evaluator.py` | ✅ | 5'/3'/full-mRNA evaluator with Kozak, uORF, MFE, ARE, PAS, miRNA, RBP, GEMORNA scores |
| `utropt/recommender.py` | ✅ | **New** `UTRRecommender`: generate UTR combinations for a fixed CDS, score with `UTREvaluator`, rank by composite objective, write FASTA + CSV report. CDS optimization intentionally excluded; use `cdsopt`. |
| `utropt/library.py` | ⏳ | Curated UTR part library (Kozak, β-globin, bGH pA, SV40 pA, etc.) |
| `utropt/assembler.py` | ⏳ | Full-mRNA assembler preserving polyA (functionality currently folded into `recommender.py`) |
| `utropt/designer.py` | ⏳ | Rule-based + GEMORNA-based candidate generation (functionality currently folded into `recommender.py`) |
| `utropt/cli.py` | ⏳ | `evaluate` / `design` / `score` commands |
| `utropt/tests/` | ⏳ | pytest suite |

### Tool Integration

| Tool | Purpose | Integration Notes |
|------|---------|-------------------|
| **GEMORNA** | Generate/predict 5'UTR and 3'UTR | Imported directly via `utropt.gemorna_adapter`; `torch` added to `[utropt]` optional deps |
| **ViennaRNA / LinearFold** | RNA folding | Reuse `cdsopt.utils.fold_tools` |
| **ATtRACT** | RBP motif scanning | Parser in `utropt.rbp_motifs`; expects `data/ATtRACT/ATtRACT_db.txt` |
| **TargetScan / miRBase data** | miRNA target prediction | Seed-match scanning with built-in fallback DB |
| **Rule-based scanners** | Kozak, uORF, ARE, PAS | Built-in |

### Design Cases

For each of the two YingXuProject mRNAs, generate candidates for:

| Group | 5' UTR | CDS | 3' UTR | Notes |
|------|--------|-----|--------|-------|
| G1 (baseline) | Original | Original | Original | Single sequence, evaluated for comparison |
| G2 | GEMORNA | Original | GEMORNA | Optimized UTR only (utropt scope) |

Implemented in `YingXuProject/run_utropt_gemorna.py` using `utropt.UTRRecommender`.
UTR lengths: 5'UTR generated in both `short` and `medium`; 3'UTR generated as `long`.
Each group recommends top 3, yielding 2 mRNAs × 2 groups × 3 = 12 candidate sequences (G1 contributes 1 baseline).

CDS optimization is handled by `cdsopt`, not utropt.

## Current State

- **All core algorithm modules implemented and tested**.
- **CLI fully functional**: optimize / resume / report commands.
- **YAML config support**: `--config` flag. CLI explicit arguments override YAML; Click defaults no longer mask YAML values.
- **Population deduplication**: Both reproduce and post-selection layers active; every generation ends with `unique=population_size`.
- **CPB reference sequences**: Real ACTB reference sequence from NCBI replaces synthetic random references.
- **Result outputs**: `pareto_front.fasta`, `fitness.csv`, `summary.json`, plus per-generation `fitness_gen_{gen:03d}.csv`.
- **LinearFold Windows build**: Compiled `linearfold_v.exe` / `linearfold_c.exe` with MinGW; `cdsopt` auto-detects via `linearfold.bat` wrapper.
- **Python 3.10 compatibility**: `pyproject.toml` now requires `>=3.10` and pins GEMORNA-compatible dependency versions (torch==2.2.0, numpy==1.26.4, pandas==2.0.2, biopython>=1.72). Windows development can still use Python 3.13, but GEMORNA generation must run under Python 3.10 (e.g. WSL conda env).
- **LinearFold submodule migration**: Switched to `DStarEpoch/LinearFork` fork with Windows MinGW support; submodule pointer updated and pushed.
- **`report` command extended**: Can now analyze full-length mRNA sequences; non-CDS sequences only receive structural metrics.
- **`utropt/` evaluation module implemented**: `UTREvaluator` supports 5'UTR (Kozak, uORF, MFE, GEMORNA TIE), 3'UTR (PAS, ARE, miRNA, RBP, MFE, GEMORNA stability), and full-mRNA metrics. Verified on RNA editing-1 5'UTR.
- **GEMORNA integration**: Registered as Git submodule; direct Python import via `utropt.gemorna_adapter` works under Python 3.13; `torch` added to `[utropt]` optional dependencies.
- **RBP motif data**: ATtRACT database added under `data/ATtRACT/` and wired into `utropt.rbp_motifs`.
- **Known active issue — Population homogenization**: `weighted_init` with `power=10` produces initial CAI ≈ 0.96+. When `target_cai` is strict (e.g. 0.97 ± 0.01), nearly all individuals start outside tolerance, causing NSGA-II to collapse into a single CAI-driven front within ~15 generations. Post-selection dedup masks the symptom by injecting random individuals, but the core convergence trend is toward monoculture. Remedies: lower `power` (e.g. 3–5), raise `target_cai`, or increase `mute_rate`.
- **No visualization**: Pareto front plots, convergence curves not implemented.

## Future Tasks

### High Priority
1. ~~Implement CLI (`main.py`)~~ ✅
2. ~~Result output & serialization~~ ✅
3. ~~Integration / smoke test~~ ✅
4. ~~YAML config file support~~ ✅
5. ~~Weighted initialization~~ ✅
6. ~~Init CDS support~~ ✅
7. ~~Target-driven objective toggling~~ ✅
8. ~~Population deduplication~~ ✅
9. ~~Fix config parameter resolution (YAML vs Click defaults)~~ ✅
10. ~~Fix environmental_selection front sorting~~ ✅
11. ~~Real reference sequences for CPB~~ ✅
12. ~~Implement `utropt/` evaluation module~~ ✅
    - Rule-based metrics (Kozak, uORF, ARE, PAS, miRNA, RBP)
    - GEMORNA 5'/3' UTR prediction integration
    - Full-mRNA structural evaluation
13. **Complete `utropt/` design workflow**
    - ~~Implement GEMORNA generation wrapper~~ ✅ (`utropt/gemorna_generator.py`)
    - ~~Implement candidate combination + scoring + recommendation~~ ✅ (`utropt/recommender.py`)
    - Implement `library.py` with curated 5'/3' UTR parts
    - Implement standalone `assembler.py` / `designer.py` (currently folded into `recommender.py`)
    - Implement `cli.py` with `evaluate` / `design` / `score` commands
    - Validate full 4-group design matrix on RNA editing-1 / RNA editing-2

### Medium Priority
14. **Address population homogenization**
    - Reduce `weighted_init` power or make it configurable
    - Consider adaptive `target_cai` or diversity-preserving selection
    - Investigate whether `existing.add(child)` + post-selection dedup is sufficient, or if environmental selection itself should penalize duplicates

15. **Convergence visualization**
    - Plot Pareto front evolution (e.g. CAI vs. avg_MFE over generations)
    - Plot hypervolume or best-front summary over time
    - Optional Matplotlib / Plotly dependency

16. **Sequence constraints**
    - Enforce specific codons at given positions (e.g. start codon, restriction sites)
    - Avoid specific motifs (e.g. cryptic splice sites, Shine-Dalgarno-like sequences)
    - `not_mutate_idx` already exists in `GAConfig`; extend to motif avoidance

### Low Priority / Nice to Have
17. **Alternative fold engines**
    - Integrate RNAfold (ViennaRNA CLI) as another engine option
    - Support LinearFold partition function mode if available

18. **Performance**
    - Profile batch evaluation bottleneck
    - Consider caching RNA fold results across generations
    - Numba / Cython for hot loops if needed

19. **Documentation**
    - API docs (sphinx / mkdocs)
    - Example notebooks

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
