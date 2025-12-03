# Repository Index

**Last Updated**: 2 December 2025

---

## Core Notebooks

### 01_null_calibration.Rmd
**Purpose**: Null calibration of the inferential score $ S_{\mathrm{Inf}} $ under the true null hypothesis where X and Y are completely independent.

**Null DGP**: $ X_i \sim N(0,1), Y_i \sim N(0,1), \text{independent} $. Despite independence, OLS models $ Y \sim 1 + X $ are fitted to establish the null reference distribution of HC estimators under correct specification but zero true explanatory power.

**Workflow**:
- Step 1: Simulate 100,000 OLS models under true null DGP across sample sizes $ N \in \{50, 100, 200, \ldots, 102400\} $
- Step 2: Compute $ S_{\mathrm{Inf}} $ for HC1, HC2, and HC3 estimators applied post-hoc
- Step 3: Extract Q95 quantiles per sample size and HC type
- Step 4: Test scaling exponent grid $ \alpha \in \{0.30, 0.35, \ldots, 0.70\} $
- Step 5: Select optimal HC type via Q95 stability analysis
- Visualization and summary statistics

**Inputs**: Configuration from `R/00_config.R`  
**Outputs**:
- `results/null_calibration_sim_results.rds` – Simulation fits and scores
- `results/null_calibration_summary.rds` – Q95 and scaling analysis
- `results/null_scaling_results.csv` – HC type ranking (stable name)

**Runtime**: ~20 minutes

---

### 02_HC2_validation.Rmd
**Purpose**: Heteroskedastic validation, power analysis, and universal lookup table construction for the selected HC estimator.

**Workflow**:
- Step 1: Load selected HC estimator from null calibration.
- Step 2: Simulate heteroskedastic DGPs ($ Y = 1 + \beta_x X + \varepsilon $, $ \sigma^2 = \sigma_0^2 + \lambda X^2 $) across a grid of $ N $, $ \lambda $, $ \beta_x $, and $ \sigma_0 $.
- Step 3: Evaluate power of $ S_{\mathrm{Inf}} $ to detect heteroskedasticity (benchmarked against Breusch-Pagan and White tests).
- Step 4: Demonstrate N-dependence of $ S_{\mathrm{Inf}} $ and N-invariance of $ S_{\mathrm{Adj}} $.
- Step 5: Construct universal lookup table mapping $ S_{\mathrm{Adj}} $ to expected coverage loss.

**Inputs**:
- `results/best_hc.rds` – Selected HC type
- Configuration: `N_GRID_HETERO`, `HETERO_STRENGTH_GRID`, `BETA_X_GRID`, `SIGMA0_GRID`

**Outputs**:
- `results/hetero_sim_results.rds` – Raw simulation data
- `results/hetero_aggregated.rds` – Aggregated cell statistics
- `results/coverage_adjusted_by_N.rds` – Wide-format lookup (Coverage × N)
- `results/universal_lookup_table.rds` – Final universal lookup table

---

## R Support Functions

### R/00_config.R
Global configuration parameters:
- `N_GRID_NULL` – Sample size grid for null calibration
- `N_SIMS_PER_N_NULL` – Simulations per sample size
- `HC_TYPES` – HC estimators to test (c("HC1", "HC2", "HC3"))
- `ALPHA_GRID` – Scaling exponent candidates (seq(0.3, 0.7, by=0.05))

### R/10_dgp_and_fits.R
Data generation and model fitting helpers. Includes `simulate_heteroskedastic_dgp` (quadratic variance) and `run_hetero_simulation_fast` (optimized matrix simulation).

### R/20_metrics.R
Metric computations:
- Classical standard error
- HC standard errors (HC1, HC2, HC3)
- Inferential score $ S_{\mathrm{Inf}} $ and $ S_{\mathrm{Adj}} $
- Scaling analysis

### R/30_null_calibration.R
**Modular null calibration functions with adaptive batching and checkpointing** (4 independent functions):

| Function | Purpose | Batching |
|----------|---------|----------|
| `simulate_fits_and_save()` | Generate homoskedastic fits in batches; parallel execution | ✓ Adaptive (init=200, max=2000) |
| `compute_scores_and_save()` | Compute $ S_{\mathrm{Inf}} $ per HC type in batches | ✓ Adaptive (init=100, max=1000) |
| `apply_scaling_and_save()` | Test $ \alpha $ grid; identify best scaling | Sequential |
| `select_best_hc()` | Rank HC types by Q95 stability | Sequential |

**Adaptive Batching**:
- Batches start small (200–100 fits/scores) and double each iteration
- Hard cap at 2000–1000 to prevent memory overload
- All batches execute in parallel via `furrr::future_map()` with `multicore` backend
- Checkpoint saves after each batch (`fits_N_{N}_batch_{num}.rds`, `scores_N_{N}_batch_{num}.rds`)
- Enables fault-tolerant restart capability

**Parallelization**:
- Global seed: `set.seed(2024)` in `R/00_config.R` (only location)
- Parallel workers: Independent seeds via `furrr_options(seed = TRUE)`
- Backend: `multicore` with `workers = availableCores() - 1`
- Reproducibility: Same global seed → identical results across runs

Each function:
- Handles one analytical step
- Saves intermediate results to `results/` (plus checkpoint saves)
- Displays progress/results
- Pure function (no side effects)

### R/40_hetero_sims.R
Heteroskedastic simulation grid routines. Implements `run_hetero_sim_grid` (parallelized over cells) and `aggregate_hetero_sims`.

### R/50_invariance_and_lookup.R
Invariance analysis and lookup table construction. Includes `run_per_coverage_regressions`, `fit_global_mixed_model`, and `build_universal_lookup_table`.

### R/60_tables_and_plots.R
Presentation utilities (tables, figures)

---

## Results Directory Structure

| File | Source | Purpose |
|------|--------|---------|
| `null_calibration_sim_results.rds` | 01_null_calibration | Fits and inferential scores |
| `null_calibration_summary.rds` | 01_null_calibration | Q95 and scaling analysis |
| `null_scaling_results.csv` | 01_null_calibration | HC type ranking (used downstream) |
| `hetero_sim_results.rds` | 02_HC2_validation | Heteroskedastic simulation outputs |
| `universal_lookup_table.rds` | 02_HC2_validation | Final calibrated lookup table |

---

## Governance

See `AGENTS.md` for repository policies:
- Minimal, stable documentation set (README, CHANGELOG, INDEX only)
- Academic tone; methods-review appropriate language
- Refactoring preferred over duplication
- Reproducibility and traceability maintained
- Configuration-driven parameters; stable result names

---

## Quick Links

- **How to run**: Edit `R/00_config.R`, then `rmarkdown::render("01_null_calibration.Rmd")`
- **View results**: After run, inspect `results/null_scaling_results.csv`
- **Next stage**: Pass selected HC to `02_HC2_validation.Rmd`
- **Repository policies**: See `AGENTS.md`
