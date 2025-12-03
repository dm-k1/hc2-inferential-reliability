# Repository Index

**Last Updated**: 3 December 2025

---

## Core Notebooks

### 00_HC_estimators_validation.Rmd
**Purpose**: Validation of the matrix-based OLS and HC covariance estimator implementation against reference implementations (`lm()` + `vcovHC()`).

**Status**: LOCKED (Part 0)

**Workflow**:
- Step 1: Generate data under homoskedastic, heteroskedastic, and collinear DGPs
- Step 2: Compute OLS coefficients and SEs using both matrix-based and reference implementations
- Step 3: Compare numerical equivalence (absolute and relative differences)
- Step 4: Validate HC1, HC2, HC3 estimators across sample sizes N ∈ {20, 50, 100, 1000}

**Inputs**: Configuration from `R/00_config.R`
**Outputs**: Validation tables demonstrating numerical equivalence (differences < 1e-10)

**Runtime**: ~2 minutes

---

### 01_null_calibration.Rmd
**Purpose**: Null calibration of the inferential score $ S_{\mathrm{Inf}} $ under the true null hypothesis where X and Y are completely independent.

**Status**: LOCKED (Part 1)

**Null DGP**: $ X_i \sim N(0,1), Y_i \sim N(0,1), \text{independent} $. Despite independence, OLS models $ Y \sim 1 + X $ are fitted to establish the null reference distribution of HC estimators under correct specification but zero true explanatory power.

**Workflow**:
- Step 1: Simulate 10,000 OLS models under true null DGP across sample sizes $ N \in \{50, 100, 200, \ldots, 102400\} $
- Step 2: Compute $ S_{\mathrm{Inf}} $ for HC1, HC2, and HC3 estimators applied post-hoc
- Step 3: Extract Q95 quantiles per sample size and HC type
- Step 4: Test scaling exponent grid $ \alpha \in \{0.30, 0.35, \ldots, 0.70\} $
- Step 5: Select optimal HC type via Q95 stability analysis
- Visualization and summary statistics

**Inputs**: Configuration from `R/00_config.R`  
**Outputs**:
- `results/scores_long.rds` – Simulation scores
- `results/q95_summary.rds` – Q95 quantiles per N and HC type
- `results/scaling_results.rds` – Scaling analysis results
- `results/best_hc.rds` – Selected HC type
- `results/hc_ranking.csv` – HC type ranking

**Runtime**: ~20 minutes

---

### 02_HC2_validation.Rmd
**Purpose**: Heteroskedastic validation, power analysis, and universal lookup table construction for the selected HC estimator.

**Status**: IN PROGRESS (Part 2)

**Workflow**:
- Step 1: Load selected HC estimator from null calibration.
- Step 2: Simulate heteroskedastic DGPs ($ Y = 1 + \beta_x X + \varepsilon $, $ \sigma^2 = \sigma_0^2 + \lambda X^2 $) across a grid of $ N $, $ \lambda $, $ \beta_x $, and $ \sigma_0 $.
- Step 3: Evaluate power of $ S_{\mathrm{Inf}} $ to detect heteroskedasticity (benchmarked against Breusch-Pagan and White tests).
- Step 4: Demonstrate N-dependence of $ S_{\mathrm{Inf}} $ and N-invariance of $ S_{\mathrm{Rel}} $.
- Step 5: Construct universal lookup table mapping $ S_{\mathrm{Rel}} $ to expected coverage loss.

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
- `N_SIMS_PER_N_NULL` – Simulations per sample size (10,000)
- `HC_TYPES` – HC estimators to test (c("HC1", "HC2", "HC3"))
- `ALPHA_GRID` – Scaling exponent candidates (seq(0.3, 0.7, by=0.05))
- Parallel execution setup via `furrr`

### R/10_dgp_and_fits.R
Data generation and model fitting:
- `simulate_homoskedastic_dgp()` – True null DGP (independent X, Y)
- `simulate_heteroskedastic_dgp()` – Quadratic variance DGP
- `fit_ols_hc()` – Matrix-based OLS + HC estimation
- `run_null_simulation_fast()` – Optimized null calibration simulation
- `run_hetero_simulation_fast()` – Optimized heteroskedastic simulation

### R/20_metrics.R
Canonical metric definitions:
- `compute_sr_inf()` – $ S_{\mathrm{Inf}} = \text{SE}_{\text{robust}} / \text{SE}_{\text{classic}} - 1 $
- `compute_sr_inf_adj()` – $ S_{\mathrm{Inf}}^* = \sqrt{N} \cdot S_{\mathrm{Inf}} $
- `compute_sr_ratio()` – $ \text{sr\_ratio} = \text{SE}_{\text{classic}} / \text{SE}_{\text{robust}} $
- `compute_sr_ratio_adj()` – $ S_{\mathrm{Rel}} = \text{sr\_ratio} - 2/(3\sqrt{N}) $
- `compute_se_metrics()` – Full metric computation from lm fit

### R/30_null_calibration.R
Modular null calibration functions with parallelized execution:

| Function | Purpose |
|----------|---------|
| `run_fast_null_calibration()` | Generate null simulations in parallel with fixed batching (N=500) |
| `apply_scaling_and_save()` | Test scaling exponents; compute Q95 statistics |
| `select_best_hc()` | Rank HC types by Q95 stability |

**Parallelization**:
- Fixed batch size of 500 simulations per parallel task
- Global seed in `R/00_config.R`; worker seeds via `furrr_options(seed = TRUE)`
- Backend: `multicore` (Unix) or `multisession` (Windows)

### R/40_hetero_sims.R
Heteroskedastic simulation grid:
- `run_hetero_sim_grid()` – Parallelized grid simulation
- `aggregate_hetero_sims()` – Cell-level aggregation

### R/45_inference_breakage_sims.R
Inference breakage analysis (exponential heteroskedasticity):
- `simulate_heteroskedastic_dgp_exp()` – Exponential variance DGP
- `run_inference_breakage_sim()` – Coverage degradation simulation
- `aggregate_breakage_results()` – Summary statistics
- `plot_inference_breakage()` – Visualization

### R/50_invariance_and_lookup.R
Invariance analysis and lookup table construction:
- `build_coverage_adjusted_table()` – Wide-format coverage table
- `compute_spread_per_coverage()` – N-invariance spread analysis
- `run_per_coverage_regressions()` – Per-stratum regression tests
- `fit_global_mixed_model()` – Mixed-effects model for N-dependence
- `build_universal_lookup_table()` – Final lookup table construction

### R/60_tables_and_plots.R
Presentation utilities:
- `print_kable_table()` – Pretty table output
- `save_table_csv()`, `save_rds()` – Result persistence
- `plot_reliability_score_vs_coverage()` – Lookup table visualization
- `plot_reliability_score_by_n()` – Reliability score vs N by coverage gap
- `plot_power_curves()` – Power analysis visualization
- `plot_n_dependence()` – N-dependence of test statistic
- `plot_n_invariance_comparison()` – N-invariance comparison plots
- `plot_benchmark_comparison()` – S_Inf vs BP/White power comparison
- `plot_universal_lookup()` – Universal lookup table visualization
- `save_plot_png()` – Figure export

---

## Results Directory Structure

| File | Source | Purpose |
|------|--------|---------|
| `scores_long.rds` | 01_null_calibration | Raw simulation scores |
| `q95_summary.rds` | 01_null_calibration | Q95 quantiles per N and HC |
| `scaling_results.rds` | 01_null_calibration | Scaling analysis |
| `best_hc.rds` | 01_null_calibration | Selected HC type |
| `hc_ranking.csv` | 01_null_calibration | HC type ranking |
| `hetero_sim_results.rds` | 02_HC2_validation | Heteroskedastic simulation outputs |
| `hetero_aggregated.rds` | 02_HC2_validation | Aggregated statistics |
| `universal_lookup_table.rds` | 02_HC2_validation | Final calibrated lookup table |
| `universal_lookup_table.csv` | 02_HC2_validation | CSV export of lookup table |

---

## Governance

See `AGENTS.MD` for repository policies:
- Minimal, stable documentation set (README, CHANGELOG, INDEX only)
- Academic tone; methods-review appropriate language
- Refactoring preferred over duplication
- Reproducibility and traceability maintained
- Configuration-driven parameters; stable result names

---

## Quick Links

- **How to run**: Edit `R/00_config.R`, then `rmarkdown::render("01_null_calibration.Rmd")`
- **View results**: After run, inspect `results/hc_ranking.csv`
- **Next stage**: Pass selected HC to `02_HC2_validation.Rmd`
- **Repository policies**: See `AGENTS.MD`
