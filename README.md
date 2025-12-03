# HC2-Based Inferential Reliability Metric for Linear Regression

## Overview

This repository implements and validates an inferential reliability metric for linear regression based on heteroskedasticity-consistent (HC) covariance matrix estimators. The core contribution concerns the development and empirical calibration of a finite-sample reliability score, $ S_{\mathrm{Inf}} $, designed to assess the inferential robustness of HC-based interval estimates under model misspecification.

## Project Structure

### Core Analysis Pipeline

1. **Null Calibration** (`01_null_calibration.Rmd`)
   - Establishes the null reference distribution of $ S_{\mathrm{Inf}} $ under the true null hypothesis where X and Y are completely independent
   - Generates data under true null DGP: $ X_i \sim N(0,1), Y_i \sim N(0,1), \text{independent} $
   - Fits 100,000 OLS models across sample sizes despite absence of true relationship
   - Applies HC1, HC2, and HC3 estimators post-hoc to establish their null behavior
   - Computes inferential scores for each HC type
   - Selects the optimal HC estimator via Q95-based stability analysis across sample sizes

2. **HC2 Validation** (`02_HC2_validation.Rmd`)
   - Evaluates the selected HC estimator under heteroskedastic data-generating processes
   - Constructs a lookup table mapping sample size and effect size to critical thresholds
   - Validates coverage properties across controlled heteroskedasticity scenarios

### Supporting Structure

- **R/** – Modular utility functions organized by analytical stage:
  - `R/00_config.R` – Global configuration (HC types, alpha grid, sample size grid)
  - `R/10_dgp_and_fits.R` – Data generation and model fitting
  - `R/20_metrics.R` – Metric computations
  - `R/30_null_calibration.R` – Null calibration routines (4 focused functions)
  - `R/40_hetero_sims.R` – Heteroskedastic simulation grid
  - `R/50_invariance_and_lookup.R` – Invariance analysis and lookup table construction
  - `R/60_tables_and_plots.R` – Presentation utilities

- **results/** – Analytical artifacts (serialized RDS and CSV outputs)
- **CHANGELOG.md** – Versioned record of substantive changes
- **AGENTS.md** – Governance policies for automated contributions

## Null Calibration Workflow

The null calibration implements the following sequential pipeline:

```
Step 1: Simulate fits      → Generate 100,000 homoskedastic models
Step 2: Compute scores     → Calculate S_Inf per HC type and sample size
Step 3: Compute Q95        → Extract 95th percentile of scores per N
Step 4: Test scaling       → Find optimal exponent α (decoupled)
Step 5: Select HC          → Rank HC types by Q95 stability
Step 6: Save results       → Persist workflow outputs
```

### Key Design Decisions

**Decoupled Alpha Estimation**: The scaling exponent $ \alpha $ is estimated independently of the HC type selection. This separation ensures that the choice of HC estimator reflects inferential score stability, not exponent estimation artifacts.

**Inferential Score Focus**: The workflow emphasizes $ S_{\mathrm{Inf}} $ exclusively, defined as the normalized deviation between the HC-based confidence interval width and the classical interval width under homoskedasticity.

**Configuration-Driven Parameters**: Sample size grids, HC type sets, and alpha search ranges are defined centrally in `R/00_config.R`, allowing systematic sensitivity analysis without code modification.

## Quick Start

### 1. Configuration

Edit `R/00_config.R` to set:
```r
N_GRID_NULL <- seq(50, 500, by = 25)      # Sample size grid
N_SIMS_PER_N_NULL <- 10000                # Fits per sample size
HC_TYPES <- c("HC1", "HC2", "HC3")        # HC estimators to test
ALPHA_GRID <- seq(0.3, 0.7, by = 0.05)    # Scaling exponent candidates
```

### 2. Run Null Calibration

```r
rmarkdown::render("01_null_calibration.Rmd")
```

This produces:
- `results/null_calibration_sim_results.rds` – Simulation and score data
- `results/null_calibration_summary.rds` – Q95 and scaling results
- `results/null_scaling_results.csv` – HC rankings for downstream use

### 3. Inspect Results

```r
scaling_results <- read.csv("results/null_scaling_results.csv")
print(scaling_results)  # Displays HC type ranking by stability
```

### 4. Run HC2 Validation (Next Stage)

```r
rmarkdown::render("02_HC2_validation.Rmd")
```

## Computational Requirements

- **Time**: ~20 minutes for null calibration (100,000 models across 12 sample sizes)
- **Memory**: ~2–3 GB for intermediate data structures
- **Storage**: ~500 MB for serialized results

## Repository Governance

This repository is maintained under policies documented in `AGENTS.md`. Key principles:

- **Minimal, stable documentation set**: Exactly one `README.md`, `CHANGELOG.md`, and `INDEX.md`
- **Academic tone**: Language appropriate for thesis reviewers and quantitative researchers
- **Non-redundancy**: Refactoring preferred over near-duplicate files
- **Reproducibility**: Seeds and configurations centralized and version-controlled

## Related Work

The HC estimators HC1, HC2, and HC3 follow the definitions in:
- Huber, P. J. (1967). "The behavior of maximum likelihood estimates under nonstandard conditions."
- Eicker, F. (1967). "Limit theorems for regressions with unequal and dependent errors."
- MacKinnon, J. G., & White, H. (1985). "Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties."

The Q95-based scaling approach draws on:
- Kiefer, N. M., & Vogelsang, T. J. (2005). "A new approach to heteroskedasticity robust testing."

## License

[Specify as appropriate for your institution]

## Contact

For methodological questions or code issues, contact the repository maintainers.
