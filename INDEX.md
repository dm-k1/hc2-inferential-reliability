# Repository Index

## Project Purpose

> **Note**: This repository contains ongoing research shared for master's program admissions review.

This project investigates a fundamental question in applied regression: **When can practitioners trust classical standard errors versus robust HC standard errors?** 

The work explores an empirically calibrated reliability score that aims to quantify heteroskedasticity severity independent of sample size, potentially enabling practitioners to make more informed decisions about inference validity.

---

## 📊 Primary Deliverables (HTML Reports)

These three HTML files contain the current analysis with interactive visualizations, statistical output, and interpretation. **These represent the current state of the research.**

### [00_HC_estimators_validation.html](00_HC_estimators_validation.html)

**Purpose**: Examine matrix-based OLS and HC implementation against R's reference functions

**What it contains**:
- Numerical validation across diverse test scenarios (homoskedastic, heteroskedastic, collinear)
- Verification that differences are at machine precision (< 10^-10)
- Distribution plots of relative errors
- Confirmation that the fast engine is algebraically equivalent to lm() + vcovHC()

**Status**: Implementation examined against reference benchmarks

---

### [01_null_calibration.html](01_null_calibration.html)

**Purpose**: Establish the reference distribution of the inferential score under homoskedasticity

**What it contains**:
- 100,000+ simulation analysis across N ∈ {50, 100, ..., 102,400}
- Comparison of HC1, HC2, and HC3 estimators
- Analysis of how 95th percentiles scale with sample size
- Selection of HC2 as the most stable estimator
- Finite-sample correction factor calibration

**Preliminary observations**:
- HC2 shows promising N-invariance properties
- Tentative cutoff value ~1.64 explored for hypothesis testing
- Scaling appears to follow 1/√N pattern

**Status**: HC2 selected for further investigation based on observed properties

---

### [02_HC2_validation.html](02_HC2_validation.html)

**Purpose**: Explore HC2 behavior under heteroskedasticity and develop N-adjusted reliability score

**What it contains**:
- Systematic heteroskedastic simulation grid (300,000+ simulations)
- Demonstration of N-dependence problem with raw inferential score
- Development of reliability score with finite-sample correction (c ≈ 0.23)
- Three-panel visualization showing 47.6% reduction in N-dependence
- Formal ANOVA testing of N-invariance properties
- Coverage degradation analysis under heteroskedasticity

**Preliminary observations**:
- Reliability score appears to reduce N-dependence substantially (though complete invariance not yet achieved)
- Remaining N-effect may suggest higher-order finite-sample effects requiring further investigation
- Performance appears comparable to standard heteroskedasticity tests in initial comparisons

**Status**: Reliability score under development with ongoing evaluation

---

## 📝 Source Notebooks (Rmd)

These R Markdown files are the source code that generated the HTML reports above.

- **00_HC_estimators_validation.Rmd** – Implementation validation source
- **01_null_calibration.Rmd** – Null calibration source  
- **02_HC2_validation.Rmd** – Reliability score development source

---

## 🔧 Supporting R Functions (R/)

Modular codebase organized by analytical stage:

**R/00_config.R**  
Central configuration: sample size grids, HC types, random seed, parallel backend setup

**R/10_dgp_and_fits.R**  
Data generation and high-performance OLS/HC fitting engine using optimized matrix algebra

**R/20_metrics.R**  
Computation of inferential scores, reliability metrics, and coverage statistics

**R/30_null_calibration.R**  
Null calibration simulation routines and parallel execution wrappers

**R/40_hetero_sims.R**  
Heteroskedastic simulation grid and aggregation functions

**R/50_invariance_and_lookup.R**  
N-invariance analysis utilities and foundations for future lookup table approaches

**R/60_tables_and_plots.R**  
Visualization functions (all plots used in HTML reports) and presentation utilities

---

## 💾 Results Directory (results/)

Cached simulation outputs for reproducibility without re-computation:

**From Part 1 (Null Calibration)**:
- `best_hc.rds` – Selected HC estimator (HC2)
- `scores_long.rds` – Full null calibration scores
- `q95_summary.rds` – 95th percentile summaries
- `scaling_results.rds` – N-scaling analysis
- `hc_ranking.csv` – HC estimator comparison table

**From Part 2 (Reliability Score)**:
- `hetero_sim_results.rds` – Full heteroskedastic simulation data
- `hetero_aggregated.rds` – Aggregated summary statistics
- `coverage_adjusted_by_N.rds` – Coverage analysis results
- `universal_lookup_table.rds` – Lookup table foundations

---

## 📚 Documentation

**[README.md](README.md)** – Project overview, quick start, and navigation guide

**[INDEX.md](INDEX.md)** – This file: detailed component descriptions

**[THEORY.md](THEORY.md)** – Mathematical foundations, methodology, and theoretical background

**[CHANGELOG.md](CHANGELOG.md)** – Development history and version notes

---

## ⚡ Quick Navigation

**Want to understand the project?** → Start with [README.md](README.md)

**Want to see the results?** → Open the HTML files (start with [Part 0](00_HC_estimators_validation.html))

**Want to understand the theory?** → Read [THEORY.md](THEORY.md)
