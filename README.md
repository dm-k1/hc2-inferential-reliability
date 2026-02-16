# HC2-Based Inferential Reliability Score for Linear Regression

> **Note**: This repository contains ongoing research shared for master's program admissions review. The work presented here represents current progress on an active research project.

## Purpose

**The Central Question**: When should practitioners trust classical standard errors versus robust (HC) standard errors in linear regression?

Heteroskedasticity violates classical regression assumptions, making standard errors unreliable. While robust HC estimators fix this, they produce wider confidence intervals. This project explores the tradeoff between validity and precision.

## Overview

This project investigates a **finite-sample reliability score** for assessing heteroskedasticity severity in linear regression. The approach involves empirical calibration to develop a sample-size-independent diagnostic that may help practitioners determine:
- When classical inference might remain reliable (mild heteroskedasticity)
- When robust inference appears essential (severe heteroskedasticity)
- How inferential reliability potentially degrades under misspecification

The analysis employs simulation (300,000+ replications) to calibrate the HC2 estimator and explore practical cutoff values.

**📊 View the analysis:** [Part 0: Implementation Validation](00_HC_estimators_validation.html) • [Part 1: Null Calibration](01_null_calibration.html) • [Part 2: Reliability Score Development](02_HC2_validation.html)

**For detailed theory and methodology**, see [THEORY.md](THEORY.md).

---

## Current Progress

The project includes three HTML reports documenting the ongoing analysis:

### [Part 0: OLS and HC Estimator Validation](00_HC_estimators_validation.html)
• Examines matrix-based implementation against R's reference functions  
• Tests numerical equivalence at machine precision across diverse test scenarios  
• Investigates computational approach for downstream analysis

### [Part 1: Null Calibration](01_null_calibration.html)
• Explores reference distribution under homoskedasticity (100,000+ simulations)  
• Compares HC1, HC2, and HC3 estimators for N-invariance properties  
• Preliminary selection of HC2 based on observed stability across sample sizes (N = 50 to 102,400)

### [Part 2: HC2 Validation & Reliability Score](02_HC2_validation.html)
• Evaluates HC2 behavior under systematic heteroskedasticity  
• Develops N-adjusted reliability score with finite-sample correction  
• Observes 47.6% reduction in N-dependence (though complete invariance not yet achieved)  
• Analyzes coverage degradation patterns under heteroskedasticity

---

## Project Structure

### Analysis Reports (HTML)

The three HTML files contain the current analysis with embedded visualizations, statistical output, and detailed interpretation. These represent the current state of the research.

### Source Notebooks (Rmd)

- `00_HC_estimators_validation.Rmd` – Implementation validation source
- `01_null_calibration.Rmd` – Null calibration source  
- `02_HC2_validation.Rmd` – Reliability score development source

### Supporting Code Modules (R/)

- **00_config.R** – Central configuration (grids, parameters, seed, parallel setup)
- **10_dgp_and_fits.R** – Data generation and high-performance OLS/HC engine
- **20_metrics.R** – Inferential score and reliability metric computation
- **30_null_calibration.R** – Null calibration simulation pipeline
- **40_hetero_sims.R** – Heteroskedastic simulation grid
- **50_invariance_and_lookup.R** – N-invariance analysis utilities
- **60_tables_and_plots.R** – Visualization and presentation functions

### Results (results/)

All simulation outputs (RDS/CSV format) provided for transparency.

---

## Viewing the Analysis

Open the HTML files in your browser:
- [00_HC_estimators_validation.html](00_HC_estimators_validation.html)
- [01_null_calibration.html](01_null_calibration.html)
- [02_HC2_validation.html](02_HC2_validation.html)

---

## Documentation

- **[README.md](README.md)** (this file) – Project overview and quick start
- **[INDEX.md](INDEX.md)** – Detailed component index
- **[THEORY.md](THEORY.md)** – Mathematical foundations and methodology
- **[CHANGELOG.md](CHANGELOG.md)** – Development history

---

## Status

🔬 **Ongoing Research**  
- Three analysis stages currently documented
- Codebase and configuration under active development
- Preliminary results available in results/ directory
- Shared for admissions committee review

**Last Updated**: February 2026
