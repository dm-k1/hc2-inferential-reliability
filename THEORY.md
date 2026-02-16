# Theoretical Framework: HC2-Based Inferential Reliability

## Introduction

### Why This Matters

Heteroskedasticity is common in applied regression but often ignored. Most practitioners use classical inference by default, while robust inference (HC) is available but produces wider confidence intervals. The critical question is: **How much wider must confidence intervals be to ensure validity under heteroskedasticity, and how can we quantify this tradeoff?**

This project addresses three fundamental gaps:

1. **Quantifies the cost of robustness**: Measures how much wider HC-based intervals are compared to classical intervals
2. **Contextualizes this cost**: Maps heteroskedasticity severity to inferential reliability loss
3. **Develops a diagnostic metric**: Creates a sample-size-independent reliability score for practitioners

### Project Contribution

I am developing an empirically calibrated reliability score based on the HC2 estimator that aims to:
- Correct for finite-sample bias using data-driven calibration (preliminary c ≈ 0.25)
- Achieve approximate N-invariance across sample sizes
- Provide guidance on when classical inference might remain reliable
- Bridge asymptotic theory (where HC and classical converge at infinite N) with finite-sample practice

This work explores calibrated cutoff values (~1.645) and potential lookup mechanisms to assess heteroskedasticity severity.

---

## 1. The Problem of Heteroskedasticity

### 1.1 Classical OLS Assumptions

In a standard linear regression model:

    y = X*beta + epsilon

where:
- y is an n × 1 vector of outcomes
- X is an n × k design matrix (including an intercept column)
- beta is a k × 1 coefficient vector  
- epsilon is an n × 1 error vector

Classical statistical inference relies on **homoskedasticity**: E[epsilon_i² | X] = sigma² for all observations i.

Under this assumption, the OLS estimator is:

    beta_hat = (X'X)^(-1) * X' * y

with covariance matrix:

    Var(beta_hat) = sigma² * (X'X)^(-1)

and standard errors can be estimated by plugging in:

    sigma_hat² = sum(residuals_i²) / (n - k)

Thus, classical confidence intervals [beta_hat ± z_critical * SE_hat] achieve their nominal coverage.

### 1.2 When Homoskedasticity Fails

In many real applications, error variance depends on the predictors:

    E[epsilon_i² | X] = sigma_i² = f(X_i)

Examples include:
- **Cross-sectional data**: Wealthier households may have more variable spending
- **Financial data**: Asset returns exhibit volatility clustering  
- **Survey data**: Income inequality increases with region size

**Critical consequence**: The OLS coefficient beta_hat remains unbiased and consistent, but the classical SE formula systematically underestimates (or overestimates) true variability. Confidence intervals become invalid—typically too narrow—leading to inflated Type I error rates.

---

## 2. Heteroskedasticity-Consistent (HC) Covariance Estimators

### 2.1 The Sandwich Formula

The solution is the **sandwich formula** for valid inference under arbitrary heteroskedasticity:

    Var_HC(beta_hat) = (X'X)^(-1) * X' * Omega_hat * X * (X'X)^(-1)

where Omega_hat is a diagonal matrix of weighted squared residuals:

    omega_i = f(residual_i, h_i)

The key insight: instead of assuming a single constant variance sigma², we estimate the covariance by allowing residuals to have different variances.

### 2.2 HC Variants

The choice of weight function omega_i defines different estimators:

| Variant | Formula | Interpretation |
|---------|---------|---|
| **HC0** | omega_i = residual_i² | Basic White estimator (biased in finite samples) |
| **HC1** | omega_i = (n/(n-k)) * residual_i² | Degrees-of-freedom correction |
| **HC2** | omega_i = residual_i² / (1 - h_i) | Leverage-adjusted (minimal bias if errors IID) |
| **HC3** | omega_i = residual_i² / (1-h_i)² | Jackknife approximation |

where h_i is the ith diagonal element of the hat matrix H = X(X'X)^(-1)X'.

**Why HC2?** Under the null hypothesis (true homoskedasticity with IID errors), HC2 is unbiased, whereas HC1 and HC3 suffer from bias. This property makes HC2 the natural choice for establishing a reference distribution.

### 2.3 Implementation

In practice, the sandwich formula is computed via:

1. Fit OLS: beta_hat = (X'X)^(-1)*X'*y
2. Compute residuals: u_hat = y - X*beta_hat  
3. Compute leverage: h_i = X_i*(X'X)^(-1)*X_i'
4. For HC2: weight each residual by 1/(1-h_i)
5. Assemble sandwich: (X'X)^(-1) [weighted cross-product] (X'X)^(-1)

This is the implementation in `R/10_dgp_and_fits.R`.

---

## 3. The Inferential Score: S_Inf

### 3.1 Definition and Interpretation

The **Inferential Score** quantifies how much wider HC-based standard errors are compared to classical ones:

    S_Inf = (SE_robust / SE_classic) - 1

This measures the relative increase in confidence interval width when using HC-based inference instead of classical inference.

**Interpretation:**
- S_Inf = 0: Robust and classical SEs are identical (no heteroskedasticity detected)
- S_Inf > 0: Robust SE is larger (indicates heteroskedasticity or misspecification)
- S_Inf < 0: Robust SE is smaller (rare; suggests anomalies)

### 3.2 Scaled Form: S_Inf*

For hypothesis testing, scaling by sqrt(N) is explored:

    S_Inf* = sqrt(N) * S_Inf

Under the null of true homoskedasticity with HC2, S_Inf* is approximately standard normal, enabling hypothesis testing:
- Test statistic: T = S_Inf*
- Rejection region: |T| > z_critical (e.g., > 1.645 for one-sided alpha=0.05)

---

## 4. Null Calibration: Establishing the Reference Distribution

### 4.1 The Null DGP

To explore the null distribution, completely independent data is generated:

    X_i ~ Normal(0, 1)
    Y_i ~ Normal(0, 1)
    X and Y independent

Despite this independence, the model Y ~ 1 + X is fit. Under repeated sampling:
- The estimated slope beta_X averages to zero (unbiased)
- Classical SE reflects sampling variability
- HC SE varies around the same baseline (under correct specification)

Why this is crucial: The null distribution tells us the "noise floor"—how much inflation the HC estimator introduces when inference is actually valid.

### 4.2 Scaling Law Discovery

Empirically, observations suggest the 95th percentile of S_Inf follows a power law:

    Q_0.95(S_Inf) ≈ c / N^alpha

Taking logarithms:

    log(Q_0.95) ≈ log(c) - alpha * log(N)

A regression of log(Q_0.95) on log(N) yields alpha ≈ 0.5, which is theoretically plausible since ratio-based statistics generally scale as O(1/sqrt(N)).

With alpha = 0.5:

    sqrt(N) * Q_0.95(S_Inf) ≈ c  (approximately invariant in N)

### 4.3 The Calibration Constant

Using HC2 across N in {50, 100, 200, ..., 102400}, preliminary estimates suggest:

    c ≈ 1.64

This closely matches the standard normal 95th percentile z_0.95 ≈ 1.645, providing external validation.

**Interpretation**: Under the null (homoskedasticity), the 95th percentile of the scaled inferential score equals the standard normal critical value. This justifies treating S_Inf* = sqrt(N) * S_Inf as a test statistic with approximately standard normal null distribution.

---

## 5. The Reliability Score: S_Rel

### 5.1 Motivation: The N-Dependence Problem

When heteroskedasticity is present, the raw inferential score exhibits strong N-dependence:

    S_Inf* = sqrt(N) * S_Inf grows with sqrt(N)

This means larger samples lead to larger test statistics, even if the *magnitude* of heteroskedasticity is constant. This N-dependence makes it difficult to create a universal reference table.

**Example**: A dataset with fixed heteroskedasticity might appear "significant" at N=1000 but "insignificant" at N=100, purely due to sample size.

### 5.2 Finite-Sample Correction

Under the null (homoskedasticity), the raw ratio of standard errors exhibits predictable finite-sample bias:

    E[SE_classic / SE_HC2] ≈ 1 + (c_emp / sqrt(N))

where c_emp ≈ 0.25 for simple regression.

To correct for this, we define the **Reliability Score**:

    S_Rel = (SE_classic / SE_HC2) - (c / sqrt(N))

where c is the empirically estimated correction factor.

### 5.3 N-Invariance Property

The key advantage: S_Rel is designed to be approximately **invariant across N**.

Across the heteroskedastic simulation grid, initial analysis examines whether:
- The spread of S_Rel across different N values is much smaller than the spread of raw ratios
- Regression of S_Rel on log(N) yields slopes close to zero
- By mixed-effects modeling, N-dependence is statistically negligible

**Consequence**: The same S_Rel value has approximately the same interpretation (expected coverage loss) regardless of N.

### 5.4 Empirical Estimation of c

To estimate c, the approach involves:

1. Subset heteroskedastic simulations to the null case (lambda = 0)
2. Compute raw ratio: ratio = SE_classic / SE_HC2
3. Regress ratio on 1/sqrt(N):

       ratio = intercept + c * (1/sqrt(N)) + error

4. Report the estimated slope: c_hat ≈ 0.25

This ensures the correction factor is data-driven, not based on theoretical assumptions.

---

## 6. The Reliability Score's Role in Future Validation

The Reliability Score S_Rel provides a sample-size-independent measure of heteroskedasticity severity. By correcting for finite-sample bias with an empirically-derived factor c ≈ 0.25, the reliability score achieves approximate N-invariance.

### 6.1 Purpose and Construction

The reliability score is constructed as:

    S_Rel = (SE_classic / SE_HC2) - (c / sqrt(N))

where c is estimated empirically by regressing the raw ratio on 1/sqrt(N) under the null hypothesis (lambda = 0).

### 6.2 Key Properties

- **N-Invariant**: The correction term c/sqrt(N) removes the primary source of finite-sample drift
- **Empirically Calibrated**: The constant c is derived from simulation under homoskedasticity, not from theory
- **Conservative for Complex Models**: Sensitivity analysis shows the correction is slightly over-aggressive for models with additional covariates, making it a safe choice across contexts

### 6.3 Future Application

The reliability score is designed to serve as a diagnostic metric in future validation analyses. The heteroskedastic simulation stage will verify that S_Rel achieves approximate N-invariance while remaining responsive to heteroskedasticity severity, eventually enabling a practical diagnostic tool for practitioners.
