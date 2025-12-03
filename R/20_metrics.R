## ============================================================
## 20_metrics.R
## 
## Purpose:
##   - Compute classic and HCk (HC1/HC2/HC3) standard errors
##   - Define core reliability metrics (Standardized Notation):
##       * sr_inf       = se_robust / se_classic - 1
##       * sr_inf_adj   = sqrt(N) * sr_inf
##       * sr_ratio     = se_classic / se_robust
##       * sr_ratio_adj = sr_ratio - c / sqrt(N)  [c estimated empirically]
##   - Provide generic scaling helpers for null-calibration work
##
## Note on sr_ratio_adj:
##   The correction factor c is estimated empirically in Notebook 02
##   (Section 5.3) by regressing sr_ratio on 1/sqrt(N) under the null.
##   The notebook computes sr_ratio_adj post-hoc using the derived c.
##   Functions here that compute sr_ratio_adj require c as a parameter.
## ============================================================

#' Get classic (homoskedastic) standard errors from an lm model
#'
#' @param fit An lm model.
#' @return Named numeric vector of standard errors for all coefficients.
get_se_classic <- function(fit) {
  sqrt(diag(vcov(fit)))
}

#' Get HCk robust standard errors from an lm model
#'
#' @param fit    An lm model.
#' @param hc_type Character, one of "HC0", "HC1", "HC2", "HC3", etc.
#' @return Named numeric vector of robust standard errors.
get_se_hc <- function(fit, hc_type = "HC2") {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' is required for get_se_hc().")
  }
  vcov_hc <- sandwich::vcovHC(fit, type = hc_type)
  sqrt(diag(vcov_hc))
}

#' Helper: extract sample size N from lm fit
#'
#' @param fit An lm model.
#' @return Integer, number of observations used.
get_n_from_fit <- function(fit) {
  nobs(fit)
}

## ============================================================
## 2. Core metrics for a single coefficient
## ============================================================

#' Compute Unscaled Inferential Score (sr_inf)
#'
#' Definition: (SE_robust / SE_classic) - 1
#'
#' @param se_classic Numeric scalar or vector.
#' @param se_robust  Numeric scalar or vector.
#' @return Numeric scalar or vector.
compute_sr_inf <- function(se_classic, se_robust) {
  se_robust / se_classic - 1
}

#' Compute Adjusted Inferential Score (sr_inf_adj)
#'
#' Definition: sqrt(N) * sr_inf
#' This behaves like a test statistic (T-like).
#'
#' @param se_classic Numeric scalar or vector.
#' @param se_robust  Numeric scalar or vector.
#' @param N          Numeric scalar or vector.
#' @return Numeric scalar or vector.
compute_sr_inf_adj <- function(se_classic, se_robust, N) {
  sqrt(N) * compute_sr_inf(se_classic, se_robust)
}

#' Compute Raw Ratio (sr_ratio)
#'
#' Definition: SE_classic / SE_robust
#'
#' @param se_classic Numeric scalar or vector.
#' @param se_robust  Numeric scalar or vector.
#' @return Numeric scalar or vector.
compute_sr_ratio <- function(se_classic, se_robust) {
  se_classic / se_robust
}

#' Compute Adjusted Reliability Score (sr_ratio_adj, S_Rel)
#'
#' Definition: sr_ratio - c / sqrt(N)
#' where c is the empirically estimated finite-sample correction factor.
#'
#' The correction factor c is estimated in Notebook 02 (Section 5.3) by
#' regressing sr_ratio on 1/sqrt(N) under the homoskedastic null.
#' For simple regression with HC2, c ≈ 0.25.
#'
#' @param se_classic Numeric scalar or vector.
#' @param se_robust  Numeric scalar or vector.
#' @param N          Numeric scalar or vector.
#' @param c_correction Numeric scalar. The empirically estimated correction
#'   factor. Must be provided explicitly (no default).
#' @return Numeric scalar or vector.
compute_sr_ratio_adj <- function(se_classic, se_robust, N, c_correction) {
  if (missing(c_correction)) {
    stop("c_correction must be provided. Estimate it empirically in the notebook.")
  }
  compute_sr_ratio(se_classic, se_robust) - c_correction / sqrt(N)
}

#' Compute reliability metrics for a single coefficient
#'
#' Note: This function computes sr_ratio but NOT sr_ratio_adj, because
#' sr_ratio_adj requires the empirically estimated correction factor c
#' which is derived in Notebook 02. Use compute_sr_ratio_adj() separately
#' after estimating c.
#'
#' @param fit        An lm model.
#' @param coef_name  Name of the coefficient to evaluate (default: COEF_OF_INTEREST from config).
#' @param hc_type    Robust SE type, e.g., "HC2".
#' @param n_override Optional integer N if you want to override nobs(fit).
#'
#' @return A one-row data.table with standardized metric names.
compute_se_metrics <- function(fit,
                               coef_name = COEF_OF_INTEREST,
                               hc_type = "HC2",
                               n_override = NULL) {
  
  # Extract N
  N <- if (is.null(n_override)) get_n_from_fit(fit) else n_override
  
  # Classic and HC SEs
  se_classic_all <- get_se_classic(fit)
  se_hc_all      <- get_se_hc(fit, hc_type = hc_type)
  
  if (!coef_name %in% names(se_classic_all)) {
    stop(sprintf("Coefficient '%s' not found in model.", coef_name))
  }
  
  se_classic <- se_classic_all[coef_name]
  se_robust  <- se_hc_all[coef_name]
  
  # Compute standardized metrics (excluding sr_ratio_adj which needs empirical c)
  sr_inf       <- compute_sr_inf(se_classic, se_robust)
  sr_inf_adj   <- compute_sr_inf_adj(se_classic, se_robust, N)
  sr_ratio     <- compute_sr_ratio(se_classic, se_robust)
  
  data.table(
    N            = N,
    coef_name    = coef_name,
    hc_type      = hc_type,
    se_classic   = as.numeric(se_classic),
    se_robust    = as.numeric(se_robust),
    sr_inf       = as.numeric(sr_inf),
    sr_inf_adj   = as.numeric(sr_inf_adj),
    sr_ratio     = as.numeric(sr_ratio)
  )
}

#' Vectorized wrapper for use inside simulations
#'
#' @param fit       An lm object.
#' @param coef_name Name of coefficient (default: COEF_OF_INTEREST from config).
#' @param hc_type   Robust SE type, default "HC2".
#'
#' @return A one-row data.table with N, metrics, etc.
compute_se_metrics_df <- function(fit,
                                  coef_name = COEF_OF_INTEREST,
                                  hc_type = "HC2") {
  compute_se_metrics(fit, coef_name, hc_type = hc_type)
}

## ============================================================
## 4. Generic scaling helpers (for null calibration)
## ============================================================

#' Scale a metric by N^alpha
#'
#' @param metric Numeric vector (e.g., sr_inf).
#' @param N      Numeric vector of sample sizes (same length as metric).
#' @param alpha  Exponent for N.
#'
#' @return Scaled metric = metric * N^alpha
scale_by_n_power <- function(metric, N, alpha) {
  metric * (N^alpha)
}

## ============================================================
## 5. Small helper to compute all metrics from raw SEs and N
## ============================================================

#' Compute metrics directly from SE vectors
#'
#' Note: This function computes sr_ratio but NOT sr_ratio_adj, because
#' sr_ratio_adj requires the empirically estimated correction factor c.
#'
#' @param se_classic Numeric scalar or vector.
#' @param se_robust  Numeric scalar or vector.
#' @param N          Numeric scalar or vector (same length).
#'
#' @return data.frame with N, sr_inf, sr_inf_adj, sr_ratio.
compute_metrics_from_ses <- function(se_classic, se_robust, N) {
  sr_inf       <- compute_sr_inf(se_classic, se_robust)
  sr_inf_adj   <- compute_sr_inf_adj(se_classic, se_robust, N)
  sr_ratio     <- compute_sr_ratio(se_classic, se_robust)
  
  data.frame(
    N            = N,
    se_classic   = se_classic,
    se_robust    = se_robust,
    sr_inf       = sr_inf,
    sr_inf_adj   = sr_inf_adj,
    sr_ratio     = sr_ratio
  )
}

#' Compute quantiles of a metric
#'
#' @param metric_vector Numeric vector of metric values
#' @param probs         Probabilities for quantiles (default c(0.025, 0.25, 0.5, 0.75, 0.975))
#'
#' @return Named numeric vector of quantiles
compute_metric_quantiles <- function(metric_vector, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  quantile(metric_vector, probs = probs, na.rm = TRUE)
}

## ============================================================
## 6. Additional Metrics and Tests
## ============================================================

#' Compute coverage indicator
#'
#' @param ci_lower Lower bound of CI
#' @param ci_upper Upper bound of CI
#' @param true_param True parameter value
#'
#' @return 1 if true_param is within [ci_lower, ci_upper], 0 otherwise
compute_coverage <- function(ci_lower, ci_upper, true_param) {

  as.integer(true_param >= ci_lower & true_param <= ci_upper)
}

## ============================================================
## 7. Post-hoc Ratio Metric Computation
## ============================================================

#' Add ratio metrics to simulation results (post-hoc)
#'
#' Computes sr_ratio from simulation results that contain se_classic and se_robust.
#' Note: sr_ratio_adj must be computed separately in the notebook using the
#' empirically derived c_correction factor.
#'
#' @param sim_results data.table with columns: N, se_classic, se_robust
#'
#' @return data.table with added column: sr_ratio
add_ratio_metrics_to_results <- function(sim_results) {
  sim_results[, sr_ratio := compute_sr_ratio(se_classic, se_robust)]
  sim_results
}

#' Add aggregated ratio metrics to aggregated results (post-hoc)
#'
#' Computes mean and SD of sr_ratio from raw simulation results and merges
#' them into aggregated results.
#'
#' Note: sr_ratio_adj aggregation must be done in the notebook after
#' computing sr_ratio_adj with the empirically derived c_correction.
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#' @param raw_results        data.table with sr_ratio column
#'
#' @return data.table with added mean/sd columns for sr_ratio
add_aggregated_ratio_metrics <- function(aggregated_results, raw_results) {
  ratio_agg <- raw_results[, .(
    mean_sr_ratio = mean(sr_ratio, na.rm = TRUE),
    sd_sr_ratio   = sd(sr_ratio, na.rm = TRUE)
  ), by = .(N, hetero_strength, beta_x, sigma0)]
  
  merge(aggregated_results, ratio_agg, 
        by = c("N", "hetero_strength", "beta_x", "sigma0"),
        all.x = TRUE)
}
