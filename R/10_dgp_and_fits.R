## ============================================================
## 10_dgp_and_fits.R
##
## Purpose:
##   - Data Generating Process (DGP) for homoskedastic null data
##   - Optimized matrix-based simulation for null calibration
##   - Core OLS + HC engine (fit_ols_hc)
##
## Model Assumptions:
##   - Linear model y = X beta + epsilon
##   - X must include an intercept (column of 1s)
##   - n > k and X must be full rank (unless check_rank=FALSE)
## ============================================================

#' Simulate homoskedastic data under the true null
#'
#' Generates independent normal data for X and Y.
#'
#' @param N Sample size.
#' @return A data.frame with columns 'x' and 'y'.
simulate_homoskedastic_dgp <- function(N) {
  data.frame(
    x = rnorm(N),
    y = rnorm(N)
  )
}

#' Simulate heteroskedastic data (Quadratic Variance)
#'
#' Generates data where the error variance is a quadratic function of X.
#' Formula: sigma^2 = sigma0^2 + lambda * x^2
#'
#' @param N Sample size
#' @param beta_x True coefficient for x
#' @param hetero_strength Strength of heteroskedasticity (lambda)
#' @param sigma0 Baseline noise level
#' @return A data.frame with columns 'x' and 'y'
simulate_heteroskedastic_dgp <- function(N, beta_x = 0.5, hetero_strength = 1, sigma0 = 1.0) {
  x <- rnorm(N)
  # Variance is a function of x
  error_variance <- sigma0^2 + hetero_strength * x^2
  u <- rnorm(N, mean = 0, sd = sqrt(error_variance))
  y <- 1 + beta_x * x + u
  data.frame(x = x, y = y)
}

#' Core OLS and HC Estimator Engine
#'
#' Computes OLS coefficients, classic standard errors, and HC0-HC3
#' robust standard errors using optimized matrix algebra.
#'
#' Implementation details:
#'   1. Cholesky decomposition: R = chol(X'X)
#'   2. Inversion: (X'X)^{-1} = chol2inv(R)
#'   3. Coefficients: beta = (X'X)^{-1} X'y via backsolve/forwardsolve
#'   4. HC covariances: sandwich form (X'X)^{-1} X' Omega X (X'X)^{-1}
#'
#' Validated against lm() + vcovHC() in 00_HC_estimators_validation.Rmd.
#'
#' @param X Numeric matrix (n x k). Must include intercept column.
#' @param y Numeric vector (n).
#' @param hc_types Character vector of HC types to compute ("HC0", "HC1", "HC2", "HC3").
#' @param check_rank Logical. If TRUE (default), checks for rank deficiency and
#'   stops with an informative error if X is not full rank. If FALSE, proceeds
#'   without checking (use only when rank is guaranteed by construction).
#'
#' @return A list containing:
#'   - beta: coefficient vector
#'   - vcov_classic: classic (homoskedastic) covariance matrix
#'   - se_classic: classic standard errors
#'   - vcov_robust: named list of HC covariance matrices
#'   - se_robust: named list of HC standard error vectors
#'   - residuals: OLS residuals
#'   - leverage: diagonal of hat matrix (h_ii)
fit_ols_hc <- function(X, y, hc_types = c("HC0", "HC1", "HC2", "HC3"), check_rank = TRUE) {
  n <- nrow(X)
  k <- ncol(X)
  
  if (check_rank) {
    if (qr(X)$rank < k) {
      stop("X is rank-deficient; OLS/HC not defined.")
    }
  }
  
  # OLS core via Cholesky
  XtX <- crossprod(X)         # k x k SPD
  R   <- chol(XtX)            # upper triangular, XtX = R'R
  XtX_inv <- chol2inv(R)      # (X'X)^(-1)
  
  xty  <- crossprod(X, y)     # X'y
  beta <- backsolve(R, forwardsolve(t(R), xty))
  
  resid <- as.vector(y - X %*% beta)
  sigma2 <- sum(resid^2) / (n - k)
  vcov_classic <- sigma2 * XtX_inv
  se_classic   <- sqrt(diag(vcov_classic))
  
  # leverage
  X_XtXinv <- X %*% XtX_inv
  h <- rowSums(X_XtXinv * X)
  
  # robust covariances
  vcov_robust_list <- list()
  se_robust_list   <- list()
  
  for (hc in hc_types) {
    u_sq <- switch(
      hc,
      "HC0" = resid^2,
      "HC1" = (n / (n - k)) * resid^2,
      "HC2" = resid^2 / (1 - h),
      "HC3" = resid^2 / ((1 - h)^2),
      stop(paste("Unsupported HC type:", hc))
    )
    Xw   <- X * sqrt(u_sq)
    meat <- crossprod(Xw)
    vcov_hc <- XtX_inv %*% meat %*% XtX_inv
    vcov_robust_list[[hc]] <- vcov_hc
    se_robust_list[[hc]]   <- sqrt(diag(vcov_hc))
  }
  
  list(
    beta         = as.vector(beta),
    vcov_classic = vcov_classic,
    se_classic   = se_classic,
    vcov_robust  = vcov_robust_list,
    se_robust    = se_robust_list,
    residuals    = resid,
    leverage     = h
  )
}

#' Optimized Null Simulation using Matrix Algebra
#'
#' Replicates the logic from `ratio hc*.ipynb` notebooks for fast simulation.
#' This function generates data, fits a model, and computes metrics for all
#' specified HC types in a single pass using efficient matrix operations,
#' avoiding the overhead of `lm()` and `sandwich()`.
#'
#' @param n_obs The sample size (N).
#' @param n_sims Number of simulations to run (default 1).
#' @param hc_types A character vector of HC estimators to compute,
#'   e.g., `c("HC1", "HC2", "HC3")`.
#'
#' @return A `data.table` with one row per HC type per simulation, containing:
#'   - sim_id: simulation index (1 to n_sims)
#'   - hc_type: HC type ("HC1", "HC2", "HC3", etc.)
#'   - sr_inf: unscaled inferential score (se_robust / se_classic - 1)
#'   - sr_ratio: raw ratio (se_classic / se_robust)
run_null_simulation_fast <- function(n_obs, n_sims = 1, hc_types = c("HC1", "HC2", "HC3")) {
  
  n_types <- length(hc_types)
  total_rows <- n_sims * n_types
  
  # Pre-allocate vectors for speed
  res_hc_type <- character(total_rows)
  res_sr_inf <- numeric(total_rows)
  res_sr_ratio <- numeric(total_rows)
  res_sim_id <- integer(total_rows)
  
  idx_counter <- 1L
  k <- 2L # Intercept + Slope
  df <- n_obs - k
  
  for (s in seq_len(n_sims)) {
    # 1. DGP - Use canonical function
    d <- simulate_homoskedastic_dgp(n_obs)
    X <- cbind(1, d$x)
    y <- d$y
    
    # 2. OLS Core (Inlined for speed)
    # Deterministic failure: if chol() fails, let it error (per AGENTS.MD)
    XtX <- crossprod(X)
    R <- chol(XtX)                # k x k upper triangular, errors if not SPD
    XtX_inv <- chol2inv(R)
    
    beta <- XtX_inv %*% crossprod(X, y)
    resid <- as.vector(y - X %*% beta)
    sigma2 <- sum(resid^2) / df
    
    # Classic SE for slope (index 2)
    se_classic_x <- sqrt(sigma2 * XtX_inv[2, 2])
    
    # Leverage for HC
    # M = X %*% XtX_inv. We need M[,2] for slope variance and rows for leverage.
    M <- X %*% XtX_inv
    h <- rowSums(M * X)
    
    # Pre-compute M[,2]^2 for robust variance
    M2_sq <- M[, 2]^2
    
    for (h_type in hc_types) {
      # Compute u_sq based on HC type
      if (h_type == "HC0") {
        u_sq <- resid^2
      } else if (h_type == "HC1") {
        u_sq <- (n_obs / df) * resid^2
      } else if (h_type == "HC2") {
        u_sq <- resid^2 / (1 - h)
      } else if (h_type == "HC3") {
        u_sq <- resid^2 / ((1 - h)^2)
      } else {
        stop("Unsupported HC type")
      }
      
      # Robust SE for slope: sqrt( sum( M[,2]^2 * u_sq ) )
      se_robust_x <- sqrt(sum(M2_sq * u_sq))
      
      # Store results
      res_hc_type[idx_counter] <- h_type
      res_sr_inf[idx_counter] <- (se_robust_x / se_classic_x) - 1
      res_sr_ratio[idx_counter] <- se_classic_x / se_robust_x
      res_sim_id[idx_counter] <- s
      
      idx_counter <- idx_counter + 1L
    }
  }
  
  data.table(
    sim_id = res_sim_id,
    hc_type = res_hc_type,
    sr_inf = res_sr_inf,
    sr_ratio = res_sr_ratio
  )
}

#' Optimized Heteroskedastic Simulation using Matrix Algebra
#'
#' Generates heteroskedastic data, fits OLS, and computes coverage and scores.
#' Coverage is computed for the CLASSICAL confidence interval (using se_classic),
#' measuring how much classical inference "breaks" under heteroskedasticity.
#'
#' @param N Sample size
#' @param hetero_strength Strength of heteroskedasticity (lambda in quadratic variance)
#' @param hc_type Robust SE type (default "HC2")
#' @param beta_x True coefficient for x (default 0.5)
#' @param sigma0 Baseline homoskedastic noise SD (default 1.0)
#' @param conf_level Confidence level (default 0.95)
#'
#' @return A data.table with:
#'   - se_classic, se_robust: standard errors
#'   - sr_inf: Inferential Score (SE_robust/SE_classic - 1)
#'   - sr_inf_adj: Scaled Inferential Score (sqrt(N) * sr_inf)
#'   - coverage: 1 if true beta_x is in classical CI (using se_classic), 0 otherwise
#'
#' Note: Ratio metrics (sr_ratio, sr_ratio_adj / Reliability Score) are NOT computed
#' here. They are derived post-hoc in Section 5 after establishing N-dependence
#' of the inferential score.
run_hetero_simulation_fast <- function(N, hetero_strength, hc_type = "HC2", beta_x = 0.5, sigma0 = 1.0, conf_level = 0.95) {
  
  # 1. DGP - Use the canonical function
  d <- simulate_heteroskedastic_dgp(N, beta_x = beta_x, hetero_strength = hetero_strength, sigma0 = sigma0)
  
  X <- cbind(1, d$x)
  y <- d$y
  
  # 2. Fit using core engine
  # check_rank=FALSE: DGP guarantees full rank (intercept + single N(0,1) regressor)
  fit <- fit_ols_hc(X, y, hc_types = hc_type, check_rank = FALSE)
  idx <- 2L  # slope for x
  
  beta_hat   <- fit$beta[idx]
  se_classic <- fit$se_classic[idx]
  se_robust  <- fit$se_robust[[hc_type]][idx]
  
  # 3. Compute Metrics
  # Classical CI coverage: measures how much classical inference degrades
  # under heteroskedasticity. This is the coverage gap we care about.
  df <- N - ncol(X)
  alpha <- 1 - conf_level
  crit_val <- qt(1 - alpha / 2, df = df)
  
  # Classical confidence interval (using se_classic, not se_robust)
  ci_lower <- beta_hat - crit_val * se_classic
  ci_upper <- beta_hat + crit_val * se_classic
  coverage <- (beta_x >= ci_lower) & (beta_x <= ci_upper)  # classical CI coverage
  
  # Scores: Only compute inferential scores during simulation
  # Ratio metrics (sr_ratio, sr_ratio_adj) are computed post-hoc in Section 5
  # after establishing N-dependence of the inferential score
  sr_inf       <- compute_sr_inf(se_classic, se_robust)
  sr_inf_adj   <- compute_sr_inf_adj(se_classic, se_robust, N)
  
  data.table(
    se_classic   = se_classic,
    se_robust    = se_robust,
    sr_inf       = sr_inf,
    sr_inf_adj   = sr_inf_adj,
    coverage     = as.integer(coverage)
  )
}

#' Fit OLS model and compute HC confidence intervals
#'
#' Convenience wrapper around fit_ols_hc() for simple regression (y ~ x).
#' Returns classic and robust SEs plus confidence interval bounds.
#'
#' @param data data.frame with columns 'x' and 'y'
#' @param hc_type Robust SE type (default "HC2")
#' @param conf_level Confidence level (default 0.95)
#' @param check_rank Logical. If TRUE (default), checks for rank deficiency.
#'   Set to FALSE only when rank is guaranteed (e.g., in controlled simulations).
#'
#' @return list with:
#'   - se_classic: Classic standard error for slope
#'   - se_robust: Robust standard error for slope
#'   - ci_lower: Lower bound of robust CI for slope
#'   - ci_upper: Upper bound of robust CI for slope
#'   - beta_x: Estimated slope coefficient
fit_ols_get_hc_ci <- function(data, hc_type = "HC2", conf_level = 0.95, check_rank = TRUE) {
  # 1. Setup matrices (simple regression y ~ x with intercept)
  y <- data$y
  X <- cbind(1, data$x) 
  
  # 2. Fit using core engine
  fit <- fit_ols_hc(X, y, hc_types = hc_type, check_rank = check_rank)
  
  idx <- 2L
  beta_x     <- fit$beta[idx]
  se_classic <- fit$se_classic[idx]
  se_robust  <- fit$se_robust[[hc_type]][idx]
  
  n <- nrow(X)
  k <- ncol(X)
  df <- n - k
  alpha <- 1 - conf_level
  crit_val <- qt(1 - alpha / 2, df = df)
  
  ci_lower <- beta_x - crit_val * se_robust
  ci_upper <- beta_x + crit_val * se_robust
  
  list(
    se_classic = as.numeric(se_classic),
    se_robust = as.numeric(se_robust),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    beta_x = as.numeric(beta_x)
  )
}
