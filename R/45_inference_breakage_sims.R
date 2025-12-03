## ============================================================
## 45_inference_breakage_sims.R
##
## Purpose:
##   - Inference breakage analysis under exponential heteroskedasticity
##   - Dedicated simulations for studying how classical inference degrades
##   - Uses canonical metrics from R/20_metrics.R
## ============================================================

#' Simulate exponential heteroskedastic DGP
#'
#' Generates data with multiplicative heteroskedasticity where the error
#' standard deviation scales exponentially with the predictor:
#'   error_sd = base_sd * (hetero_strength ^ x)
#'
#' @param N Sample size
#' @param beta_x True slope coefficient
#' @param base_sd Baseline error standard deviation
#' @param hetero_strength Exponential scaling factor (>1 for increasing variance)
#'
#' @return data.frame with columns x, y
simulate_heteroskedastic_dgp_exp <- function(N,
                                              beta_x = 0.5,
                                              base_sd = 1.0,
                                              hetero_strength = 2.0) {
  x <- rnorm(N)
  error_sd <- base_sd * (hetero_strength ^ x)
  y <- 1 + beta_x * x + rnorm(N, mean = 0, sd = error_sd)
  data.frame(x = x, y = y)
}

#' Run inference breakage simulation
#'
#' Evaluates how inference quality degrades under exponential heteroskedasticity.
#' Unlike run_hetero_simulation_fast (which computes only robust CI coverage),
#' this function computes BOTH classical and robust CI coverage to directly
#' study the "breakage" of classical inference.
#'
#' @param N_grid Vector of sample sizes
#' @param hetero_strength_grid Vector of exponential heteroskedasticity strengths
#' @param n_sims Number of simulations per cell
#' @param beta_x True slope coefficient
#' @param hc_type HC estimator type (default "HC2")
#' @param conf_level Confidence level for intervals
#' @param verbose Print progress
#'
#' @return data.table with columns:
#'   - N, hetero_strength, sim_id
#'   - se_classic, se_robust
#'   - sr_inf, sr_inf_adj, sr_ratio (canonical metrics)
#'   - coverage_classic: 1 if true beta in classical CI (using se_classic)
#'   - coverage_robust: 1 if true beta in robust CI (using se_robust)
#'
#' Note: sr_ratio_adj is NOT returned because it requires the empirically
#' estimated correction factor c_correction from Notebook 02 Section 5.3.
#' Use compute_sr_ratio_adj() post-hoc after estimating c_correction.
run_inference_breakage_sim <- function(N_grid,
                                        hetero_strength_grid,
                                        n_sims = 1000,
                                        beta_x = 0.5,
                                        hc_type = "HC2",
                                        conf_level = 0.95,
                                        verbose = TRUE) {
  
  # Build parameter grid
  grid <- expand.grid(N = N_grid, hetero_strength = hetero_strength_grid)
  n_cells <- nrow(grid)
  
  if (verbose) {
    cat(sprintf("Running inference breakage simulation:\n"))
    cat(sprintf("  Grid cells: %d\n", n_cells))
    cat(sprintf("  Sims per cell: %d\n", n_sims))
    cat(sprintf("  Total: %d simulations\n\n", n_cells * n_sims))
  }
  
  # Run in parallel across grid cells
  results <- furrr::future_map_dfr(
    seq_len(n_cells),
    function(cell_idx) {
      N <- grid$N[cell_idx]
      hs <- grid$hetero_strength[cell_idx]
      
      # Run n_sims for this cell
      cell_results <- lapply(seq_len(n_sims), function(sim_id) {
        # Generate data
        df <- simulate_heteroskedastic_dgp_exp(
          N = N, beta_x = beta_x, hetero_strength = hs
        )
        
        # Fit model
        fit <- lm(y ~ x, data = df)
        
        # Extract SEs
        se_classic <- sqrt(diag(vcov(fit)))["x"]
        se_robust <- sqrt(diag(sandwich::vcovHC(fit, type = hc_type)))["x"]
        
        # Compute canonical metrics
        # Note: sr_ratio_adj is NOT computed here because it requires the
        # empirically estimated correction factor c_correction, which is
        # derived in Notebook 02 Section 5.3. Use compute_sr_ratio_adj()
        # post-hoc after estimating c_correction.
        sr_inf <- compute_sr_inf(se_classic, se_robust)
        sr_inf_adj <- compute_sr_inf_adj(se_classic, se_robust, N)
        sr_ratio <- compute_sr_ratio(se_classic, se_robust)
        
        # Compute coverage
        beta_hat <- coef(fit)["x"]
        crit_val <- qt(1 - (1 - conf_level) / 2, df = N - 2)
        
        ci_classic_lower <- beta_hat - crit_val * se_classic
        ci_classic_upper <- beta_hat + crit_val * se_classic
        ci_robust_lower <- beta_hat - crit_val * se_robust
        ci_robust_upper <- beta_hat + crit_val * se_robust
        
        coverage_classic <- as.integer(beta_x >= ci_classic_lower & beta_x <= ci_classic_upper)
        coverage_robust <- as.integer(beta_x >= ci_robust_lower & beta_x <= ci_robust_upper)
        
        data.table(
          N = N,
          hetero_strength = hs,
          sim_id = sim_id,
          se_classic = as.numeric(se_classic),
          se_robust = as.numeric(se_robust),
          sr_inf = sr_inf,
          sr_inf_adj = sr_inf_adj,
          sr_ratio = sr_ratio,
          coverage_classic = coverage_classic,
          coverage_robust = coverage_robust
        )
      })
      
      rbindlist(cell_results)
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  
  if (verbose) cat("Inference breakage simulation complete.\n")
  
  results
}

#' Aggregate inference breakage results
#'
#' Computes summary statistics per (N, hetero_strength) cell, including
#' separate coverage gaps for classical and robust CIs.
#'
#' @param results data.table from run_inference_breakage_sim
#'
#' @return data.table with aggregated statistics:
#'   - coverage_classic, coverage_robust: mean coverage rates
#'   - coverage_gap_classic_pct: 95 - 100 * mean(coverage_classic)
#'     Measures how much classical inference has "broken" under heteroskedasticity.
#'   - coverage_gap_robust_pct: 95 - 100 * mean(coverage_robust)
#'     Measures deviation of robust CI from nominal (should be small if HC works).
aggregate_breakage_results <- function(results) {

  aggregated <- results[, .(
    n_sims = .N,
    mean_sr_inf = mean(sr_inf),
    mean_sr_inf_adj = mean(sr_inf_adj),
    mean_sr_ratio = mean(sr_ratio),
    coverage_classic = mean(coverage_classic),
    coverage_robust = mean(coverage_robust),
    coverage_gap_classic_pct = 95 - mean(coverage_classic) * 100,
    coverage_gap_robust_pct = 95 - mean(coverage_robust) * 100
  ), by = .(N, hetero_strength)]
  
  setorder(aggregated, N, hetero_strength)
  
  aggregated
}

#' Plot inference breakage: classical coverage vs heteroskedasticity strength
#'
#' Visualizes how classical CI coverage degrades as heteroskedasticity increases.
#' This directly shows the "breakage" of classical inference.
#'
#' @param aggregated_results data.table from aggregate_breakage_results
#' @param title Plot title
#'
#' @return ggplot object
plot_inference_breakage <- function(aggregated_results,
                                     title = "Classical Coverage Degradation under Exponential Heteroskedasticity") {
  
  p <- ggplot(aggregated_results, aes(x = hetero_strength, y = coverage_classic * 100, 
                                       color = as.factor(N))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      x = "Heteroskedasticity Strength (exponential factor)",
      y = "Classical Coverage (%)",
      color = "Sample Size (N)",
      caption = "Dashed line: nominal 95% coverage"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  p
}
