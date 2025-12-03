## ============================================================
## 40_hetero_sims.R
## 
## Purpose:
##   - Run heteroskedastic simulations over grid of N and hetero_strength
##   - Compute coverage + all metrics for each simulation
##   - Parallel execution for efficiency
## ============================================================

#' Run heteroskedastic simulation grid
#'
#' For each (N, hetero_strength, beta_x, sigma0) combination, run n_sims_per_cell simulations
#' and compute metrics. Grid cells are processed in parallel via `furrr::future_map_dfr()`
#' using the global `future::plan()` configuration set in `R/00_config.R`.
#' Simulations within each cell executed serially.
#' Random seeds derived from global seed via `furrr_options(seed = TRUE)`.
#'
#' @param N_grid                Vector of sample sizes
#' @param hetero_strength_grid  Vector of heteroskedasticity strengths
#' @param beta_x_grid           Vector of slope coefficients
#' @param sigma0_grid           Vector of baseline noise SDs
#' @param n_sims_per_cell       Number of sims per cell
#' @param hc_type              Robust SE type (default "HC2")
#' @param conf_level           Confidence level
#' @param verbose              Print progress
#'
#' @return data.table with columns:
#'   - N, hetero_strength, beta_x, sigma0, sim_id
#'   - se_classic, se_robust, sr_inf, sr_inf_adj
#'   - coverage (classical CI coverage)
#'
#' Note: Ratio metrics (sr_ratio, sr_ratio_adj) are NOT computed during simulation.
#' They are derived post-hoc in Section 5 after establishing N-dependence.
run_hetero_sim_grid <- function(N_grid,
                                hetero_strength_grid,
                                beta_x_grid,
                                sigma0_grid,
                                n_sims_per_cell = 500,
                                hc_type = "HC2",
                                conf_level = 0.95,
                                verbose = TRUE) {
  
  # Create parameter grid
  param_grid <- expand.grid(
    N = N_grid, 
    hetero_strength = hetero_strength_grid,
    beta_x = beta_x_grid,
    sigma0 = sigma0_grid
  )
  
  total_cells <- nrow(param_grid)
  total_sims <- total_cells * n_sims_per_cell
  
  if (verbose) {
    cat(sprintf("Running heteroskedastic simulation grid:\n"))
    cat(sprintf("  N values: %s\n", paste(N_grid, collapse = ", ")))
    cat(sprintf("  Hetero strengths: %s\n", paste(range(hetero_strength_grid), collapse = " to ")))
    cat(sprintf("  Beta X values: %s\n", paste(beta_x_grid, collapse = ", ")))
    cat(sprintf("  Sigma0 values: %s\n", paste(sigma0_grid, collapse = ", ")))
    cat(sprintf("  Sims per cell: %d\n", n_sims_per_cell))
    cat(sprintf("  Total grid cells: %d\n", total_cells))
    cat(sprintf("  Total simulations: %d\n\n", total_sims))
  }
  
  # Parallel execution: process each grid cell in parallel
  # Within each worker, run n_sims_per_cell serially
  results_all <- furrr::future_map_dfr(
    seq_len(nrow(param_grid)),
    function(i) {
      N       <- param_grid$N[i]
      lambda  <- param_grid$hetero_strength[i]
      beta_x  <- param_grid$beta_x[i]
      sigma0  <- param_grid$sigma0[i]
      
      if (verbose && (i %% max(1, floor(total_cells / 20)) == 0)) {
         # Simple progress indicator (approximate, as workers run in parallel)
         # cat(sprintf("Processing cell %d/%d...\n", i, total_cells))
      }

      # Run simulations for this cell serially
      sim_list <- vector("list", n_sims_per_cell)
      for (k in seq_len(n_sims_per_cell)) {
        res <- run_hetero_simulation_fast(
          N               = N,
          hetero_strength = lambda,
          hc_type        = hc_type,
          beta_x         = beta_x,
          sigma0         = sigma0,
          conf_level     = conf_level
        )
        res[, sim_id := k]
        sim_list[[k]] <- res
      }
      
      cell_dt <- rbindlist(sim_list)
      
      # Add cell identifiers
      cell_dt[, `:=`(
        N               = N,
        hetero_strength = lambda,
        beta_x          = beta_x,
        sigma0          = sigma0
      )]
      
      # DETERMINISTIC FAILURE CHECKS (per AGENTS.MD Section 6.4)
      # Any NA or Inf indicates a bug or numerical instability. Fail loudly.
      if (anyNA(cell_dt)) {
        stop(sprintf("Deterministic Failure: NA values detected in cell N=%d, lambda=%.2f, beta=%.2f, sigma0=%.2f", 
                     N, lambda, beta_x, sigma0))
      }
      
      numeric_cols <- names(cell_dt)[sapply(cell_dt, is.numeric)]
      for (col in numeric_cols) {
        if (any(is.infinite(cell_dt[[col]]))) {
          stop(sprintf("Deterministic Failure: Inf values detected in column '%s' for cell N=%d, lambda=%.2f, beta=%.2f, sigma0=%.2f", 
                       col, N, lambda, beta_x, sigma0))
        }
      }
      
      cell_dt
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  
  if (verbose) {
    cat(sprintf("\n✓ Heteroskedastic grid simulations complete.\n"))
    cat(sprintf("  Total simulations completed: %d\n", nrow(results_all)))
  }
  
  results_all
}

#' Aggregate heteroskedastic simulation results
#'
#' For each (N, hetero_strength, beta_x, sigma0) cell, computes:
#'   - Mean and SD of se_classic, se_robust
#'   - Mean and SD of sr_inf, sr_inf_adj
#'   - Mean coverage (for classical CI using se_classic)
#'   - coverage_gap_pct = 95 - 100 * mean(coverage)
#'
#' Note: The "coverage" column from run_hetero_simulation_fast measures
#' whether the CLASSICAL confidence interval (using se_classic) covers the true
#' parameter. The coverage_gap_pct therefore measures how much classical
#' inference has "broken" under heteroskedasticity.
#'
#' Note: Ratio metrics (sr_ratio, sr_ratio_adj) are added post-hoc via
#' add_ratio_metrics_to_results() after this aggregation.
#'
#' @param hetero_sim_results data.table from run_hetero_sim_grid
#'
#' @return data.table with aggregated statistics per cell
aggregate_hetero_sims <- function(hetero_sim_results) {
  
  # Core metrics from simulation: inferential scores and coverage
  # Note: 'coverage' is for classical CI (se_classic), so coverage_gap_pct
  # measures how much classical inference degrades under heteroskedasticity
  aggregated <- hetero_sim_results[, .(
    n_sims = .N,
    mean_se_classic = mean(se_classic, na.rm = TRUE),
    sd_se_classic = sd(se_classic, na.rm = TRUE),
    mean_se_robust = mean(se_robust, na.rm = TRUE),
    sd_se_robust = sd(se_robust, na.rm = TRUE),
    mean_sr_inf = mean(sr_inf, na.rm = TRUE),
    sd_sr_inf = sd(sr_inf, na.rm = TRUE),
    mean_sr_inf_adj = mean(sr_inf_adj, na.rm = TRUE),
    sd_sr_inf_adj = sd(sr_inf_adj, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    coverage_gap_pct = 95 - mean(coverage, na.rm = TRUE) * 100  # classical CI coverage gap
  ), by = .(N, hetero_strength, beta_x, sigma0)]
  
  # Ratio metrics (sr_ratio, sr_ratio_adj) are added later via add_ratio_metrics_to_results()
  
  setorder(aggregated, N, hetero_strength, beta_x, sigma0)
  aggregated
}

#' Build wide format table: Coverage_Gap_Pct rows, N columns, Reliability Score values
#'
#' Creates the main lookup table: for each classical coverage gap level, shows the
#' corresponding reliability score (sr_ratio_adj) at each sample size N.
#'
#' Note: coverage_gap_pct is derived from classical CI coverage (using se_classic),
#' as computed by run_hetero_simulation_fast and aggregated by aggregate_hetero_sims.
#' This measures how much classical inference has degraded under heteroskedasticity.
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#'
#' @return list with:
#'   - wide_table: data.table (Coverage_Gap_Pct as rows, N as columns, values are sr_ratio_adj)
#'   - coverage_gaps: vector of unique classical coverage gap levels
#'   - N_values: vector of sample sizes
build_coverage_adjusted_table <- function(aggregated_results) {
  
  # Round coverage_gap_pct to nearest integer for cleaner display
  aggregated_results[, coverage_gap_rounded := round(coverage_gap_pct, 0)]
  
  # Get unique coverage gaps (sorted)
  coverage_gaps <- sort(unique(aggregated_results$coverage_gap_rounded))
  N_values <- sort(unique(aggregated_results$N))
  
  # Create matrix
  wide_matrix <- matrix(NA, 
                        nrow = length(coverage_gaps),
                        ncol = length(N_values))
  
  for (i in seq_along(coverage_gaps)) {
    for (j in seq_along(N_values)) {
      subset <- aggregated_results[
        coverage_gap_rounded == coverage_gaps[i] & N == N_values[j]
      ]
      if (nrow(subset) > 0) {
        # Average if multiple rows (shouldn't happen, but just in case)
        wide_matrix[i, j] <- mean(subset$mean_sr_ratio_adj, na.rm = TRUE)
      }
    }
  }
  
  # Convert to data.table
  wide_dt <- as.data.table(wide_matrix)
  setnames(wide_dt, as.character(N_values))
  wide_dt[, Coverage_Gap_Pct := coverage_gaps]
  
  # Reorder columns
  setcolorder(wide_dt, c("Coverage_Gap_Pct", as.character(N_values)))
  
  list(
    wide_table = wide_dt,
    coverage_gaps = coverage_gaps,
    N_values = N_values
  )
}

#' Interpolate to fill missing values in coverage-adjusted table
#'
#' Use linear interpolation across hetero levels for smooth mapping
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#' @param N_val             Specific N to interpolate for
#' @param coverage_gap_seq   Sequence of coverage gaps to interpolate to
#'
#' @return data.table with interpolated results
interpolate_coverage_gap_for_N <- function(aggregated_results,
                                           N_val,
                                           coverage_gap_seq = seq(0, 40, by = 1)) {
  
  subset_N <- aggregated_results[N == N_val]
  subset_N <- subset_N[order(coverage_gap_pct)]
  
  # Linear interpolation
  interp_scores <- approx(
    x = subset_N$coverage_gap_pct,
    y = subset_N$mean_sr_ratio_adj,
    xout = coverage_gap_seq,
    method = "linear",
    rule = 2  # use boundary values if out of range
  )
  
  data.table(
    N = N_val,
    coverage_gap_pct = interp_scores$x,
    sr_ratio_adj_interp = interp_scores$y
  )
}
