## ============================================================
## 30_null_calibration.R
##
## Purpose: OPTIMIZED MODULAR NULL CALIBRATION WORKFLOW
##
## Each function handles ONE step of the analysis:
##   Step 1: run_fast_null_calibration() -> run high-speed, matrix-based simulations
##   Step 2: apply_scaling_and_save()    -> scale and select HC (sr_inf and sr_ratio)
##   Step 3: select_best_hc()            -> final HC selection
##
## ============================================================

#' Step 1: Run fast null calibration simulations
#'
#' @param N_grid Vector of sample sizes
#' @param n_sims_per_N Number of simulations per sample size
#' @param hc_types Vector of HC types to test
#' @param verbose Print progress messages
#' @param output_dir Directory to save results
#' @return A data.table with the simulation results
run_fast_null_calibration <- function(N_grid,
                                      n_sims_per_N = 10000,
                                      hc_types = c("HC1", "HC2", "HC3"),
                                      verbose = TRUE,
                                      output_dir = "results") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (verbose) cat("Step 1: Running fast null calibration...\n")

  # Initialize results list
  results_all <- vector("list", length(N_grid))
  
  # Process each N sequentially to allow for progress updates and checks
  for (i in seq_along(N_grid)) {
    N <- N_grid[i]
    
    if (verbose) {
      cat(sprintf("\nProcessing N = %d (%d/%d)...\n", N, i, length(N_grid)))
    }
    
    # Define batches to avoid massive future overhead
    # Batch size of 500 is a reasonable balance
    batch_size <- 500
    n_batches <- ceiling(n_sims_per_N / batch_size)
    
    if (verbose) cat(sprintf("  - Running %d simulations in %d batches...\n", n_sims_per_N, n_batches))
    
    # Run batches in parallel
    batch_results <- furrr::future_map(
      1:n_batches,
      function(b) {
        # Determine how many sims in this batch
        n_in_batch <- if(b == n_batches) n_sims_per_N - (n_batches-1)*batch_size else batch_size
        
        # Run n_in_batch simulations sequentially within the worker
        # This reduces the overhead of spawning futures for every single simulation
        # Use vectorized n_sims for performance (avoids creating thousands of data.tables)
        res <- run_null_simulation_fast(N, n_sims = n_in_batch, hc_types = hc_types)
        
        # Adjust sim_id to be global across batches
        res[, sim_id := sim_id + (b-1)*batch_size]
        
        res
      },
      .options = furrr::furrr_options(seed = TRUE, packages = c("data.table", "sandwich"))
    )
    
    # Combine batch results for this N
    N_results <- rbindlist(batch_results)
    N_results[, N := N]
    
    # CHECK: Deterministic failure on NA/Inf
    # We strictly enforce that no results can be NA or Inf
    if (anyNA(N_results)) {
      stop(sprintf("Deterministic Failure: NA values detected in simulation results for N=%d", N))
    }
    
    # Check for Inf values in numeric columns
    numeric_cols <- names(N_results)[sapply(N_results, is.numeric)]
    for (col in numeric_cols) {
      if (any(is.infinite(N_results[[col]]))) {
        stop(sprintf("Deterministic Failure: Inf values detected in column '%s' for N=%d", col, N))
      }
    }
    
    if (verbose) cat(sprintf("  - Completed N=%d. Checks passed.\n", N))
    
    results_all[[i]] <- N_results
  }

  # Combine results
  scores_long <- rbindlist(results_all)


  # Save results
  saveRDS(scores_long, file.path(output_dir, "scores_long.rds"))

  if (verbose) cat(sprintf("\n✓ Step 1 complete: results saved to %s\n", file.path(output_dir, "scores_long.rds")))

  return(scores_long)
}


#' Step 2: Apply scaling analysis
#'
#' Computes Q95 for sr_inf and sr_ratio deviation, tests scaling across alpha grid.
#'
#' @param scores_long  data.table with scores
#' @param alpha_grid   Vector of exponents to test
#' @param verbose      Print progress
#' @param output_dir   Directory to save results
#'
#' @return list with summaries
apply_scaling_and_save <- function(scores_long,
                                    alpha_grid = seq(0.3, 0.7, by = 0.05),
                                    verbose = TRUE,
                                    output_dir = "results") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (verbose) cat("Step 2a: Computing Q95 per N per HC...\n")

  # Compute Q95 per N and HC type for sr_inf and sr_ratio deviation
  # We look at abs(sr_ratio - 1) for scaling
  scores_long[, sr_ratio_dev := abs(sr_ratio - 1)]

  q95_summary <- scores_long[, .(
    Q95_SR_Inf = quantile(sr_inf, probs = 0.95, na.rm = TRUE),
    Q95_SR_Ratio = quantile(sr_ratio_dev, probs = 0.95, na.rm = TRUE)
  ), by = .(N, hc_type)]

  if (verbose) cat(sprintf("  Q95 computed: %d rows\n\n", nrow(q95_summary)))

  # Test scaling across alpha grid
  if (verbose) cat("Step 2b: Testing alpha grid...\n")

  hc_types <- unique(q95_summary$hc_type)
  best_results <- vector("list", 0)
  all_alphas <- vector("list", 0)

  for (hc in hc_types) {
    q95_hc <- q95_summary[hc_type == hc]

    slopes_inf <- numeric(length(alpha_grid))
    slopes_ratio <- numeric(length(alpha_grid))

    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]

      # Scale Q95 by N^alpha (Multiplication to cancel N^-0.5 decay)
      q95_test <- copy(q95_hc)
      q95_test[, Q95_SR_Inf_Scaled := Q95_SR_Inf * (N^alpha)]
      q95_test[, Q95_SR_Ratio_Scaled := Q95_SR_Ratio * (N^alpha)]

      # Fit log(Q95_Scaled) ~ log(N)
      fit_inf <- lm(log(Q95_SR_Inf_Scaled) ~ log(N), data = q95_test)
      fit_ratio <- lm(log(Q95_SR_Ratio_Scaled) ~ log(N), data = q95_test)

      slopes_inf[i] <- coef(fit_inf)[2]
      slopes_ratio[i] <- coef(fit_ratio)[2]

      # Store for visualization
      q95_test[, alpha := alpha]
      q95_test[, slope_inf := slopes_inf[i]]
      q95_test[, slope_ratio := slopes_ratio[i]]
      all_alphas[[length(all_alphas) + 1]] <- q95_test
    }

    # Find best alpha for Inf Score
    best_idx_inf <- which.min(abs(slopes_inf))
    best_alpha_inf <- alpha_grid[best_idx_inf]

    # Find best alpha for Ratio
    best_idx_ratio <- which.min(abs(slopes_ratio))
    best_alpha_ratio <- alpha_grid[best_idx_ratio]

    # Compute stability range at best alpha (Inf Score)
    q95_best <- copy(q95_hc)
    q95_best[, Q95_SR_Inf_Scaled := Q95_SR_Inf * (N^best_alpha_inf)]
    stability_range <- max(q95_best$Q95_SR_Inf_Scaled, na.rm = TRUE) -
                      min(q95_best$Q95_SR_Inf_Scaled, na.rm = TRUE)

    best_results[[length(best_results) + 1]] <- data.table(
      hc_type = hc,
      best_alpha_inf = best_alpha_inf,
      best_slope_inf = slopes_inf[best_idx_inf],
      stability_range_inf = stability_range,
      best_alpha_ratio = best_alpha_ratio,
      best_slope_ratio = slopes_ratio[best_idx_ratio]
    )

    if (verbose) {
      cat(sprintf("  %s: Inf alpha=%.2f (stab=%.4f), Ratio alpha=%.2f\n",
                  hc, best_alpha_inf, stability_range, best_alpha_ratio))
    }
  }

  best_results_dt <- rbindlist(best_results)
  all_alphas_dt <- rbindlist(all_alphas)

  # Save all results
  saveRDS(q95_summary, file.path(output_dir, "q95_summary.rds"))
  saveRDS(best_results_dt, file.path(output_dir, "scaling_results.rds"))
  saveRDS(all_alphas_dt, file.path(output_dir, "all_alphas.rds"))

  if (verbose) {
    cat(sprintf("\n✓ Step 2 complete:\n"))
    cat(sprintf("  Q95 summary: %d rows\n", nrow(q95_summary)))
    cat(sprintf("  Scaling results: %d HC types\n", nrow(best_results_dt)))
  }

  list(
    q95_summary = q95_summary,
    scaling_results = best_results_dt,
    all_alphas = all_alphas_dt
  )
}

#' Step 3: Select best HC type
#'
#' Ranks HC types by stability_range of Inf Score.
#'
#' @param scaling_results Output from apply_scaling_and_save()$scaling_results
#' @param verbose        Print results
#' @param output_dir     Directory to save results
#'
#' @return Selected HC type (character)
select_best_hc <- function(scaling_results,
                           verbose = TRUE,
                           output_dir = "results") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  setkey(scaling_results, stability_range_inf)
  best_hc <- scaling_results[1, hc_type]

  if (verbose) {
    cat("\n", strrep("=", 80), "\n", sep = "")
    cat("STEP 3: HC TYPE SELECTION\n")
    cat(strrep("=", 80), "\n", sep = "")
    cat("\nRanking by stability_range_inf (smaller = better):\n\n")

    for (i in seq_len(nrow(scaling_results))) {
      row <- scaling_results[i]
      
      if (i == 1) rank_str <- "✓ BEST"
      else if (i == 2) rank_str <- "  2nd "
      else if (i == 3) rank_str <- "  3rd "
      else rank_str <- sprintf("  %dth ", i)
      
      cat(sprintf("%s %s: stability=%.6f, alpha_inf=%.2f, alpha_ratio=%.2f\n",
                  rank_str, row$hc_type, row$stability_range_inf,
                  row$best_alpha_inf, row$best_alpha_ratio))
    }

    cat(sprintf("\n→ SELECTED: %s\n", best_hc))
  }

  # Save results
  saveRDS(best_hc, file.path(output_dir, "best_hc.rds"))
  write.csv(scaling_results, file.path(output_dir, "hc_ranking.csv"), row.names = FALSE)

  best_hc
}
