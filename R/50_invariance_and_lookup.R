## ============================================================
## 50_invariance_and_lookup.R
## 
## Purpose:
##   - Build wide tables (Coverage_Gap × N × Score)
##   - Compute spreads per coverage gap
##   - Run per-coverage Score ~ log(N) regressions
##   - Run global mixed-effects model
##   - Compute universal lookup table (averaged over N)
## ============================================================

#' Compute spread of adjusted score across N for each coverage gap
#'
#' Spread = max(sr_ratio_adj across N) - min(sr_ratio_adj across N)
#'
#' @param wide_table data.table from build_coverage_adjusted_table
#'
#' @return data.table with Coverage_Gap_Pct and Spread
compute_spread_per_coverage <- function(wide_table) {
  
  # Get N column names (all except Coverage_Gap_Pct)
  N_cols <- names(wide_table)[names(wide_table) != "Coverage_Gap_Pct"]
  N_values <- as.numeric(N_cols)
  
  spreads <- vector("numeric", nrow(wide_table))
  
  for (i in 1:nrow(wide_table)) {
    row_values <- as.numeric(wide_table[i, ..N_cols])
    row_values <- row_values[!is.na(row_values)]
    
    if (length(row_values) > 1) {
      spreads[i] <- max(row_values, na.rm = TRUE) - min(row_values, na.rm = TRUE)
    } else {
      spreads[i] <- NA
    }
  }
  
  data.table(
    Coverage_Gap_Pct = wide_table$Coverage_Gap_Pct,
    Spread = spreads
  )
}

#' Run per-coverage regression: Adjusted ~ log(N)
#'
#' For each coverage gap, fit:
#'   SR_Ratio_Adj ~ log(N)
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#'
#' @return data.table with:
#'   - coverage_gap_pct, slope, intercept, p_value, r_squared, sig_flag, n_points
run_per_coverage_regressions <- function(aggregated_results) {
  
  coverage_gaps <- sort(unique(aggregated_results$coverage_gap_pct))
  results <- vector("list", length(coverage_gaps))
  
  for (i in seq_along(coverage_gaps)) {
    gap <- coverage_gaps[i]
    
    # Subset to this coverage gap (may need to round/group)
    gap_rounded <- round(gap, 1)
    subset_gap <- aggregated_results[
      abs(coverage_gap_pct - gap_rounded) < 0.5
    ]
    
    if (nrow(subset_gap) > 2) {
      # Fit regression
      fit <- lm(mean_sr_ratio_adj ~ log(N), data = subset_gap)
      
      slope <- coef(fit)[2]
      intercept <- coef(fit)[1]
      
      # p-value for slope
      summary_fit <- summary(fit)
      if (nrow(summary_fit$coefficients) >= 2) {
        p_val <- summary_fit$coefficients[2, 4]
      } else {
        p_val <- NA
      }
      r_sq <- summary_fit$r.squared
      
      # Significance flag
      sig_flag <- if (!is.na(p_val) && p_val < 0.05) "YES" else "NO"
      
      results[[i]] <- data.table(
        coverage_gap_pct = gap_rounded,
        slope = slope,
        intercept = intercept,
        p_value = p_val,
        r_squared = r_sq,
        significant_n_dependence = sig_flag,
        n_points = nrow(subset_gap)
      )
    }
  }
  
  rbindlist(results, fill = TRUE)
}

#' Fit global mixed-effects model
#'
#' Fits: mean_sr_ratio_adj ~ log(N) + (1 | cov_gap_factor),
#' where cov_gap_factor is a factor built from rounded coverage_gap_pct.
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#'
#' @return list with:
#'   - model: the lmerTest model object
#'   - slope: fixed effect slope for log(N)
#'   - se_slope: standard error of the slope
#'   - p_value: p-value for the slope
#'   - total_effect: estimated effect over the range of N
#'   - summary_text: text summary of key findings
fit_global_mixed_model <- function(aggregated_results) {
  
  # Convert coverage_gap to factor for random intercept
  aggregated_results[, cov_gap_factor := as.factor(round(coverage_gap_pct, 1))]
  
  # Fit mixed model
  model <- lmerTest::lmer(
    mean_sr_ratio_adj ~ log(N) + (1 | cov_gap_factor),
    data = aggregated_results,
    REML = TRUE
  )
  
  # Extract key results
  summary_model <- summary(model)
  
  slope <- summary_model$coefficients[2, 1]
  se_slope <- summary_model$coefficients[2, 2]
  p_val <- summary_model$coefficients[2, 5]
  
  # Compute effect over range of N
  N_min <- min(aggregated_results$N, na.rm = TRUE)
  N_max <- max(aggregated_results$N, na.rm = TRUE)
  
  log_delta <- log(N_max) - log(N_min)
  total_effect <- slope * log_delta
  
  summary_text <- sprintf(
    "Global mixed model:\n  Slope (log(N)): %.6f (SE: %.6f, p: %.4f)\n  Effect over N=[%.0f, %.0f]: %.4f",
    slope, se_slope, p_val, N_min, N_max, total_effect
  )
  
  list(
    model = model,
    slope = slope,
    se_slope = se_slope,
    p_value = p_val,
    total_effect = total_effect,
    summary_text = summary_text
  )
}

#' Build universal lookup table
#'
#' Averages the reliability score (sr_ratio_adj) across all N values for each
#' classical coverage gap level, producing a single "universal" mapping from
#' coverage gap to reliability score.
#'
#' Note: coverage_gap_pct is derived from classical CI coverage (using se_classic).
#' This measures how much classical inference has degraded under heteroskedasticity.
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#'
#' @return data.table with:
#'   - coverage_gap_pct: classical CI coverage gap (95 - actual classical coverage)
#'   - true_coverage: actual classical CI coverage (95 - gap)
#'   - universal_sr_ratio_adj: averaged reliability score
#'   - n_cells: number of cells averaged
build_universal_lookup_table <- function(aggregated_results) {
  
  # Group by coverage gap and average over N
  universal <- aggregated_results[,
    .(
      universal_sr_ratio_adj = mean(mean_sr_ratio_adj, na.rm = TRUE),
      n_cells = .N
    ),
    by = .(coverage_gap_pct)
  ]
  
  # Sort by coverage gap
  setorder(universal, coverage_gap_pct)
  
  # Compute true coverage from gap
  universal[, true_coverage := 95 - coverage_gap_pct]
  
  # Reorder columns
  setcolorder(universal, c("coverage_gap_pct", "true_coverage", "universal_sr_ratio_adj", "n_cells"))
  
  universal
}

#' Run full invariance and lookup analysis
#'
#' Main orchestration function for:
#'  1. Spread computation
#'  2. Per-coverage regressions
#'  3. Global mixed model
#'  4. Universal lookup table
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#'
#' @return list with:
#'   - spread_table
#'   - regression_table
#'   - mixed_model_result
#'   - universal_lookup
run_full_invariance_analysis <- function(aggregated_results) {
  
  cat("Starting invariance and lookup analysis...\n\n")
  
  # 0. Build wide table (needed for spread)
  cat("Building coverage-adjusted wide table...\n")
  wide_table_list <- build_coverage_adjusted_table(aggregated_results)
  
  # 1. Spread computation
  cat("Computing spreads per coverage gap...\n")
  spread_table <- compute_spread_per_coverage(wide_table_list$wide_table)
  
  # 2. Per-coverage regressions
  cat("Running per-coverage regressions...\n")
  regression_table <- run_per_coverage_regressions(aggregated_results)
  
  # 3. Global mixed model
  cat("Fitting global mixed-effects model...\n")
  mixed_result <- fit_global_mixed_model(aggregated_results)
  cat(mixed_result$summary_text, "\n\n")
  
  # 4. Universal lookup table
  cat("Building universal lookup table...\n")
  universal_lookup <- build_universal_lookup_table(aggregated_results)
  
  cat("Invariance analysis complete.\n\n")
  
  list(
    spread_table = spread_table,
    regression_table = regression_table,
    mixed_model_result = mixed_result,
    universal_lookup = universal_lookup
  )
}

#' Prepare data for lookup table output
#'
#' Format the universal lookup table for publication/use
#'
#' @param universal_lookup data.table from build_universal_lookup_table
#'
#' @return data.table with cleaned up names and formatted values
format_universal_lookup <- function(universal_lookup) {
  
  output <- copy(universal_lookup)
  
  # Round for presentation
  output[, coverage_gap_pct := round(coverage_gap_pct, 1)]
  output[, true_coverage := round(true_coverage, 1)]
  output[, universal_sr_ratio_adj := round(universal_sr_ratio_adj, 4)]
  
  setcolorder(output, c("coverage_gap_pct", "true_coverage", "universal_sr_ratio_adj"))
  
  # Drop n_cells for publication
  output[, n_cells := NULL]
  
  output
}
