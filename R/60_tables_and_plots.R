## ============================================================
## 60_tables_and_plots.R
## 
## Purpose:
##   - Pretty table outputs (kable)
##   - Plots and visualizations
##   - Saving CSV/PNG outputs
## ============================================================

#' Pretty-print a data.table or data.frame as a kable table
#'
#' @param dt           data.table or data.frame
#' @param caption      Table caption
#' @param digits       Number of decimal places
#' @param full_width   Use full width
#'
#' @return kable object (for display in notebook)
print_kable_table <- function(dt, 
                              caption = "",
                              digits = 3,
                              full_width = FALSE) {
  
  # Convert to data.frame if needed
  if (is.data.table(dt)) {
    df <- as.data.frame(dt)
  } else {
    df <- dt
  }
  
  kable(
    df,
    caption = caption,
    digits = digits,
    format = "html"
  )
}

#' Save results to CSV
#'
#' @param dt       data.table or data.frame
#' @param filename Filename (will be saved in results/ directory)
#' @param verbose  Print confirmation message
save_table_csv <- function(dt, filename, verbose = TRUE) {
  
  output_path <- file.path(RESULTS_PATH, filename)
  
  if (is.data.table(dt)) {
    fwrite(dt, output_path)
  } else {
    write.csv(dt, output_path, row.names = FALSE)
  }
  
  if (verbose) {
    cat(sprintf("Saved: %s\n", output_path))
  }
}

#' Save results to RDS (binary format for R)
#'
#' @param obj       R object to save
#' @param filename  Filename (will be saved in results/ directory)
#' @param verbose   Print confirmation message
save_rds <- function(obj, filename, verbose = TRUE) {
  
  output_path <- file.path(RESULTS_PATH, filename)
  saveRDS(obj, output_path)
  
  if (verbose) {
    cat(sprintf("Saved: %s\n", output_path))
  }
}

#' Create a plot of Reliability Score vs True Coverage
#'
#' @param universal_lookup data.table with coverage_gap_pct, true_coverage, universal_sr_ratio_adj
#' @param title            Plot title
#'
#' @return ggplot object
plot_reliability_score_vs_coverage <- function(universal_lookup, 
                                               title = "Universal Reliability Score vs True Coverage") {
  
  p <- ggplot(universal_lookup, 
              aes(x = true_coverage, y = universal_sr_ratio_adj)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2) +
    labs(
      title = title,
      x = "True Coverage (%)",
      y = "Universal Reliability Score",
      caption = "Based on simulations averaged across N"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  p
}

#' Create a plot of Spread vs Coverage Gap
#'
#' @param spread_table data.table with Coverage_Gap_Pct, Spread
#' @param title        Plot title
#'
#' @return ggplot object
plot_spread_vs_coverage <- function(spread_table,
                                    title = "Spread (Range across N) vs Coverage Gap") {
  
  p <- ggplot(spread_table, aes(x = Coverage_Gap_Pct, y = Spread)) +
    geom_line(color = "darkred", size = 1) +
    geom_point(color = "darkred", size = 2) +
    labs(
      title = title,
      x = "Coverage Gap (%)",
      y = "Spread (max - min)",
      caption = "Spread computed as range of reliability scores across N values"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  p
}

#' Create a plot of Reliability Score across N for selected coverage gaps
#'
#' @param aggregated_results data.table from aggregate_hetero_sims
#' @param coverage_gaps_to_plot Vector of coverage gaps to include
#' @param title                 Plot title
#'
#' @return ggplot object
plot_reliability_score_by_n <- function(aggregated_results,
                                        coverage_gaps_to_plot = c(0, 10, 20, 30, 40),
                                        title = "Reliability Score vs Sample Size by Coverage Gap") {
  
  # Subset to selected coverage gaps
  subset_data <- aggregated_results[
    round(coverage_gap_pct, 0) %in% coverage_gaps_to_plot
  ]
  subset_data[, cov_gap_label := paste0("Gap: ", round(coverage_gap_pct, 0), "%")]
  
  p <- ggplot(subset_data, aes(x = N, y = mean_sr_ratio_adj, color = cov_gap_label)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(
      title = title,
      x = "Sample Size (log scale)",
      y = "Mean Reliability Score",
      color = "Coverage Gap",
      caption = "Mean Reliability Score by coverage gap"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  p
}

#' Create a plot of null distribution quantiles
#'
#' @param quantile_summary data.table with Q95_SR_Inf and Q95_SR_Ratio columns,
#'   e.g. the q95_summary output from apply_scaling_and_save().
#' @param hc_types_to_plot Vector of HC types to include
#' @param metric           "sr_inf" or "sr_ratio"
#' @param title            Plot title
#'
#' @return ggplot object
plot_null_quantiles <- function(quantile_summary,
                                hc_types_to_plot = c("HC1", "HC2", "HC3"),
                                metric = "sr_inf",
                                title = "Null Distribution Quantiles (Inferential Score)") {
  
  # Filter data
  # Note: The column names in quantile_summary are Q95_SR_Inf and Q95_SR_Ratio
  # We need to map the input 'metric' to the column name
  
  col_name <- if (metric == "sr_inf") "Q95_SR_Inf" else "Q95_SR_Ratio"
  
  subset_data <- quantile_summary[hc_type %in% hc_types_to_plot]
  
  p <- ggplot(subset_data, aes(x = N, y = .data[[col_name]], color = hc_type)) +
    geom_line(size = 1) +
    scale_x_log10() +
    labs(
      title = title,
      x = "Sample Size (log scale)",
      y = "Quantile Value (Q95)",
      color = "HC Type",
      caption = "Null simulations (homoskedastic)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  p
}

#' Save a ggplot to PNG
#'
#' @param plot      ggplot object
#' @param filename  Filename (will be saved in results/ directory)
#' @param width     Width in inches (default 10)
#' @param height    Height in inches (default 6)
#' @param verbose   Print confirmation message
save_plot_png <- function(plot, filename, width = 10, height = 6, verbose = TRUE) {
  
  output_path <- file.path(RESULTS_PATH, filename)
  
  ggsave(
    filename = output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  if (verbose) {
    cat(sprintf("Saved: %s\n", output_path))
  }
}

#' Print comprehensive summary of all results
#'
#' @param null_results      List returned by apply_scaling_and_save(),
#'   containing scaling_results and related summaries.
#' @param hetero_results    data.table from run_hetero_sim_grid
#' @param aggregated        data.table from aggregate_hetero_sims
#' @param invariance_results List from run_full_invariance_analysis
print_full_summary <- function(null_results,
                               hetero_results,
                               aggregated,
                               invariance_results) {
  
  title <- "FULL ANALYSIS SUMMARY"
  cat("\n")
  cat(strrep("=", nchar(title)), "\n")
  cat(title, "\n")
  cat(strrep("=", nchar(title)), "\n\n")
  
  # Null calibration summary
  cat("1. NULL CALIBRATION RESULTS\n")
  cat(strrep("-", 30), "\n")
  print(null_results$scaling_results)
  cat("\n\n")
  
  # Spread summary
  cat("2. SPREAD TABLE (N-Invariance)\n")
  cat(strrep("-", 30), "\n")
  print(invariance_results$spread_table)
  cat("\n\n")
  
  # Regression summary
  cat("3. PER-COVERAGE REGRESSION RESULTS\n")
  cat(strrep("-", 30), "\n")
  print(invariance_results$regression_table)
  cat("\n\n")
  
  # Mixed model summary
  cat("4. GLOBAL MIXED MODEL RESULTS\n")
  cat(strrep("-", 30), "\n")
  cat(invariance_results$mixed_model_result$summary_text)
  cat("\n\n")
  
  # Universal lookup table
  cat("5. UNIVERSAL LOOKUP TABLE\n")
  cat(strrep("-", 30), "\n")
  print(invariance_results$universal_lookup)
  cat("\n\n")
}

## ============================================================
## Centralized Visualization Functions (from 02_HC2_validation.Rmd)
## ============================================================

#' Plot power curves for heteroskedasticity detection
#'
#' @param power_summary data.table with columns: N, hetero_strength, rejection_rate
#' @param title Plot title
#' @param subtitle Plot subtitle
#'
#' @return ggplot object
plot_power_curves <- function(power_summary,
                               title = "Power of Inferential Score to Detect Heteroskedasticity",
                               subtitle = "Rejection Rate of T_SInf > 1.645 (alpha = 0.05)") {
  
  p <- ggplot(power_summary, aes(x = hetero_strength, y = rejection_rate, color = as.factor(N))) +
    geom_line(size = 1) +
    geom_point() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Heteroskedasticity Strength (lambda)",
      y = "Rejection Rate (Power)",
      color = "Sample Size (N)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  p
}

#' Plot N-dependence of test statistic
#'
#' @param sinf_n_summary data.table with columns: N, hetero_strength, mean_T_sinf
#' @param title Plot title
#' @param subtitle Plot subtitle
#'
#' @return ggplot object
plot_n_dependence <- function(sinf_n_summary,
                               title = "N-Dependence of Test Statistic T_SInf",
                               subtitle = "Mean T_SInf increases with N for fixed misspecification") {
  
  p <- ggplot(sinf_n_summary, aes(x = N, y = mean_T_sinf, color = as.factor(hetero_strength))) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample Size (log scale)",
      y = "Mean T_SInf",
      color = "Lambda"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  p
}

#' Plot N-invariance comparison (Raw Ratio vs Reliability Score)
#'
#' @param metrics_n_summary data.table with columns: N, hetero_strength, mean_sr_ratio, mean_sr_ratio_adj
#'
#' @return arranged ggplot (grid of 2 plots)
plot_n_invariance_comparison <- function(metrics_n_summary) {
  
  # Plot Raw Ratio vs N
  p1 <- ggplot(metrics_n_summary, aes(x = N, y = mean_sr_ratio, color = as.factor(hetero_strength))) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(title = "Raw Ratio vs N", y = "Mean Raw Ratio", x = "N (log)") +
    theme_minimal() + 
    theme(legend.position = "none")
  
  # Plot Reliability Score vs N
  p2 <- ggplot(metrics_n_summary, aes(x = N, y = mean_sr_ratio_adj, color = as.factor(hetero_strength))) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(title = "Reliability Score vs N", y = "Mean Reliability Score", x = "N (log)") +
    theme_minimal() + 
    theme(legend.position = "right")
  
  gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(0.45, 0.55))
}

#' Plot benchmark comparison (S_Inf vs BP/White tests)
#'
#' @param bench_long data.table in long format with columns: N, lambda, Test, Power
#' @param title Plot title
#'
#' @return ggplot object
plot_benchmark_comparison <- function(bench_long,
                                       title = "Power Comparison: S_Inf vs Standard Tests") {
  
  p <- ggplot(bench_long, aes(x = lambda, y = Power, color = Test)) +
    geom_line() +
    facet_wrap(~N) +
    labs(
      title = title,
      x = "Heteroskedasticity Strength (lambda)",
      y = "Rejection Rate (alpha=0.05)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  p
}

#' Plot universal lookup table
#'
#' @param universal_lookup data.table with coverage_gap_pct, universal_sr_ratio_adj
#' @param hc_type Selected HC type for title
#'
#' @return ggplot object
plot_universal_lookup <- function(universal_lookup, hc_type = "HC2") {
  
  p <- ggplot(
    universal_lookup,
    aes(x = coverage_gap_pct, y = universal_sr_ratio_adj)
  ) +
    geom_line(size = 1, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    labs(
      title = sprintf("Universal %s Lookup Table", hc_type),
      x = "Empirical Coverage Gap (%)",
      y = "Reliability Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12)
    )
  
  p
}

#' Plot inferential score distribution
#'
#' @param hetero_results data.table with sr_inf and hetero_strength columns
#' @param lambda_subset Vector of lambda values to include in density plot
#'
#' @return arranged ggplot grid (histogram + density by lambda)
plot_sinf_distribution <- function(hetero_results, lambda_subset = c(0, 0.5, 1.0, 1.5, 2.0)) {
  
  # Overall histogram
  p1 <- ggplot(hetero_results, aes(x = sr_inf)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Distribution of Inferential Score (All Simulations)",
      x = expression(S[Inf]),
      y = "Count"
    ) +
    theme_minimal()
  
  # Density by lambda
  subset_data <- hetero_results[hetero_strength %in% lambda_subset]
  p2 <- ggplot(subset_data, aes(x = sr_inf, fill = as.factor(hetero_strength))) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = "Inferential Score by Heteroskedasticity Level",
      x = expression(S[Inf]),
      y = "Density",
      fill = "Lambda"
    ) +
    theme_minimal()
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

#' Plot finite-sample bias under null
#'
#' Visualizes the relationship between 1/sqrt(N) and mean sr_ratio under null.
#'
#' @param null_summary data.table with N, inv_sqrt_N, mean_sr_ratio, se_sr_ratio columns
#' @param intercept_est Estimated intercept from regression
#' @param slope_est Estimated slope (c_empirical) from regression
#'
#' @return ggplot object
plot_finite_sample_bias <- function(null_summary, intercept_est, slope_est) {
  
  p <- ggplot(null_summary, aes(x = inv_sqrt_N, y = mean_sr_ratio)) +
    geom_point(size = 3, color = "steelblue") +
    geom_errorbar(aes(ymin = mean_sr_ratio - 1.96*se_sr_ratio, 
                      ymax = mean_sr_ratio + 1.96*se_sr_ratio), 
                  width = 0.002, color = "steelblue") +
    geom_abline(intercept = intercept_est, slope = slope_est, 
                color = "red", linetype = "solid", size = 1) +
    labs(
      title = "Finite-Sample Bias in Raw Ratio Under Null",
      subtitle = sprintf("Empirical fit: c = %.3f", slope_est),
      x = expression(1/sqrt(N)),
      y = "Mean Raw Ratio"
    ) +
    theme_minimal() +
    annotate("text", x = 0.10, y = intercept_est + slope_est * 0.05, 
             label = sprintf("c = %.3f", slope_est), color = "red", size = 5)
  
  p
}

#' Plot sensitivity analysis: simple vs multiple regression
#'
#' @param df_summary Aggregated data.table with N, model, inv_sqrt_N, mean_sr_ratio, se_sr_ratio
#' @param c_simple Scaling factor for simple regression
#' @param c_multiple Scaling factor for multiple regression
#' @param intercept_simple Intercept for simple regression fit
#' @param intercept_multiple Intercept for multiple regression fit
#'
#' @return ggplot object
plot_sensitivity_comparison <- function(df_summary, c_simple, c_multiple, 
                                         intercept_simple, intercept_multiple) {
  
  p <- ggplot(df_summary, aes(x = inv_sqrt_N, y = mean_sr_ratio, color = model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_sr_ratio - 1.96*se_sr_ratio, 
                      ymax = mean_sr_ratio + 1.96*se_sr_ratio), 
                  width = 0.003) +
    geom_abline(intercept = intercept_simple, slope = c_simple, 
                color = "steelblue", linetype = "solid", size = 1) +
    geom_abline(intercept = intercept_multiple, slope = c_multiple, 
                color = "coral", linetype = "solid", size = 1) +
    scale_color_manual(values = c("Simple (p=2)" = "steelblue", "Multiple (p=3)" = "coral")) +
    labs(
      title = "Finite-Sample Bias Decreases with Model Complexity",
      subtitle = sprintf("Simple (p=2): c = %.3f | Multiple (p=3): c = %.3f", c_simple, c_multiple),
      x = expression(1/sqrt(N)),
      y = "Mean Raw Ratio",
      color = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  p
}

#' Plot T_SInf vs Raw Ratio comparison
#'
#' Side-by-side comparison showing N-dependence of T_SInf vs near-invariance of Raw Ratio.
#'
#' @param comparison_summary data.table with N, hetero_strength, mean_T_sinf, mean_sr_ratio
#'
#' @return arranged ggplot grid
plot_tsinf_vs_ratio <- function(comparison_summary) {
  
  # T_SInf (N-dependent)
  p1 <- ggplot(comparison_summary, aes(x = N, y = mean_T_sinf, color = as.factor(hetero_strength))) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    scale_color_viridis_d(option = "plasma", end = 0.9) +
    labs(
      title = expression("Scaled Inferential Score " * T[SInf] * " vs N"),
      subtitle = "Strong N-dependence: lines diverge with increasing N",
      y = expression("Mean " * T[SInf]),
      x = "N (log scale)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Raw Ratio (near-invariant)
  p2 <- ggplot(comparison_summary, aes(x = N, y = mean_sr_ratio, color = as.factor(hetero_strength))) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    scale_color_viridis_d(option = "plasma", end = 0.9) +
    labs(
      title = "Raw Ratio vs N",
      subtitle = "Much weaker N-dependence: lines remain approximately parallel",
      y = "Mean Raw Ratio",
      x = "N (log scale)",
      color = expression(lambda)
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(0.45, 0.55))
}

#' Plot scaling factor comparison (Raw vs Adjusted)
#'
#' @param scaling_long Melted data.table with N, hetero_strength, Metric, Value
#'
#' @return ggplot object with facets
plot_scaling_comparison <- function(scaling_long) {
  
  p <- ggplot(scaling_long, aes(x = N, y = Value, color = as.factor(hetero_strength))) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    scale_color_viridis_d(option = "plasma", end = 0.9) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(
      title = "Comparison: Raw Ratio vs Empirically-Adjusted Reliability Score",
      subtitle = "Flatter lines indicate better N-invariance",
      x = "Sample Size (log scale)",
      y = "Score Value",
      color = expression(lambda)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p
}

#' Plot coverage degradation under heteroskedasticity
#'
#' @param coverage_summary data.table with N, hetero_strength, mean_coverage columns
#'
#' @return ggplot object
plot_coverage_degradation <- function(coverage_summary) {
  
  p <- ggplot(coverage_summary, aes(x = hetero_strength, y = mean_coverage, color = as.factor(N))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
    labs(
      title = "Classical CI Coverage Degradation Under Heteroskedasticity",
      subtitle = "Coverage of classical Wald intervals (using SE_classic)",
      x = "Heteroskedasticity Strength (lambda)",
      y = "Coverage (%)",
      color = "Sample Size (N)"
    ) +
    theme_minimal() +
    coord_cartesian(ylim = c(50, 100))
  
  p
}
