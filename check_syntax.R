source_files <- c('R/00_config.R', 'R/10_dgp_and_fits.R', 'R/20_metrics.R', 
                  'R/30_null_calibration.R', 'R/40_hetero_sims.R', 
                  'R/50_invariance_and_lookup.R', 'R/60_tables_and_plots.R')

for (f in source_files) {
  tryCatch({
    parse(file = f)
    cat(paste0("✓ ", f, " syntax OK\n"))
  }, error = function(e) {
    cat(paste0("✗ ERROR in ", f, ": ", e$message, "\n"))
  })
}
cat("\nAll R files parsed successfully!\n")
