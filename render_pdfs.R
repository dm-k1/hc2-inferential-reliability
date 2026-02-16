#!/usr/bin/env Rscript
library(rmarkdown)

cat("=== HC2 Reliability Score: PDF Rendering Script ===\n")
cat("This script renders the R Markdown notebooks to PDF.\n")
cat("Rendering can take 30-60 minutes depending on your system.\n\n")

files <- list(
  list(
    file = "00_HC_estimators_validation.Rmd",
    desc = "Part 0: OLS and HC Estimator Validation",
    time = "~10 minutes"
  ),
  list(
    file = "01_null_calibration.Rmd",
    desc = "Part 1: Null Calibration of Inferential Score",
    time = "~20 minutes"
  ),
  list(
    file = "02_HC2_validation.Rmd",
    desc = "Part 2: HC2 Validation Under Heteroskedasticity",
    time = "~30 minutes"
  )
)

for (item in files) {
  f <- item$file
  desc <- item$desc
  
  cat(sprintf(">>> Rendering: %s\n", desc))
  cat(sprintf("    File: %s | Est. time: %s\n", f, item$time))
  
  tryCatch(
    {
      render(
        f,
        output_format = "pdf_document",
        output_options = list(latex_engine = "pdflatex"),
        quiet = TRUE
      )
      pdf_file <- sub("\\.Rmd$", ".pdf", f)
      cat(sprintf("    ✓ Successfully created: %s\n\n", pdf_file))
    },
    error = function(e) {
      cat(sprintf("    ✗ Error: %s\n\n", e$message))
    }
  )
}

cat("PDF rendering complete!\n")
cat("Check the project root directory for .pdf files.\n")
