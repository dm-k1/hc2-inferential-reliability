## ============================================================
## 00_config.R
## 
## Purpose:
##   - Global configuration for all simulations
##   - Paths, seeds, parameter grids, number of simulations
##   - Package loading and dependencies
## ============================================================

# ---- Packages ----
required_packages <- c(
  "data.table",
  "dplyr",
  "lmerTest",
  "sandwich",
  "future",      # For parallel execution plan
  "furrr",       # For future_lapply (parallel map)
  "parallelly",  # For availableCores()
  "knitr",
  "ggplot2",
  "lmtest",      # For benchmarking (BP/White tests)
  "gridExtra",   # For arranging plots
  "reshape2"     # For melting data
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ---- Seed & reproducibility ----
# Global seed for reproducibility; parallel seeds handled via future.seed = TRUE in furrr functions
set.seed(2024)

# ---- File paths ----
BASE_PATH <- getwd()
R_PATH <- file.path(BASE_PATH, "R")
RESULTS_PATH <- file.path(BASE_PATH, "results")

# Create results directory if it doesn't exist
if (!dir.exists(RESULTS_PATH)) {
  dir.create(RESULTS_PATH, recursive = TRUE)
}

# ---- Part 1: Null Calibration Grids ----
# N grid for null calibration (geometric progression)
N_GRID_NULL <- c(50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400)
N_SIMS_PER_N_NULL <- 10000  # Reduced to 10,000 for performance (was 100,000)

# HC types to compare
HC_TYPES <- c("HC1", "HC2", "HC3")

# Scaling exponents to test
ALPHA_GRID <- seq(0.3, 0.7, by = 0.05)

# ---- Part 2: Heteroskedastic Validation Grids ----
# N grid for hetero validation
N_GRID_HETERO <- c(50, 100, 200, 400, 800, 1600)

# Heteroskedasticity strength (lambda)
HETERO_STRENGTH_GRID <- seq(0, 2, by = 0.2)

# Effect size (beta_x)
BETA_X_GRID <- c(0.2, 0.5, 0.8)

# Baseline noise (sigma0)
SIGMA0_GRID <- c(0.5, 1.0, 2.0)

# Simulations per cell (N x lambda x beta x sigma0)
N_SIMS_PER_CELL_HETERO <- 500

# ---- Regression specification ----
# Y = beta_0 + beta_x * X + epsilon
COEF_OF_INTEREST <- "x"  # Standardized to lowercase 'x'

# ---- Parallel computing ----
# Set up future backend for parallelization
# Use multicore on Unix-like systems and multisession on Windows for portability
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = max(1, parallelly::availableCores() - 1))
} else {
  plan(multicore, workers = max(1, parallelly::availableCores() - 1))
}

# ---- Confidence level for coverage studies ----
CONF_LEVEL <- 0.95  # 95% confidence intervals

# ---- Output settings ----
SAVE_INTERMEDIATE <- TRUE  # save intermediate results
VERBOSE <- TRUE  # print progress messages
