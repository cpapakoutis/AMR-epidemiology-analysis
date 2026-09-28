############################################################
# Project: WHO GLASS Antimicrobial Resistance Analysis
# Script: 00_setup_packages.R
#
# Purpose:
#   Ensure that the R packages directly required by the
#   analysis scripts are installed and loaded.
#
# Notes:
#   - Only direct project dependencies are listed here.
#   - Optional / unused packages are intentionally excluded.
#   - Explicit library() calls are used so dependency tools
#     such as renv can detect the project requirements.
############################################################


# ----------------------------------------------------------
# 1. Define direct project dependencies
# ----------------------------------------------------------

required_packages <- c(
  "tidyverse",   # data manipulation, import, and ggplot2
  "here",        # project-relative file paths
  "viridis",     # colour scales used in exploratory figures
  "sandwich",    # cluster-robust covariance estimators
  "lmtest",      # robust coefficient tests
  "lme4",        # binomial mixed-effects models
  "glmmTMB"      # beta-binomial mixed-effects models
)


# ----------------------------------------------------------
# 2. Install only packages that are missing
# ----------------------------------------------------------

installed_packages <- rownames(installed.packages())

missing_packages <- setdiff(
  required_packages,
  installed_packages
)

if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}


# ----------------------------------------------------------
# 3. Load direct project dependencies
# ----------------------------------------------------------
#
# These explicit library() calls are intentional:
# renv can detect them when building the project lockfile.

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(here)
  library(viridis)
  library(sandwich)
  library(lmtest)
  library(lme4)
  library(glmmTMB)
  
})


# ----------------------------------------------------------
# 4. Confirmation message
# ----------------------------------------------------------

message(
  "All required project packages are installed and loaded."
)

