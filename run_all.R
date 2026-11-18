# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using
# cross-national adult skills data
# Script: run_all.R
#
# Runs the whole analysis in order. Run from the repository root,
# after placing data/raw/PIAAC_Cycle2.rds (see data/README.md).

dir.create("data/processed",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures",   recursive = TRUE, showWarnings = FALSE)

scripts <- c(
  "R/01_data_cleaning.R",       # cleaning, PPP conversion, train/test split
  "R/02_descriptives.R",        # descriptive statistics (Chapter 3)
  "R/03_ols_assumptions.R",     # OLS diagnostics, VIF, HC3 standard errors
  "R/04_rq1_specifications.R",  # RQ1: nested specifications S0-S3
  "R/05_rq1_pv_weights.R",      # RQ1 robustness: 10 plausible values, weights
  "R/06_rq3_regularized.R",     # RQ3: Ridge, LASSO, Elastic Net
  "R/07_multilevel.R",          # RQ2-RQ3: multilevel models
  "R/08_rq2_sensitivity.R"      # RQ2 robustness: leave-one-country-out
)

for (s in scripts) {
  message("Running ", s)
  source(s, local = new.env())   # each script in a clean environment
}
