# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 04_rq1_specifications.R
#
# RQ1: contribution of cognitive skills and Big Five personality
# traits to wage prediction (Section 5.2).
#   S0 = controls (including country dummies)
#   S1 = S0 + cognitive skills
#   S2 = S0 + Big Five
#   S3 = S0 + cognitive skills + Big Five
# The contribution of each block is the change in the test-set R²,
# RMSE and MAE when the block is added.
# The "_nojob" specifications exclude the job-characteristics
# controls (required education, sector, occupational skill level,
# firm size, weekly hours).
# Every change is accompanied by a 95% paired-bootstrap interval on
# the test set (Section C2). Cognitive skills use plausible value 1;
# standard errors are HC3.

# Run from the repository root (see README).
# ---- PACKAGES ----

library(tidyverse)
library(sandwich)   # vcovHC()
library(lmtest)     # coeftest()


# ---- LOAD DATA AND VARIABLE BLOCKS ----
# Blocks (personal_controls, job_controls, controls, cognitive,
# bigfive) and helpers (make_formula, check_blocks) are defined once
# in 00_variable_blocks.R.

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
test  <- readRDS("data/processed/test.rds")

# Stops the script if train contains a column not assigned to any
# block, so that no specification can change silently.
check_blocks(train)

cat("Training set:", nrow(train), "| Test set:", nrow(test), "\n\n")


# ---- A. SPECIFICATIONS ----

specs <- list(
  S0       = controls,
  S1       = c(controls, cognitive),
  S2       = c(controls, bigfive),
  S3       = c(controls, cognitive, bigfive),
  S0_nojob = personal_controls,
  S1_nojob = c(personal_controls, cognitive),
  S2_nojob = c(personal_controls, bigfive),
  S3_nojob = c(personal_controls, cognitive, bigfive)
)

models <- map(specs, ~ lm(make_formula(.x), data = train))

# Check: S3 must equal the model with every non-auxiliary column,
# estimated independently of the blocks (the 03 baseline).
baseline_check <- lm(LOG_WAGE ~ ., data = train %>% select(-all_of(aux_vars)))
stopifnot(isTRUE(all.equal(unname(predict(baseline_check, newdata = test)),
                           unname(predict(models$S3, newdata = test)))))
cat("Check passed: S3 reproduces the 03 baseline.\n\n")


# ---- B. OUT-OF-SAMPLE PERFORMANCE ----

oos_metrics <- function(model, newdata) {
  pred <- predict(model, newdata = newdata)
  y    <- newdata$LOG_WAGE
  list(
    RMSE = sqrt(mean((y - pred)^2)),
    MAE  = mean(abs(y - pred)),
    R2   = 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  )
}

perf <- imap_dfr(models, function(m, s) {
  bind_cols(
    tibble(spec = s,
           n_coef = length(coef(m)),
           adj_R2_train = summary(m)$adj.r.squared),
    as_tibble(oos_metrics(m, test))
  )
})

cat("=== RQ1 — OUT-OF-SAMPLE PERFORMANCE BY SPECIFICATION (test set) ===\n")
perf %>%
  mutate(across(c(adj_R2_train, RMSE, MAE, R2), ~ round(.x, 4))) %>%
  print(n = Inf)

write_csv(perf, "output/tables/table_3x_rq1_performance.csv")


# ---- C. CONTRIBUTION OF EACH BLOCK ----
# Change in each metric when a block is added. A positive dR2 and a
# negative dRMSE / dMAE mean the block improves prediction.

get_metric <- function(s, metric) perf[[metric]][perf$spec == s]

delta_row <- function(label, from, to) {
  tibble(
    comparison = label,
    from = from,
    to   = to,
    dR2   = get_metric(to, "R2")   - get_metric(from, "R2"),
    dRMSE = get_metric(to, "RMSE") - get_metric(from, "RMSE"),
    dMAE  = get_metric(to, "MAE")  - get_metric(from, "MAE")
  )
}

deltas <- bind_rows(
  # Main specifications (all controls)
  delta_row("Cognitive skills alone",           "S0", "S1"),
  delta_row("Big Five alone",                   "S0", "S2"),
  delta_row("Big Five beyond cognitive skills", "S1", "S3"),
  delta_row("Cognitive skills beyond Big Five", "S2", "S3"),
  delta_row("Both blocks together",             "S0", "S3"),
  # Comparison: without job-characteristics controls
  delta_row("Cognitive skills alone",           "S0_nojob", "S1_nojob"),
  delta_row("Big Five alone",                   "S0_nojob", "S2_nojob"),
  delta_row("Big Five beyond cognitive skills", "S1_nojob", "S3_nojob"),
  delta_row("Cognitive skills beyond Big Five", "S2_nojob", "S3_nojob"),
  delta_row("Both blocks together",             "S0_nojob", "S3_nojob")
)

cat("\n=== RQ1 — CONTRIBUTION OF EACH BLOCK (test set) ===\n")
deltas %>%
  mutate(across(c(dR2, dRMSE, dMAE), ~ round(.x, 4))) %>%
  print(n = Inf)

write_csv(deltas, "output/tables/table_3x_rq1_deltas.csv")


# ---- C2. BOOTSTRAP CONFIDENCE INTERVALS FOR THE CONTRIBUTIONS ----
# Paired bootstrap on the test set:
#   - the models are not re-estimated: their test-set predictions
#     are computed once and kept fixed;
#   - in each replication, test observations are resampled with
#     replacement within each country;
#   - the same resampled observations are used for every model, so
#     each difference compares the models on the same workers;
#   - 95% percentile intervals from B = 2,000 replications.
# The intervals reflect the variability of the test set only.

pred_mat <- sapply(models, function(m) predict(m, newdata = test))  # n_test x 8
y_test   <- test$LOG_WAGE

from_idx <- match(deltas$from, colnames(pred_mat))
to_idx   <- match(deltas$to,   colnames(pred_mat))

# Test-set row indices by country (for stratified resampling)
idx_by_country <- split(seq_len(nrow(test)), test$CNTRYID)

B <- 2000
set.seed(123)

boot_array <- array(NA_real_,
                    dim = c(B, nrow(deltas), 3),
                    dimnames = list(NULL, NULL, c("dR2", "dRMSE", "dMAE")))

for (b in seq_len(B)) {
  idx <- unlist(lapply(idx_by_country,
                       function(i) i[sample.int(length(i), replace = TRUE)]),
                use.names = FALSE)
  yb  <- y_test[idx]
  err <- yb - pred_mat[idx, , drop = FALSE]   # n_test x 8 errors
  sse  <- colSums(err^2)
  rmse <- sqrt(sse / length(yb))
  mae  <- colMeans(abs(err))
  r2   <- 1 - sse / sum((yb - mean(yb))^2)
  boot_array[b, , "dR2"]   <- r2[to_idx]   - r2[from_idx]
  boot_array[b, , "dRMSE"] <- rmse[to_idx] - rmse[from_idx]
  boot_array[b, , "dMAE"]  <- mae[to_idx]  - mae[from_idx]
}

# Note: the argument is called metric_name, not metric. Inside
# tibble(), a column created earlier in the same call ("metric") would
# otherwise mask the function argument in the following lines.
boot_ci <- map_dfr(c("dR2", "dRMSE", "dMAE"), function(metric_name) {
  tibble(
    comparison = deltas$comparison,
    from       = deltas$from,
    to         = deltas$to,
    metric     = metric_name,
    estimate   = deltas[[metric_name]],
    ci_low     = apply(boot_array[, , metric_name], 2, quantile, probs = 0.025),
    ci_high    = apply(boot_array[, , metric_name], 2, quantile, probs = 0.975)
  )
}) %>%
  mutate(
    # Improvement = higher R², lower RMSE / MAE
    improves = if_else(metric == "dR2", ci_low > 0, ci_high < 0)
  )

cat("\n=== RQ1 — CONTRIBUTION OF EACH BLOCK WITH 95% BOOTSTRAP CI ===\n")
cat("(paired bootstrap on the test set, B =", B, ", stratified by country)\n\n")
boot_ci %>%
  filter(metric %in% c("dR2", "dRMSE")) %>%
  mutate(value = sprintf("%+.4f [%+.4f, %+.4f]%s",
                         estimate, ci_low, ci_high, if_else(improves, " *", ""))) %>%
  select(comparison, from, to, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  print(n = Inf, width = Inf)
cat("\n* = 95% interval excludes zero in the direction of improvement\n")

write_csv(boot_ci, "output/tables/table_3x_rq1_deltas_bootstrap.csv")


# ---- D. COEFFICIENTS PER STANDARD DEVIATION ----
# Cognitive skills are on the PIAAC 0-500 point scale: coefficients
# are multiplied by the training-sample SD, so they read as the
# change in log wage per 1 SD of skill. This is a pure rescaling:
# predictions and p-values do not change.
# Big Five scores are already expressed in within-country SD units
# (standardised by the OECD in the Public Use File), so they are
# reported as they are.
# pct_wage = 100 * (exp(coef) - 1): approximate % change in the
# hourly wage per 1 SD.

scale_factor <- c(
  setNames(map_dbl(cognitive, ~ sd(train[[.x]])), cognitive),
  setNames(rep(1, length(bigfive)), bigfive)
)

cat("\nSD of cognitive skills in the training sample (points):\n")
print(round(scale_factor[cognitive], 1))

coef_specs <- c("S1", "S2", "S3", "S1_nojob", "S2_nojob", "S3_nojob")

coef_table <- imap_dfr(models[coef_specs], function(m, s) {
  ct    <- coeftest(m, vcov. = vcovHC(m, type = "HC3"))
  terms <- intersect(rownames(ct), c(cognitive, bigfive))
  tibble(
    spec     = s,
    term     = terms,
    coef_sd  = ct[terms, "Estimate"]   * scale_factor[terms],
    se_sd    = ct[terms, "Std. Error"] * scale_factor[terms],
    p_value  = ct[terms, "Pr(>|t|)"],
    pct_wage = 100 * (exp(coef_sd) - 1)
  )
})

stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*",
            p < 0.1 ~ ".", TRUE ~ "")
}

cat("\n=== RQ1 — % WAGE CHANGE PER 1 SD (HC3 p-values, PV1 only) ===\n")
cat("Signif.: *** 0.001  ** 0.01  * 0.05  . 0.1\n\n")
coef_table %>%
  mutate(cell = sprintf("%+.1f%%%s", pct_wage, stars(p_value))) %>%
  select(term, spec, cell) %>%
  pivot_wider(names_from = spec, values_from = cell, values_fill = "") %>%
  print(n = Inf)

write_csv(coef_table, "output/tables/table_3x_rq1_coefficients.csv")


# END OF SCRIPT 04 ----
