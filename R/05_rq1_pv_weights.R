# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 05_rq1_pv_weights.R
#
# RQ1: plausible values (Rubin's rules) and sampling weights
# (Sections 4.4 and 5.2).
#
# PART A - PLAUSIBLE VALUES
# S3 (and S3_nojob) are estimated ten times, once per plausible
# value, and combined with Rubin's rules:
#   estimate      = mean of the 10 estimates
#   within var.   = mean of the 10 squared HC3 standard errors
#   between var.  = variance of the 10 estimates
#   total var.    = within + (1 + 1/10) * between
# The same loop reports how the test-set metrics and the
# contribution of the cognitive skills (S1 vs S0) change across the
# ten plausible values.
#
# PART B - SAMPLING WEIGHTS (robustness check)
# S3 is re-estimated with (1) the final sampling weight SPFWT0 and
# (2) SPFWT0 rescaled so that every country has the same total
# weight.
#
# Coefficients are reported per 1 SD (cognitive skills: SD of the
# plausible value used, in the training sample; Big Five: already in
# within-country SD units).


# Run from the repository root (see README).
# ---- PACKAGES ----

library(tidyverse)
library(sandwich)   # vcovHC()
library(lmtest)     # coeftest()


# ---- LOAD DATA AND VARIABLE BLOCKS ----

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
test  <- readRDS("data/processed/test.rds")
check_blocks(train)

y_test <- test$LOG_WAGE
M <- 10   # number of plausible values


# ---- A. HELPERS ----

# Out-of-sample metrics (same definitions as in 03 and 04)
oos <- function(pred, y) {
  tibble(
    RMSE = sqrt(mean((y - pred)^2)),
    MAE  = mean(abs(y - pred)),
    R2   = 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  )
}

# Predictors of each specification, with plausible value m
spec_vars <- function(spec, m) {
  cog <- cognitive_pv(m)
  switch(spec,
         S1       = c(controls, cog),
         S3       = c(controls, cog, bigfive),
         S3_nojob = c(personal_controls, cog, bigfive))
}

# Readable labels (the same skill has a different column name per PV)
term_label <- function(x) {
  case_when(
    str_starts(x, "PVLIT") ~ "Literacy",
    str_starts(x, "PVNUM") ~ "Numeracy",
    str_starts(x, "PVAPS") ~ "Adaptive problem solving",
    TRUE ~ x
  )
}

stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*",
            p < 0.1 ~ ".", TRUE ~ "")
}



# ---- B PLAUSIBLE VALUES ----

# S0 has no cognitive skills, so it does not depend on the PV
s0_model <- lm(make_formula(controls), data = train)
s0_R2    <- oos(predict(s0_model, newdata = test), y_test)$R2

pv_results <- map(seq_len(M), function(m) {
  
  fits <- map(c(S1 = "S1", S3 = "S3", S3_nojob = "S3_nojob"),
              ~ lm(make_formula(spec_vars(.x, m)), data = train))
  
  # Test-set metrics with plausible value m (train and test)
  metrics <- imap_dfr(fits, function(f, s) {
    bind_cols(tibble(pv = m, spec = s),
              oos(predict(f, newdata = test), y_test))
  })
  
  # Coefficients and HC3 variances, per 1 SD
  coefs <- imap_dfr(fits[c("S3", "S3_nojob")], function(f, s) {
    V     <- vcovHC(f, type = "HC3")
    vars  <- c(cognitive_pv(m), bigfive)
    sds   <- c(map_dbl(cognitive_pv(m), ~ sd(train[[.x]])), rep(1, length(bigfive)))
    tibble(
      pv      = m,
      spec    = s,
      term    = term_label(vars),
      q       = unname(coef(f)[vars] * sds),      # estimate per 1 SD
      u       = unname(diag(V)[vars] * sds^2)     # squared SE per 1 SD
    )
  })
  
  cat("Plausible value", m, "done\n")
  list(metrics = metrics, coefs = coefs)
})

pv_metrics <- map_dfr(pv_results, "metrics")
pv_coefs   <- map_dfr(pv_results, "coefs")


# ---- B1. Rubin's rules ----

rubin <- pv_coefs %>%
  group_by(spec, term) %>%
  summarise(
    estimate  = mean(q),
    var_within  = mean(u),
    var_between = var(q),
    .groups = "drop"
  ) %>%
  mutate(
    var_total = var_within + (1 + 1 / M) * var_between,
    se        = sqrt(var_total),
    # Rubin's degrees of freedom (Inf if the PVs give identical estimates)
    df        = if_else(var_between > 0,
                        (M - 1) * (1 + var_within / ((1 + 1 / M) * var_between))^2,
                        Inf),
    p_value   = 2 * pt(-abs(estimate / se), df),
    # Share of the total variance due to the uncertainty in the skill
    # measurement (differences between PVs)
    share_pv  = (1 + 1 / M) * var_between / var_total,
    pct_wage  = 100 * (exp(estimate) - 1)
  )

# PV1-only results, for comparison with 04
pv1_only <- pv_coefs %>%
  filter(pv == 1) %>%
  transmute(spec, term, estimate_pv1 = q, se_pv1 = sqrt(u))

rubin_table <- rubin %>%
  left_join(pv1_only, by = c("spec", "term")) %>%
  mutate(se_ratio = se / se_pv1)   # > 1: PV1 alone understated the SE

term_order <- c("Literacy", "Numeracy", "Adaptive problem solving", bigfive)

cat("\n=== A1. RUBIN'S RULES — % WAGE CHANGE PER 1 SD (10 PVs, HC3) ===\n")
cat("Signif.: *** 0.001  ** 0.01  * 0.05  . 0.1\n")
cat("se_ratio = Rubin SE / PV1-only SE (> 1: PV1 alone understated uncertainty)\n")
cat("share_pv = share of total variance due to skill-measurement uncertainty\n\n")
rubin_table %>%
  mutate(term = factor(term, levels = term_order)) %>%
  arrange(spec, term) %>%
  transmute(
    spec, term,
    pct_PV1   = sprintf("%+.1f%%", 100 * (exp(estimate_pv1) - 1)),
    pct_rubin = sprintf("%+.1f%%%s", pct_wage, stars(p_value)),
    se_pv1    = round(se_pv1, 4),
    se_rubin  = round(se, 4),
    se_ratio  = round(se_ratio, 2),
    share_pv  = round(share_pv, 3)
  ) %>%
  print(n = Inf, width = Inf)

write_csv(rubin_table, "output/tables/table_3x_rq1_rubin.csv")


# ---- B2. Test-set metrics across the 10 PVs ----

cat("\n=== A2. TEST-SET METRICS ACROSS THE 10 PLAUSIBLE VALUES ===\n")
pv_metrics %>%
  group_by(spec) %>%
  summarise(across(c(RMSE, MAE, R2),
                   list(PV1 = ~ .x[pv == 1], mean = mean, min = min, max = max))) %>%
  pivot_longer(-spec, names_to = c("metric", "stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(across(c(PV1, mean, min, max), ~ round(.x, 4))) %>%
  print(n = Inf)

# Contribution of cognitive skills (S1 vs S0) for each PV
dR2_skills <- pv_metrics %>%
  filter(spec == "S1") %>%
  transmute(pv, dR2_skills = R2 - s0_R2)

cat("\nContribution of cognitive skills (dR2, S1 vs S0) across the 10 PVs:\n")
cat("PV1:", round(dR2_skills$dR2_skills[dR2_skills$pv == 1], 4),
    "| mean:", round(mean(dR2_skills$dR2_skills), 4),
    "| range: [", round(min(dR2_skills$dR2_skills), 4), ",",
    round(max(dR2_skills$dR2_skills), 4), "]\n")
cat("(Compare the range with the 95% bootstrap interval for PV1 in 04:\n",
    "if the range lies well inside it, the choice of PV1 does not matter.)\n")

write_csv(pv_metrics, "output/tables/table_3x_rq1_pv_metrics.csv")



# ---- C SAMPLING WEIGHTS (robustness check) ----

stopifnot(all(train[[weight_var]] > 0))

# Both weights rescaled to mean 1 (only relative weights matter)
train_w <- train %>%
  group_by(CNTRYID) %>%
  mutate(w_country_equal = .data[[weight_var]] / sum(.data[[weight_var]])) %>%
  ungroup() %>%
  mutate(
    w_raw           = .data[[weight_var]] / mean(.data[[weight_var]]),
    w_country_equal = w_country_equal / mean(w_country_equal)
  )

cat("\n=== B0. SHARE OF EACH COUNTRY IN THE ESTIMATION ===\n")
cat("(unweighted = share of observations; raw = share of SPFWT0;\n",
    " equal = identical share for every country)\n")
train_w %>%
  group_by(CNTRYID) %>%
  summarise(unweighted = n(), raw = sum(w_raw), equal = sum(w_country_equal)) %>%
  mutate(across(c(unweighted, raw, equal), ~ round(100 * .x / sum(.x), 1))) %>%
  arrange(desc(raw)) %>%
  print(n = Inf)

f_s3 <- make_formula(c(controls, cognitive, bigfive))

fit_unw <- lm(f_s3, data = train_w)
fit_raw <- lm(f_s3, data = train_w, weights = w_raw)
fit_eq  <- lm(f_s3, data = train_w, weights = w_country_equal)

sd_pv1 <- c(map_dbl(cognitive, ~ sd(train[[.x]])), rep(1, length(bigfive)))
names(sd_pv1) <- c(cognitive, bigfive)

weighted_table <- imap_dfr(
  list(`Unweighted`            = fit_unw,
       `Raw SPFWT0`            = fit_raw,
       `Equal country weights` = fit_eq),
  function(f, s) {
    ct    <- coeftest(f, vcov. = vcovHC(f, type = "HC3"))
    vars  <- c(cognitive, bigfive)
    met   <- oos(predict(f, newdata = test), y_test)   # unweighted test metrics
    tibble(
      model    = s,
      term     = term_label(vars),
      coef_sd  = unname(ct[vars, "Estimate"]   * sd_pv1[vars]),
      se_sd    = unname(ct[vars, "Std. Error"] * sd_pv1[vars]),
      p_value  = unname(ct[vars, "Pr(>|t|)"]),
      pct_wage = 100 * (exp(coef_sd) - 1),
      test_RMSE = met$RMSE,
      test_R2   = met$R2
    )
  }
)

cat("\n=== B1. S3 UNWEIGHTED vs WEIGHTED — % WAGE CHANGE PER 1 SD (PV1, HC3) ===\n")
cat("Signif.: *** 0.001  ** 0.01  * 0.05  . 0.1\n\n")
weighted_table %>%
  mutate(cell = sprintf("%+.1f%%%s", pct_wage, stars(p_value)),
         term = factor(term, levels = term_order)) %>%
  select(term, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell) %>%
  arrange(term) %>%
  print(n = Inf)

cat("\nTest-set metrics (unweighted, same test set as 03/04):\n")
weighted_table %>%
  distinct(model, test_RMSE, test_R2) %>%
  mutate(across(c(test_RMSE, test_R2), ~ round(.x, 4))) %>%
  print()

write_csv(weighted_table, "output/tables/table_3x_rq1_weighted.csv")


# END OF SCRIPT 05 ----
