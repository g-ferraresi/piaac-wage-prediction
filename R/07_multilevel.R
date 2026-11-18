# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 07_multilevel.R
# RQ2 and RQ3 — multilevel models (Sections 4.6, 5.3 and 5.4)
#
# Reads train.rds / test.rds (01_data_cleaning.R) and the variable
# blocks in 00_variable_blocks.R.
#
# Models:
#   ML0  random intercept for country            (analogue of OLS S3)
#   ML1  random intercept + uncorrelated random slopes on the three
#        cognitive skills, (1 + LIT + NUM + APS || country)
#                                                (analogue of OLS S4)
#   ML1 + one Big Five trait: five variants, each adding a random
#        slope on ONE trait to ML1
#
# Rules fixed in advance:
#   - The country dummies of OLS are replaced by the random intercept
#     (a model cannot contain both for the same grouping).
#   - Only variances of the random slopes are estimated, not their
#     correlations ("||"): 12 countries cannot support more.
#   - If ML1 is singular (a slope variance estimated at zero), the
#     slope with the smallest variance is removed and the model is
#     refitted, until the fit is not singular.
#   - Estimation by REML (lme4 default). Likelihood-ratio tests with
#     anova(), which refits the models by ML (Bates et al., 2015).
#     A test of a variance = 0 lies on the boundary of the parameter
#     space: the chi-square p-value is conservative, and for a single
#     variance the correct p-value is half of the chi-square(1)
#     p-value (Section 4.6.4).
#   - Test-set predictions are CONDITIONAL on the estimated country
#     effects (predict() default, re.form = NULL): every test country
#     is also in the training set.
#
# Cognitive skills (PV1) are standardized with the training mean and
# SD, exactly as in 06_rq3_regularized.R. With "||" this matters: the
# zero of the skill scale is the training mean (Bates et al., 2015,
# section on uncorrelated random effects). Big Five are within-country
# z-scores (OECD). R² is the out-of-sample R² used in 03-06.


# Run from the repository root (see README).
# ---- PACKAGES ----

library(tidyverse)
library(lme4)
library(sandwich)   # vcovHC()
library(lmtest)     # waldtest()


# ---- LOAD DATA AND VARIABLE BLOCKS ----

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
test  <- readRDS("data/processed/test.rds")
check_blocks(train)

country_labels <- c(
  "152" = "Chile", "191" = "Croatia", "203" = "Czechia",
  "250" = "France", "372" = "Ireland", "376" = "Israel",
  "380" = "Italy", "428" = "Latvia", "440" = "Lithuania",
  "554" = "New Zealand", "616" = "Poland", "724" = "Spain"
)


# ---- STANDARDIZE COGNITIVE SKILLS (as in 06) ----

skill_mean <- setNames(map_dbl(cognitive, ~ mean(train[[.x]])), cognitive)
skill_sd   <- setNames(map_dbl(cognitive, ~ sd(train[[.x]])),   cognitive)

standardize_skills <- function(data) {
  for (v in cognitive) data[[v]] <- (data[[v]] - skill_mean[v]) / skill_sd[v]
  data
}

train_s <- standardize_skills(train)
test_s  <- standardize_skills(test)
y_test  <- test_s$LOG_WAGE


# ---- HELPERS (oos and paired_bootstrap as in 06) ----

oos <- function(pred, y) {
  tibble(
    RMSE = sqrt(mean((y - pred)^2)),
    MAE  = mean(abs(y - pred)),
    R2   = 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  )
}

paired_bootstrap <- function(pred_mat, y, strata, comparisons, B = 2000, seed = 123) {
  from_idx <- match(comparisons$from, colnames(pred_mat))
  to_idx   <- match(comparisons$to,   colnames(pred_mat))
  stopifnot(!anyNA(from_idx), !anyNA(to_idx))
  idx_by_stratum <- split(seq_along(y), strata)
  out <- array(NA_real_, dim = c(B, nrow(comparisons), 3),
               dimnames = list(NULL, NULL, c("dR2", "dRMSE", "dMAE")))
  set.seed(seed)
  for (b in seq_len(B)) {
    idx <- unlist(lapply(idx_by_stratum,
                         function(i) i[sample.int(length(i), replace = TRUE)]),
                  use.names = FALSE)
    yb  <- y[idx]
    err <- yb - pred_mat[idx, , drop = FALSE]
    sse  <- colSums(err^2)
    rmse <- sqrt(sse / length(yb))
    mae  <- colMeans(abs(err))
    r2   <- 1 - sse / sum((yb - mean(yb))^2)
    out[b, , "dR2"]   <- r2[to_idx]   - r2[from_idx]
    out[b, , "dRMSE"] <- rmse[to_idx] - rmse[from_idx]
    out[b, , "dMAE"]  <- mae[to_idx]  - mae[from_idx]
  }
  full <- apply(pred_mat, 2, function(p) unlist(oos(p, y)))
  map_dfr(c("dR2", "dRMSE", "dMAE"), function(metric_name) {
    base_name <- sub("^d", "", metric_name)
    tibble(
      comparison = comparisons$label,
      from       = comparisons$from,
      to         = comparisons$to,
      metric     = metric_name,
      estimate   = full[base_name, to_idx] - full[base_name, from_idx],
      ci_low     = apply(out[, , metric_name, drop = FALSE], 2, quantile, probs = 0.025),
      ci_high    = apply(out[, , metric_name, drop = FALSE], 2, quantile, probs = 0.975)
    )
  }) %>%
    mutate(improves = if_else(metric == "dR2", ci_low > 0, ci_high < 0),
           worsens  = if_else(metric == "dR2", ci_high < 0, ci_low > 0))
}

# Multilevel formula: fixed part = S3 without the country dummies,
# random part = country intercept (+ uncorrelated slopes if given)
fixed_vars <- c(setdiff(controls, "CNTRYID"), cognitive, bigfive)

ml_formula <- function(slopes = character(0)) {
  random_part <- if (length(slopes) == 0) {
    "(1 | CNTRYID)"
  } else {
    paste0("(1 + ", paste(slopes, collapse = " + "), " || CNTRYID)")
  }
  as.formula(paste("LOG_WAGE ~", paste(fixed_vars, collapse = " + "), "+", random_part))
}

# Fit with REML, recording convergence warnings and singularity
fit_ml <- function(slopes = character(0)) {
  warnings_seen <- character(0)
  m <- withCallingHandlers(
    lmer(ml_formula(slopes), data = train_s, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(model = m, slopes = slopes, singular = isSingular(m), warnings = warnings_seen)
}

# Variance components of a fitted model, as a tidy table
var_components <- function(m, label) {
  as.data.frame(VarCorr(m)) %>%
    as_tibble() %>%   # a tibble, so that print(n = Inf) works
    transmute(model = label,
              component = if_else(grp == "Residual", "Residual", var1),
              variance = vcov, sd = sdcor)
}


# ---- OLS BENCHMARKS (same as 06) ----

f_ols_s3 <- make_formula(c(controls, cognitive, bigfive))
f_ols_s4 <- make_formula(c(controls, cognitive, bigfive, paste0("CNTRYID:", cognitive)))

ols_s3 <- lm(f_ols_s3, data = train_s)
ols_s4 <- lm(f_ols_s4, data = train_s)


# ---- A. ML0 — RANDOM INTERCEPT ----

ml0 <- fit_ml()
cat("=== ML0: random intercept ===\n")
cat("Singular:", ml0$singular,
    "| Warnings:", if (length(ml0$warnings)) paste(ml0$warnings, collapse = "; ") else "none", "\n\n")


# ---- B. ML1 — RANDOM SLOPES ON COGNITIVE SKILLS (with the fallback rule) ----

slopes <- cognitive
repeat {
  ml1 <- fit_ml(slopes)
  if (!ml1$singular || length(slopes) == 0) break
  vc <- as.data.frame(VarCorr(ml1$model))
  slope_var <- vc %>% filter(var1 %in% slopes)
  drop_slope <- slope_var$var1[which.min(slope_var$vcov)]
  cat("Singular fit with slopes on", paste(slopes, collapse = ", "),
      "-> removing the random slope with the smallest variance:", drop_slope, "\n")
  slopes <- setdiff(slopes, drop_slope)
}

cat("\n=== ML1: random slopes retained on:",
    if (length(slopes)) paste(slopes, collapse = ", ") else "none", "===\n")
cat("Singular:", ml1$singular,
    "| Warnings:", if (length(ml1$warnings)) paste(ml1$warnings, collapse = "; ") else "none", "\n")


# ---- B1. Variance components ----

vc_table <- bind_rows(var_components(ml0$model, "ML0"),
                      var_components(ml1$model, "ML1"))

cat("\n=== VARIANCE COMPONENTS (REML) ===\n")
vc_table %>% mutate(across(c(variance, sd), ~ signif(.x, 4))) %>% print(n = Inf)

write_csv(vc_table, "output/tables/table_3x_ml_variance_components.csv")


# ---- B2. Likelihood-ratio test: are the random slopes needed? ----

cat("\n=== LR TEST: ML0 vs ML1 (refitted by ML) ===\n")
lrt_slopes <- anova(ml0$model, ml1$model)
print(lrt_slopes)
cat("(Boundary test: the chi-square p-value is conservative.)\n")


# ---- B3. Joint test on the 33 country x skill interactions of S4 ----
# Fixed-effects analogue of the LR test above: H0 = all interactions
# are zero (same returns in all countries). HC3-robust Wald test.

cat("\n=== JOINT TEST ON THE COUNTRY x SKILL INTERACTIONS (OLS S4) ===\n")
joint_s4 <- waldtest(ols_s3, ols_s4, vcov = function(m) vcovHC(m, type = "HC3"))
print(joint_s4)


# ---- C. RQ2 — COUNTRY-SPECIFIC RETURNS TO COGNITIVE SKILLS ----
# Multilevel: country return = common slope + country deviation
# (conditional modes, i.e. shrunk towards the common slope).
# OLS S4:     country return = main effect + country interaction
# (no shrinkage). Comparing the two shows how much the multilevel
# model pulls each country towards the common return. The multilevel
# returns are DESCRIPTIVE: with 12 countries they are not precise
# estimates of single-country returns.

country_codes <- rownames(coef(ml1$model)$CNTRYID)

ml_slopes <- coef(ml1$model)$CNTRYID[, cognitive, drop = FALSE] %>%
  as_tibble(rownames = "CNTRYID") %>%
  pivot_longer(-CNTRYID, names_to = "skill", values_to = "ml")

b_ols4 <- coef(ols_s4)
ref_country <- levels(train_s$CNTRYID)[1]

ols_slopes <- expand_grid(CNTRYID = country_codes, skill = cognitive) %>%
  mutate(ols_s4 = map2_dbl(CNTRYID, skill, function(cn, sk) {
    inter <- paste0("CNTRYID", cn, ":", sk)
    b_ols4[[sk]] + if (cn == ref_country) 0 else b_ols4[[inter]]
  }))

country_returns <- ml_slopes %>%
  left_join(ols_slopes, by = c("CNTRYID", "skill")) %>%
  mutate(
    country = unname(country_labels[CNTRYID]),
    skill   = dplyr::recode(skill, PVLIT1 = "Literacy", PVNUM1 = "Numeracy",
                     PVAPS1 = "Adaptive problem solving"),
    pct_ml     = 100 * (exp(ml) - 1),
    pct_ols_s4 = 100 * (exp(ols_s4) - 1)
  )

cat("\n=== RQ2: % WAGE CHANGE PER 1 SD OF SKILL, BY COUNTRY ===\n")
cat("ML1 (shrunk towards the common return) vs OLS S4 (no shrinkage)\n\n")
country_returns %>%
  mutate(cell = sprintf("%+.1f%% (%+.1f%%)", pct_ml, pct_ols_s4)) %>%
  select(country, skill, cell) %>%
  pivot_wider(names_from = skill, values_from = cell) %>%
  arrange(country) %>%
  print(n = Inf, width = Inf)
cat("\nCells: ML1 (OLS S4 in brackets)\n")

cat("\nCommon (fixed) slopes in ML1, % per 1 SD:\n")
print(round(100 * (exp(fixef(ml1$model)[cognitive]) - 1), 2))

write_csv(country_returns, "output/tables/table_3x_rq2_country_returns.csv")


# ---- D. BIG FIVE — ONE ADDITIONAL RANDOM SLOPE AT A TIME ----

bf_variants <- map(set_names(bigfive), function(trait) fit_ml(c(slopes, trait)))

bf_table <- imap_dfr(bf_variants, function(fit, trait) {
  vc  <- as.data.frame(VarCorr(fit$model))
  lrt <- anova(ml1$model, fit$model)
  tibble(
    trait        = trait,
    slope_var    = vc$vcov[vc$var1 == trait & !is.na(vc$var1)],
    slope_sd     = vc$sdcor[vc$var1 == trait & !is.na(vc$var1)],
    singular     = fit$singular,
    LR_chisq     = lrt$Chisq[2],
    p_value      = lrt$`Pr(>Chisq)`[2],
    # p-value halved for the boundary (Section 4.6.4); 1 if LR = 0
    p_boundary   = if_else(lrt$Chisq[2] < 1e-6, 1, lrt$`Pr(>Chisq)`[2] / 2),
    warnings     = if (length(fit$warnings)) paste(fit$warnings, collapse = "; ") else ""
  )
})

cat("\n=== BIG FIVE: RANDOM SLOPE ON ONE TRAIT, ADDED TO ML1 ===\n")
cat("slope_sd: SD across countries of the trait's return (log points per 1 SD)\n")
cat("LR test vs ML1; p_boundary: p-value halved for the boundary\n\n")
bf_table %>%
  mutate(across(c(slope_var, slope_sd, LR_chisq), ~ signif(.x, 3)),
         p_value = signif(p_value, 3),
         p_boundary = signif(p_boundary, 3)) %>%
  print(n = Inf, width = Inf)

write_csv(bf_table, "output/tables/table_3x_rq2_bigfive_slopes.csv")


# ---- E. RQ3 — TEST-SET PERFORMANCE AND BOOTSTRAP COMPARISON ----

pred_list <- list(
  OLS_S3 = unname(predict(ols_s3, newdata = test_s)),
  OLS_S4 = unname(predict(ols_s4, newdata = test_s)),
  ML0    = unname(predict(ml0$model, newdata = test_s)),
  ML1    = unname(predict(ml1$model, newdata = test_s))
)
for (trait in bigfive) {
  pred_list[[paste0("ML1_", trait)]] <- unname(predict(bf_variants[[trait]]$model, newdata = test_s))
}
pred_mat <- do.call(cbind, pred_list)

performance <- map_dfr(colnames(pred_mat), function(mname) {
  bind_cols(tibble(model = mname), oos(pred_mat[, mname], y_test))
})

cat("\n=== RQ3: TEST-SET PERFORMANCE ===\n")
performance %>% mutate(across(c(RMSE, MAE, R2), ~ round(.x, 4))) %>% print(n = Inf)

write_csv(performance, "output/tables/table_3x_rq3_multilevel_performance.csv")

comparisons <- bind_rows(
  tribble(
    ~label,                                         ~from,    ~to,
    "Random intercept vs country dummies (S3)",     "OLS_S3", "ML0",
    "Random slopes vs fixed interactions (S4)",     "OLS_S4", "ML1",
    "Random slopes vs random intercept only",       "ML0",    "ML1"
  ),
  tibble(label = paste("Random slope on", bigfive, "added to ML1"),
         from  = "ML1",
         to    = paste0("ML1_", bigfive))
)

boot_ml <- paired_bootstrap(pred_mat, y_test, test_s$CNTRYID, comparisons)

cat("\n=== RQ3: DIFFERENCES WITH 95% PAIRED-BOOTSTRAP CI (test set, B = 2000) ===\n")
cat("dR2 > 0 and dRMSE < 0 mean that 'to' predicts better than 'from'\n\n")
boot_ml %>%
  filter(metric %in% c("dR2", "dRMSE")) %>%
  mutate(flag  = case_when(improves ~ " better", worsens ~ " worse", TRUE ~ ""),
         value = sprintf("%+.4f [%+.4f, %+.4f]%s", estimate, ci_low, ci_high, flag)) %>%
  select(comparison, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  print(n = Inf, width = Inf)
cat("\nbetter / worse = 95% interval excludes zero in that direction;",
    "no label = difference not distinguishable from zero\n")

write_csv(boot_ml, "output/tables/table_3x_rq3_multilevel_bootstrap.csv")


# ---- F. DIAGNOSTIC CHECKS ON ML1 (Sections 4.6.6 and 5.3) ----

res1 <- resid(ml1$model)

# Level-1 residuals: skewness and excess kurtosis (as for OLS in 03)
cat("\n=== DIAGNOSTIC CHECKS ON ML1 ===\n")
cat("Skewness:", round(mean((res1 - mean(res1))^3) / sd(res1)^3, 3),
    "| Excess kurtosis:", round(mean((res1 - mean(res1))^4) / sd(res1)^4 - 3, 3), "\n")

# Level-1 residuals: standard deviation by country
tibble(country = unname(country_labels[as.character(train_s$CNTRYID)]), res = res1) %>%
  group_by(country) %>%
  summarise(sd_res = round(sd(res), 3)) %>%
  arrange(sd_res) %>%
  print(n = Inf)

# Level 2: predicted country deviations (12 values per random effect)
print(round(ranef(ml1$model)$CNTRYID, 3))


# END OF SCRIPT 07 ----
