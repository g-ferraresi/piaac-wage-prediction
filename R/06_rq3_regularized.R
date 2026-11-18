# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 06_rq3_regularized.R
# RQ3 — regularized estimators (Ridge, LASSO, Elastic Net)
# compared with OLS on the same specifications (Sections 4.5, 5.4)
#
# Reads train.rds / test.rds (01_data_cleaning.R) and the variable
# blocks in 00_variable_blocks.R.
#
# Design:
#   S3 = controls + cognitive skills + Big Five (+ country dummies)
#   S4 = S3 + interactions country x cognitive skill
# Every estimator uses the same information (same design matrix):
# differences in test-set accuracy come from the estimation method.
#   - Country dummies are NOT penalized (penalty.factor = 0): all
#     estimators keep the full country information, as OLS does.
#   - In S4 the country x skill interactions ARE penalized: Ridge /
#     LASSO / Elastic Net shrink the country-specific skill returns
#     towards the common return. This is the "fixed-effects" analogue
#     of the multilevel random slopes (07_multilevel.R).
#
# The script keeps the TWO ROLES of regularization separate:
#   PART A — regularization as a COMPETING ESTIMATOR (RQ3):
#            test-set RMSE, MAE, R² vs OLS, with paired-bootstrap CIs
#   PART B — LASSO as a VARIABLE-SELECTION tool: which variables
#            are set to zero
#
# Tuning (training set only, the test set is never used for tuning):
#   - 10-fold cross-validation, folds stratified by country, the SAME
#     folds for every alpha and specification
#   - Ridge: alpha = 0; LASSO: alpha = 1; Elastic Net: best alpha in
#     {0.1, ..., 0.9} by cross-validated MSE
#   - lambda.min (lowest CV error) for prediction; lambda.1se (the
#     sparsest model within 1 SE of the minimum) reported for LASSO
#     selection as the more conservative alternative
#
# Cognitive skills (PV1) are standardized with the TRAINING mean and
# SD. For OLS this changes nothing in the predictions (checked below);
# for the S4 interactions it makes the country main effects refer to a
# worker of average skill, and it is the scaling that the multilevel
# models also use. Big Five are already within-country z-scores.
# R² is the out-of-sample R² used in 03-05.


# Run from the repository root (see README).
# ---- PACKAGES ----

library(tidyverse)
library(glmnet)     # cv.glmnet()


# ---- LOAD DATA AND VARIABLE BLOCKS ----

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
test  <- readRDS("data/processed/test.rds")
check_blocks(train)


# ---- A. STANDARDIZE COGNITIVE SKILLS (training mean and SD) ----

skill_mean <- setNames(map_dbl(cognitive, ~ mean(train[[.x]])), cognitive)
skill_sd   <- setNames(map_dbl(cognitive, ~ sd(train[[.x]])),   cognitive)

standardize_skills <- function(data) {
  for (v in cognitive) data[[v]] <- (data[[v]] - skill_mean[v]) / skill_sd[v]
  data
}

train_s <- standardize_skills(train)
test_s  <- standardize_skills(test)

y_train <- train_s$LOG_WAGE
y_test  <- test_s$LOG_WAGE


# ---- B. SPECIFICATIONS AND DESIGN MATRICES ----

interactions <- paste0("CNTRYID:", cognitive)

formulas <- list(
  S3 = make_formula(c(controls, cognitive, bigfive)),
  S4 = make_formula(c(controls, cognitive, bigfive, interactions))
)

# Design matrix without the intercept column (glmnet adds its own),
# plus the name of the original variable behind each column
build_x <- function(f, data) {
  mm <- model.matrix(f, data)
  term_labels <- attr(terms(f), "term.labels")
  list(X = mm[, -1, drop = FALSE],
       col_term = term_labels[attr(mm, "assign")[-1]])
}

design_train <- map(formulas, ~ build_x(.x, train_s))
design_test  <- map(formulas, ~ build_x(.x, test_s))

# Train and test must have exactly the same columns in the same order
walk(names(formulas), function(s) {
  stopifnot(identical(colnames(design_train[[s]]$X), colnames(design_test[[s]]$X)))
})

# Penalty factors: 0 for the country dummies (main effects only),
# 1 for everything else, interactions included
penalty <- map(design_train, function(d) {
  ifelse(grepl("^CNTRYID[0-9]+$", colnames(d$X)), 0, 1)
})

cat("Columns in the design matrix: S3 =", ncol(design_train$S3$X),
    "| S4 =", ncol(design_train$S4$X), "\n")
cat("Unpenalized columns (country dummies):", sum(penalty$S3 == 0), "\n\n")


# ---- C. HELPERS ----

oos <- function(pred, y) {
  tibble(
    RMSE = sqrt(mean((y - pred)^2)),
    MAE  = mean(abs(y - pred)),
    R2   = 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  )
}

# Paired bootstrap on the test set, stratified by country (same method
# as 04_rq1_specifications.R, section C2). pred_mat: one column per
# model; comparisons: tibble with columns label, from, to.
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
  # Point estimates on the full test set
  full <- apply(pred_mat, 2, function(p) unlist(oos(p, y)))   # 3 x models
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


# ---- D. OLS BENCHMARKS ----

ols <- map(formulas, ~ lm(.x, data = train_s))

# Check: standardizing the skills does not change the OLS predictions
# (S3 must reproduce the 03 baseline, estimated on unstandardized data)
ols_raw_s3 <- lm(formulas$S3, data = train)
stopifnot(isTRUE(all.equal(unname(predict(ols_raw_s3, newdata = test)),
                           unname(predict(ols$S3, newdata = test_s)))))
cat("Check passed: OLS S3 reproduces the 03 baseline.\n\n")



# ---- E. REGULARIZATION AS A COMPETING ESTIMATOR (RQ3) ----

# ---- E1. Cross-validation folds (stratified by country) ----

set.seed(2026)
foldid <- integer(nrow(train_s))
for (idx in split(seq_len(nrow(train_s)), train_s$CNTRYID)) {
  foldid[idx] <- sample(rep_len(1:10, length(idx)))
}


# ---- E2. Tuning: alpha grid x lambda path ----

alphas <- c(0, seq(0.1, 0.9, by = 0.1), 1)

cv_fits <- map(names(formulas), function(s) {
  fits <- map(alphas, function(a) {
    cv.glmnet(design_train[[s]]$X, y_train,
              alpha = a, foldid = foldid,
              penalty.factor = penalty[[s]],
              standardize = TRUE)
  })
  cat("Cross-validation done for", s, "\n")
  set_names(fits, paste0("alpha_", alphas))
}) %>% set_names(names(formulas))

alpha_grid <- imap_dfr(cv_fits, function(fits, s) {
  tibble(
    spec       = s,
    alpha      = alphas,
    cv_mse_min = map_dbl(fits, ~ min(.x$cvm)),
    lambda_min = map_dbl(fits, "lambda.min"),
    lambda_1se = map_dbl(fits, "lambda.1se"),
    nonzero_min = map_int(fits, ~ as.integer(.x$nzero[which(.x$lambda == .x$lambda.min)]))
  )
})

cat("\n=== A2. CROSS-VALIDATED MSE BY ALPHA (training set) ===\n")
alpha_grid %>%
  mutate(across(c(cv_mse_min, lambda_min, lambda_1se), ~ signif(.x, 4))) %>%
  print(n = Inf)

write_csv(alpha_grid, "output/tables/table_3x_rq3_alpha_grid.csv")

# Selected models: Ridge (alpha 0), LASSO (alpha 1), Elastic Net (best
# alpha strictly between 0 and 1)
best_en_alpha <- alpha_grid %>%
  filter(alpha > 0, alpha < 1) %>%
  group_by(spec) %>%
  slice_min(cv_mse_min, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(spec, alpha) %>%
  deframe()

cat("\nElastic Net alpha selected by CV: S3 =", best_en_alpha["S3"],
    "| S4 =", best_en_alpha["S4"], "\n")

selected <- map(names(formulas), function(s) {
  list(
    Ridge = cv_fits[[s]][["alpha_0"]],
    LASSO = cv_fits[[s]][["alpha_1"]],
    EN    = cv_fits[[s]][[paste0("alpha_", best_en_alpha[s])]]
  )
}) %>% set_names(names(formulas))


# ---- E3. Test-set performance ----

pred_list <- list()
for (s in names(formulas)) {
  pred_list[[paste0("OLS_", s)]] <- unname(predict(ols[[s]], newdata = test_s))
  for (est in c("Ridge", "LASSO", "EN")) {
    pred_list[[paste0(est, "_", s)]] <-
      as.numeric(predict(selected[[s]][[est]], newx = design_test[[s]]$X, s = "lambda.min"))
  }
}
pred_mat <- do.call(cbind, pred_list)

model_info <- tibble(model = colnames(pred_mat)) %>%
  separate(model, into = c("estimator", "spec"), sep = "_", remove = FALSE) %>%
  mutate(
    alpha = case_when(estimator == "Ridge" ~ 0,
                      estimator == "LASSO" ~ 1,
                      estimator == "EN"    ~ unname(best_en_alpha[spec]),
                      TRUE ~ NA_real_),
    lambda = map2_dbl(estimator, spec, function(e, s) {
      if (e == "OLS") NA_real_ else selected[[s]][[e]]$lambda.min
    }),
    n_nonzero = map2_int(estimator, spec, function(e, s) {
      if (e == "OLS") length(coef(ols[[s]])) - 1L
      else {
        f <- selected[[s]][[e]]
        as.integer(f$nzero[which(f$lambda == f$lambda.min)])
      }
    })
  )

performance <- bind_cols(
  model_info,
  map_dfr(seq_len(ncol(pred_mat)), ~ oos(pred_mat[, .x], y_test))
)

cat("\n=== A3. TEST-SET PERFORMANCE BY ESTIMATOR AND SPECIFICATION ===\n")
performance %>%
  mutate(across(c(RMSE, MAE, R2), ~ round(.x, 4)),
         lambda = signif(lambda, 3)) %>%
  print(n = Inf)

write_csv(performance, "output/tables/table_3x_rq3_regularized_performance.csv")


# ---- E4. Paired bootstrap: are the differences real? ----

comparisons <- tribble(
  ~label,                               ~from,       ~to,
  "Ridge vs OLS (S3)",                  "OLS_S3",    "Ridge_S3",
  "LASSO vs OLS (S3)",                  "OLS_S3",    "LASSO_S3",
  "Elastic Net vs OLS (S3)",            "OLS_S3",    "EN_S3",
  "Ridge vs OLS (S4)",                  "OLS_S4",    "Ridge_S4",
  "LASSO vs OLS (S4)",                  "OLS_S4",    "LASSO_S4",
  "Elastic Net vs OLS (S4)",            "OLS_S4",    "EN_S4",
  "Country x skill interactions (OLS)", "OLS_S3",    "OLS_S4",
  "Country x skill interactions (Ridge)", "Ridge_S3", "Ridge_S4",
  "Country x skill interactions (LASSO)", "LASSO_S3", "LASSO_S4",
  "Country x skill interactions (EN)",  "EN_S3",     "EN_S4"
)

boot_rq3 <- paired_bootstrap(pred_mat, y_test, test_s$CNTRYID, comparisons)

cat("\n=== A4. DIFFERENCES WITH 95% PAIRED-BOOTSTRAP CI (test set, B = 2000) ===\n")
cat("dR2 > 0 and dRMSE < 0 mean that 'to' predicts better than 'from'\n\n")
boot_rq3 %>%
  filter(metric %in% c("dR2", "dRMSE")) %>%
  mutate(flag  = case_when(improves ~ " better", worsens ~ " worse", TRUE ~ ""),
         value = sprintf("%+.4f [%+.4f, %+.4f]%s", estimate, ci_low, ci_high, flag)) %>%
  select(comparison, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  print(n = Inf, width = Inf)
cat("\nbetter / worse = 95% interval excludes zero in that direction;",
    "no label = difference not distinguishable from zero\n")

write_csv(boot_rq3, "output/tables/table_3x_rq3_bootstrap.csv")


# ---- E5. Shrinkage of the key coefficients (S3) ----
# Coefficients of cognitive skills (per 1 SD, standardized above) and
# Big Five (per 1 within-country SD), as % wage change, for OLS and
# the three regularized estimators at lambda.min.

key_vars <- c(cognitive, bigfive)

coef_s3 <- tibble(term = key_vars,
                  OLS = unname(coef(ols$S3)[key_vars]))
for (est in c("Ridge", "LASSO", "EN")) {
  b <- coef(selected$S3[[est]], s = "lambda.min")
  coef_s3[[est]] <- as.numeric(b[key_vars, 1])
}

cat("\n=== A5. KEY COEFFICIENTS IN S3 — % WAGE CHANGE PER 1 SD ===\n")
coef_s3 %>%
  mutate(across(c(OLS, Ridge, LASSO, EN), ~ sprintf("%+.1f%%", 100 * (exp(.x) - 1)))) %>%
  print(n = Inf)

write_csv(coef_s3, "output/tables/table_3x_rq3_key_coefficients.csv")



# ---- F.  LASSO AS A VARIABLE-SELECTION TOOL ----
# For each original variable (a factor spans several dummy columns),
# how many of its columns LASSO keeps (non-zero) at lambda.min and at
# the more conservative lambda.1se. Specification S3.

lasso_s3 <- selected$S3$LASSO
cols     <- colnames(design_train$S3$X)

b_min <- as.numeric(coef(lasso_s3, s = "lambda.min"))[-1]
b_1se <- as.numeric(coef(lasso_s3, s = "lambda.1se"))[-1]

lasso_selection <- tibble(
  variable = design_train$S3$col_term,
  column   = cols,
  coef_min = b_min,
  coef_1se = b_1se
) %>%
  group_by(variable) %>%
  summarise(
    n_columns   = n(),
    kept_min    = sum(coef_min != 0),
    kept_1se    = sum(coef_1se != 0),
    .groups = "drop"
  ) %>%
  mutate(
    block = case_when(
      variable %in% cognitive         ~ "Cognitive skills",
      variable %in% bigfive           ~ "Big Five",
      variable %in% job_controls      ~ "Job controls",
      variable == "CNTRYID"           ~ "Country (unpenalized)",
      TRUE                            ~ "Personal controls"
    ),
    status_1se = case_when(kept_1se == 0 ~ "dropped",
                           kept_1se < n_columns ~ "partly kept",
                           TRUE ~ "kept")
  ) %>%
  arrange(block, variable)

cat("\n=== B. LASSO VARIABLE SELECTION (S3) ===\n")
cat("lambda.min =", signif(lasso_s3$lambda.min, 3),
    "| lambda.1se =", signif(lasso_s3$lambda.1se, 3), "\n")
cat("Columns kept: at lambda.min", sum(b_min != 0), "of", length(b_min),
    "| at lambda.1se", sum(b_1se != 0), "of", length(b_1se), "\n\n")
print(lasso_selection, n = Inf)

write_csv(lasso_selection, "output/tables/table_3x_rq3_lasso_selection.csv")


# END OF SCRIPT 06 ----
