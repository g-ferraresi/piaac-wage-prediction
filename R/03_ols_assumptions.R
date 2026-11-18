# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 03_ols_assumptions.R
#
# OLS baseline (specification S3) and diagnostic checks
# (Sections 4.3 and 5.1).
#   A1. Linearity          - residuals vs fitted values
#   A2. Homoscedasticity   - Breusch-Pagan test; HC3 standard errors
#   A3. Normality          - Q-Q plot, skewness and excess kurtosis,
#                            Kolmogorov-Smirnov test
#   A4. Multicollinearity  - variance inflation factors
#   A5. Intraclass correlation coefficient (null multilevel model)

# Run from the repository root (see README).
# ---- PACKAGES ----

library(tidyverse)
library(lmtest)   # bptest()
library(car)      # vif()
library(lme4)     # lmer() — null model for ICC


# ---- LOAD CLEANED DATA (from 01_data_cleaning.R) ----

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
test  <- readRDS("data/processed/test.rds")

check_blocks(train)

cat("Training set:", nrow(train), "observations,", nlevels(train$CNTRYID), "countries\n")
cat("Test set:    ", nrow(test),  "observations\n\n")


# ---- A. OLS BASELINE ----

# Specification S3: controls + cognitive skills (PV1) + Big Five.
# CNTRYID is a factor, so it enters as 11 country dummies.

ols_model <- lm(make_formula(c(controls, cognitive, bigfive)), data = train)

ols_pred <- predict(ols_model, newdata = test)
y_test   <- test$LOG_WAGE

ols_rmse <- sqrt(mean((y_test - ols_pred)^2))
ols_mae  <- mean(abs(y_test - ols_pred))

# Out-of-sample R² = 1 - SSE/SST on the test set. Unlike the squared
# correlation between observed and predicted values, it also
# penalizes systematically biased predictions.
ols_r2 <- 1 - sum((y_test - ols_pred)^2) / sum((y_test - mean(y_test))^2)

cat("=== OLS BASELINE PERFORMANCE (test set) ===\n")
cat("RMSE:", round(ols_rmse, 4), "\n")
cat("MAE: ", round(ols_mae,  4), "\n")
cat("R²:  ", round(ols_r2,   4), "\n\n")


# ---- A2. ASSUMPTION CHECKS ----


# A1. LINEARITY — Residuals vs Fitted Values
# Visual check: residuals should scatter randomly around zero
# with no systematic pattern. A curved or funnel shape signals
# non-linearity or heteroscedasticity.

ols_resid  <- residuals(ols_model)
ols_fitted <- fitted(ols_model)

p_resid_fitted <- data.frame(fitted = ols_fitted, resid = ols_resid) %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.08, size = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, color = "firebrick",
             linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE,
              color = "darkorange", linewidth = 0.8) +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Orange line = LOESS smoother (should be flat at zero)",
       x = "Fitted values (log wage)", y = "Residuals") +
  theme_minimal()
print(p_resid_fitted)
ggsave("output/figures/fig_3x_residuals_vs_fitted.png", p_resid_fitted, width = 7, height = 4.5, dpi = 300)


# A2. HOMOSCEDASTICITY — Breusch-Pagan Test
# H0: residuals are homoscedastic (constant variance).
# A significant result (p < 0.05) indicates heteroscedasticity.
# Note: with large N, even minor deviations become significant —
# inspect the scale-location plot alongside the test result.

bp_test <- bptest(ols_model)

cat("=== A2. BREUSCH-PAGAN TEST (Homoscedasticity) ===\n")
print(bp_test)
cat("\n")

p_scale_location <- data.frame(
    fitted = ols_fitted,
    sqrt_std_resid = sqrt(abs(scale(ols_resid)))
  ) %>%
  ggplot(aes(x = fitted, y = sqrt_std_resid)) +
  geom_point(alpha = 0.08, size = 0.5, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "darkorange", linewidth = 0.8) +
  labs(title = "Scale-Location Plot",
       subtitle = "Flat orange line = homoscedasticity",
       x = "Fitted values (log wage)",
       y = expression(sqrt("|Standardised residuals|"))) +
  theme_minimal()
print(p_scale_location)
ggsave("output/figures/fig_3x_scale_location.png", p_scale_location, width = 7, height = 4.5, dpi = 300)


# A3. NORMALITY OF RESIDUALS
# Shapiro-Wilk is NOT used here: it is invalid/unreliable for
# N > 5000 (our training set has ~20,000+ observations), and with
# very large N even trivial deviations from normality become
# "significant". Instead we rely on:
#   (i)   Q-Q plot (visual — the standard diagnostic in applied work)
#   (ii)  Skewness and excess kurtosis (descriptive — should be
#         close to 0 for approximate normality)
#   (iii) Kolmogorov-Smirnov test against a normal distribution
#         with the residuals' own mean/sd (informative but, like
#         Shapiro-Wilk, over-powered at large N — interpret
#         alongside the Q-Q plot, not in isolation)

cat("=== A3. NORMALITY OF RESIDUALS ===\n")

# (i) Q-Q plot
p_qq <- data.frame(resid = ols_resid) %>%
  ggplot(aes(sample = resid)) +
  stat_qq(alpha = 0.3, size = 0.6, color = "steelblue") +
  stat_qq_line(color = "firebrick", linewidth = 0.8) +
  labs(title = "Normal Q-Q Plot of OLS Residuals",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()
print(p_qq)
ggsave("output/figures/fig_3x_qq_residuals.png", p_qq, width = 6, height = 5, dpi = 300)

# (ii) Skewness and excess kurtosis (manual — avoids adding a
# dependency just for two descriptive moments)
resid_skewness <- mean((ols_resid - mean(ols_resid))^3) / sd(ols_resid)^3
resid_kurtosis <- mean((ols_resid - mean(ols_resid))^4) / sd(ols_resid)^4 - 3  # excess kurtosis

cat("Skewness:        ", round(resid_skewness, 3), " (0 = symmetric)\n")
cat("Excess kurtosis: ", round(resid_kurtosis, 3), " (0 = normal-tailed; >0 = heavy-tailed)\n\n")

# (iii) Kolmogorov-Smirnov test against N(mean, sd) of the residuals
ks_test <- ks.test(ols_resid, "pnorm", mean = mean(ols_resid), sd = sd(ols_resid))
cat("Kolmogorov-Smirnov test vs. Normal(resid mean, resid sd):\n")
print(ks_test)
cat("NOTE: with N =", length(ols_resid), ", this test will almost certainly\n",
    "reject exact normality even for minor deviations. Judge practical\n",
    "relevance from the Q-Q plot and skewness/kurtosis above, not from\n",
    "the p-value alone.\n\n")


# A2b. HETEROSCEDASTICITY-CONSISTENT STANDARD ERRORS (HC3)
# The Breusch-Pagan test rejects homoscedasticity. OLS coefficients
# remain unbiased, but the classical standard errors are unreliable:
# they are replaced by the HC3 estimator (MacKinnon & White, 1985;
# Long & Ervin, 2000).

library(sandwich)

cat("=== A2b. HC3-ROBUST STANDARD ERRORS ===\n")

vcov_hc3 <- vcovHC(ols_model, type = "HC3")

# Classical vs. robust standard errors, side by side, for comparison
se_classical <- sqrt(diag(vcov(ols_model)))
se_robust    <- sqrt(diag(vcov_hc3))

se_comparison <- data.frame(
  term          = names(se_classical),
  estimate      = coef(ols_model),
  se_classical  = round(se_classical, 4),
  se_hc3        = round(se_robust, 4),
  pct_change    = round((se_robust - se_classical) / se_classical * 100, 1)
) %>%
  arrange(desc(abs(pct_change)))

cat("\nStandard errors: classical (OLS) vs. HC3-robust, sorted by %% change\n")
cat("(large positive %% change = classical SEs were understating uncertainty)\n\n")
print(head(se_comparison, 15), row.names = FALSE)

write_csv(se_comparison, "output/tables/table_3x_se_classical_vs_hc3.csv")

# Regression table with HC3 standard errors (coefficients unchanged,
# only SE/t/p recomputed)
robust_coeftest <- coeftest(ols_model, vcov. = vcov_hc3)

cat("\n=== FULL MODEL WITH HC3-ROBUST INFERENCE ===\n")
cat("(coefficients identical to classical OLS above; only SE/t/p differ)\n\n")
print(robust_coeftest)

# Summary: how many coefficients change significance status
# (classical p < 0.05 vs. robust p < 0.05) — useful one-line check
# for whether robust SEs would substantively change the story.
classical_p <- summary(ols_model)$coefficients[, "Pr(>|t|)"]
robust_p    <- robust_coeftest[, "Pr(>|t|)"]

n_flip <- sum((classical_p < 0.05) != (robust_p < 0.05))
cat("\nCoefficients whose significance (at 5%) flips between classical",
    "and HC3-robust SEs:", n_flip, "out of", length(classical_p), "\n")
cat("(0 or very few flips = the heteroscedasticity found by Breusch-Pagan,\n",
    "while statistically significant, has limited practical consequence\n",
    "for which coefficients you'd call 'significant' in the thesis.)\n\n")


# A4. MULTICOLLINEARITY — Variance Inflation Factor
# VIF > 10 signals problematic collinearity; VIF > 5 worth inspecting.
# Dummy-expanded factors are assessed via Generalised VIF (GVIF),
# which car::vif() returns automatically for factor predictors.

vif_results <- vif(ols_model)

cat("=== A4. VIF — MULTICOLLINEARITY CHECK ===\n")
print(round(vif_results, 2))
cat("\n")

vif_df <- as.data.frame(vif_results)

# car::vif() returns either a plain vector (all-numeric predictors)
# or a matrix with GVIF/Df/GVIF^(1/(2*Df)) columns (mix with factors)
if (ncol(vif_df) == 1) {
  vif_df <- vif_df %>%
    rownames_to_column("variable") %>%
    rename(VIF = 2)
} else {
  vif_df <- vif_df %>%
    rownames_to_column("variable") %>%
    rename(GVIF = 2, Df = 3, GVIF_adj = 4)
  vif_df$VIF <- vif_df$GVIF_adj^2   # square to get a VIF-comparable scale
}

vif_flagged <- vif_df %>% filter(VIF > 10)

if (nrow(vif_flagged) > 0) {
  cat("Variables with VIF > 10 (potentially problematic):\n")
  print(vif_flagged)
} else {
  cat("No variable exceeds the VIF threshold of 10.\n")
}
cat("\n")

write_csv(vif_df, "output/tables/table_3x_vif.csv")


# A5. INTRACLASS CORRELATION COEFFICIENT (ICC)
# Fits a null (random-intercept-only, no predictors) multilevel
# model of LOG_WAGE on country, to quantify how much of the total
# wage variance lies BETWEEN countries before any predictor is
# added. This is the standard first step before any multilevel
# modelling (Raudenbush & Bryk, 2002) and directly motivates the
# move beyond pooled OLS (RQ2/RQ3).
#
# ICC = between-country variance / (between-country + within-country variance)

null_model <- lmer(LOG_WAGE ~ 1 + (1 | CNTRYID), data = train)

cat("=== A5. NULL MULTILEVEL MODEL (for ICC) ===\n")
print(summary(null_model))

var_components <- as.data.frame(VarCorr(null_model))
var_between <- var_components$vcov[var_components$grp == "CNTRYID"]
var_within  <- var_components$vcov[var_components$grp == "Residual"]

icc <- var_between / (var_between + var_within)

cat("\nVariance between countries:", round(var_between, 4), "\n")
cat("Variance within countries: ", round(var_within, 4), "\n")
cat("ICC:", round(icc, 4),
    "-> approximately", paste0(round(icc * 100, 1), "%"),
    "of total wage variance lies between countries\n")
cat("(Share of wage variance located between countries in the null model,\n",
    "i.e. before any predictor is added. It motivates modelling the country\n",
    "level explicitly: country dummies in OLS, random effects in the\n",
    "multilevel models.)\n\n")


# SUMMARY: ALL ASSUMPTION CHECK RESULTS

cat("=== SUMMARY: ASSUMPTION CHECKING RESULTS ===\n")
cat("A1. Linearity:        see fig_3x_residuals_vs_fitted.png\n")
cat("A2. Homoscedasticity: Breusch-Pagan p =", format.pval(bp_test$p.value, digits = 3),
    "| see fig_3x_scale_location.png\n")
cat("A3. Normality:        skewness =", round(resid_skewness, 3),
    ", excess kurtosis =", round(resid_kurtosis, 3),
    "| KS p =", format.pval(ks_test$p.value, digits = 3),
    "| see fig_3x_qq_residuals.png\n")
cat("A4. Multicollinearity:", nrow(vif_flagged), "variable(s) with VIF > 10\n")
cat("A5. ICC:              ", round(icc, 4), "(", round(icc * 100, 1), "% between-country variance)\n")


# END OF SCRIPT 03 ----