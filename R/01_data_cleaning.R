# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 01_data_cleaning.R
#
# Data preparation: variable selection, recoding of missing codes,
# country selection, PPP conversion, sample restriction, trimming,
# listwise deletion and training/test split (Chapter 3).
# All the other scripts read the objects saved at the end of this
# script (df_model.rds, train.rds, test.rds).
#
# Country selection:
#   - countries in which hourly earnings (EARNHRC2) are entirely
#     missing are excluded;
#   - Estonia is excluded: four variables (EDCAT8_TC1, C2_Q10,
#     D2_Q10, D2_Q11) are entirely missing;
#   - Spain is retained by dropping D2_Q10 (contract type), which
#     is entirely missing for Spain.
# Final sample: 12 countries (Chile, Croatia, Czechia, France,
# Ireland, Israel, Italy, Latvia, Lithuania, New Zealand, Poland,
# Spain); N = 26,336 (training set 21,074, test set 5,262).


# ---- LIBRARIES AND DATASET ----

library(tidyverse)
library(caret)

# Run from the repository root (see README).
# The PIAAC Cycle 2 Public Use File is not included in this repository:
# see data/README.md for how to download it and create PIAAC_Cycle2.rds.
PIAAC_Cycle2 <- readRDS("data/raw/PIAAC_Cycle2.rds")
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ---- COUNTRIES IN THE PUBLIC USE FILE ----

# Number of countries in the original dataset
n_distinct(PIAAC_Cycle2$CNTRYID)

# Observations and % of missing wages by country
PIAAC_Cycle2 %>%
  mutate(EARNHRC2 = suppressWarnings(as.numeric(EARNHRC2))) %>%
  group_by(CNTRYID) %>%
  summarise(
    n_obs        = n(),
    pct_wage_na  = round(mean(is.na(EARNHRC2)) * 100, 1)
  ) %>%
  arrange(CNTRYID) %>%
  print(n = Inf)

# Number of countries with and without the wage variable
PIAAC_Cycle2 %>%
  mutate(EARNHRC2 = suppressWarnings(as.numeric(EARNHRC2))) %>%
  group_by(CNTRYID) %>%
  summarise(pct_wage_na = mean(is.na(EARNHRC2)) * 100) %>%
  summarise(
    tot     = n(),
    with_wage = sum(pct_wage_na < 100),
    no_wage = sum(pct_wage_na == 100)
  )

# Countries with the wage variable
PIAAC_Cycle2 %>%
  mutate(EARNHRC2 = suppressWarnings(as.numeric(EARNHRC2))) %>%
  group_by(CNTRYID) %>%
  summarise(
    n_obs       = n(),
    pct_wage_na = round(mean(is.na(EARNHRC2)) * 100, 1)
  ) %>%
  filter(pct_wage_na < 100) %>%
  arrange(CNTRYID) %>%
  print(n = Inf)

# Variable blocks, plausible-value and weight names
source("R/00_variable_blocks.R")

# ---- A. VARIABLE SELECTION ----
# Extended Mincer (1974) framework, with measured cognitive skills
# (Hanushek et al., 2015) and Big Five personality traits (Mueller
# & Plug, 2004; Nyhus & Pons, 2005).
# Not used: AGE_R and AGEG5LFS (replaced by AGEG10LFS, available for
# all countries), D2_Q10 (contract type, entirely missing for Spain)
# and EARNHRDCLC2 (wage deciles).

df <- PIAAC_Cycle2 %>%
  select(
    # Country identifier
    CNTRYID,
    # Dependent variable
    EARNHRC2,        # Hourly earnings in local currency (continuous)
    # Demographics
    AGEG10LFS,       # Age group (10-year bands, available for all countries)
    GENDER_R,        # Gender
    # Education and training
    EDCAT8_TC1,      # Educational attainment (8 ISCED categories)
    YRSGETC2,        # Years of education required to get current job
    NFE12C2,         # Non-formal education in last 12 months (0/1)
    # Work experience
    C2_Q10,          # Total years of paid work over lifetime
    # PAIDWORK12 is not included: it is constant (= 1) once the sample
    # is restricted to positive wage earners.
    # Measured cognitive skills: all 10 plausible values for
    # literacy, numeracy and adaptive problem solving (PV1 is used
    # in the main models; PV2-PV10 for Rubin's rules)
    all_of(pv_all),
    # Final full sample weight (weighted robustness check only)
    all_of(weight_var),
    # Occupation and job characteristics
    # ISCO1C is not included: 22.2% of the sample has missing codes
    # (995-999) and would be lost with listwise deletion.
    ISIC1C,          # Economic sector (ISIC Rev.4, 1-digit)
    ISCOSKIL4,       # Skill level required by occupation (4 categories)
    # D2_Q04 (employee vs. self-employed) is not included: it is
    # constant in the final sample.
    D2_Q07a,         # Firm size (number of employees)
    D2_Q11,          # Weekly working hours
    # Personality — Big Five: 3-item domain scores common to BFI-2-S
    # and BFI-2-XS (instrument differs by country), provided by the
    # OECD as within-country z-scores
    AGRE, CONS, EMOS, EXTR, OPEM,
    # Migration and language background
    NATIVELANG,      # Native language same as test language (0/1)
    IMGENC2,         # First/second generation immigrant
    # Family background (origin)
    PAREDC2,         # Highest parental education level
    # Current family characteristics
    PARTNER          # Living with spouse or partner (0/1)
  ) %>%
  mutate(
    # Generic missing codes (international codebook) 
    #   .  = Not administered
    #   .d = Don't know
    #   .n = Not stated or inferred
    #   .r = Refused
    #   .u = Unknown
    #   .v = Valid skip
    across(everything(), ~na_if(.x, ".")),
    across(everything(), ~na_if(.x, ".d")),
    across(everything(), ~na_if(.x, ".n")),
    across(everything(), ~na_if(.x, ".r")),
    across(everything(), ~na_if(.x, ".u")),
    across(everything(), ~na_if(.x, ".v")),

    # Variable-specific meaning of ".a"
    #   - AGEG10LFS: ".a" = younger than 16 (outside the target
    #     population) -> NA
    #   - ISCOSKIL4: ".a" = has not worked in the past 5 years -> NA
    #   - IMGENC2:   ".a" = non-immigrant with one foreign-born parent,
    #     a substantive category: kept as level "4"
    AGEG10LFS = na_if(AGEG10LFS, ".a"),
    ISCOSKIL4 = na_if(ISCOSKIL4, ".a"),
    # IMGENC2: ".a" relabelled as level "4"
    IMGENC2 = if_else(IMGENC2 == ".a", "4", IMGENC2),

    # AGEG10LFS ".b" = ">65" (also out of the 16-65 target population)
    AGEG10LFS = na_if(AGEG10LFS, ".b"),

    # ISIC1C: the numeric codes 995-999 are missing codes (the
    # valid categories are the letters A-U)
    ISIC1C = na_if(ISIC1C, "995"),
    ISIC1C = na_if(ISIC1C, "996"),
    ISIC1C = na_if(ISIC1C, "997"),
    ISIC1C = na_if(ISIC1C, "998"),
    ISIC1C = na_if(ISIC1C, "999")
  ) %>%
  mutate(
    across(c(YRSGETC2, C2_Q10, D2_Q11,
             all_of(pv_all), all_of(weight_var),
             AGRE, CONS, EMOS, EXTR, OPEM,
             EARNHRC2), as.numeric),
    across(c(CNTRYID, AGEG10LFS, GENDER_R, EDCAT8_TC1,
             ISIC1C, ISCOSKIL4, D2_Q07a,
             NFE12C2, NATIVELANG, IMGENC2,
             PAREDC2, PARTNER), as.factor)
  )


# ---- B. FINAL COUNTRY LIST + FILTER ----

final_countries <- c(
  "152", # Chile
  "191", # Croatia
  "203", # Czechia
  "250", # France
  "372", # Ireland
  "376", # Israel
  "380", # Italy
  "428", # Latvia
  "440", # Lithuania
  "554", # New Zealand
  "616", # Poland
  "724"  # Spain
)

df <- df %>%
  filter(CNTRYID %in% final_countries) %>%
  mutate(CNTRYID = droplevels(CNTRYID))

cat("Countries retained:", nlevels(df$CNTRYID), "\n")
cat("Observations after country filter:", nrow(df), "\n\n")


# ---- C. PPP CONVERSION TO 2022 USD ----
# Source: OECD Data Explorer (DSD_NAMAIN10, DF_TABLE4, PPP_B1GQ,
# reference year 2022): local currency units per USD.
# Croatia: wages are recorded in HRK; the factor combines the fixed
# HRK/EUR conversion rate (7.53450, Council Regulation (EU)
# 2022/1207) with the EUR/USD PPP.

ppp_2022 <- data.frame(
  CNTRYID = c("152", "191", "203", "250",
              "372", "376", "380", "428", "440",
              "554", "616", "724"),
  country = c("Chile", "Croatia", "Czechia", "France",
              "Ireland", "Israel", "Italy", "Latvia", "Lithuania",
              "New Zealand", "Poland", "Spain"),
  ppp     = c(434.95, 3.28, 12.47, 0.68,
              0.72,   3.44,  0.59,  0.47, 0.47,
              1.43,   1.80,  0.56)
)
# Note: for Eurozone countries (France, Ireland, Italy, Latvia,
# Lithuania, Spain) and Croatia (fixed-rate EUR since 2023) the
# PPP is expressed in EUR per USD. For non-Eurozone countries
# (Chile, Czechia, Israel, New Zealand, Poland) it is in national
# currency per USD.

df <- df %>%
  left_join(ppp_2022 %>% select(CNTRYID, ppp), by = "CNTRYID") %>%
  mutate(
    EARNHRC2_PPP = EARNHRC2 / ppp,   # hourly wage in 2022 USD PPP
    # left_join silently coerces CNTRYID from factor to character
    # when joining against ppp_2022 (character CNTRYID) — re-factor
    # explicitly so downstream models treat it as a proper factor
    # (country fixed effects), not as a string.
    CNTRYID = factor(CNTRYID)
  )

# Sanity check: every retained country should have a non-missing ppp
stopifnot(all(!is.na(df$ppp)))


# ---- D. SAMPLE RESTRICTION + DEPENDENT VARIABLE ----
# Restrict to respondents with positive PPP-adjusted hourly
# earnings, required for the log transformation.

df <- df %>%
  filter(!is.na(EARNHRC2_PPP), EARNHRC2_PPP > 0) %>%
  mutate(LOG_WAGE = log(EARNHRC2_PPP))

cat("Observations after wage restriction:", nrow(df), "\n\n")

# ---- D2. TRIM EXTREME WAGE OUTLIERS ----
# The PPP-adjusted wage distribution contains implausible extremes
# (hourly wages below $0.01 or above $10,000). The sample is trimmed
# at the pooled 1st-99th percentile of LOG_WAGE; a single threshold
# is used because wages are on a common scale after PPP conversion.

log_wage_p01 <- quantile(df$LOG_WAGE, 0.01, na.rm = TRUE)
log_wage_p99 <- quantile(df$LOG_WAGE, 0.99, na.rm = TRUE)

cat("LOG_WAGE trim thresholds: [", round(log_wage_p01, 3), ",",
    round(log_wage_p99, 3), "] — implied hourly wage (USD PPP): [",
    round(exp(log_wage_p01), 2), ",", round(exp(log_wage_p99), 2), "]\n")

n_before_trim <- nrow(df)

df <- df %>%
  filter(LOG_WAGE >= log_wage_p01, LOG_WAGE <= log_wage_p99)

cat("Observations removed by wage trimming:", n_before_trim - nrow(df),
    "(", round((n_before_trim - nrow(df)) / n_before_trim * 100, 2), "% )\n")
cat("Observations after wage trimming:", nrow(df), "\n\n")

# Check trimming didn't disproportionately hit any single country
cat("Sample size by country after trimming:\n")
df %>% count(CNTRYID) %>% arrange(desc(n)) %>% print(n = 15)
cat("\n")


# ---- E. OUTLIER HANDLING ----
# Cap weekly working hours at 60 (approx. 95th percentile —
# plausible maximum; values above are treated as data-entry
# noise rather than genuine variation).

df <- df %>%
  mutate(D2_Q11 = ifelse(D2_Q11 > 60, NA, D2_Q11))


# ---- F. MISSINGNESS DIAGNOSTICS (POST-CLEANING) ----

missing_pct <- df %>%
  summarise(across(everything(), ~round(mean(is.na(.)) * 100, 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

cat("=== MISSINGNESS BY VARIABLE (retained countries only) ===\n")
print(missing_pct, n = Inf)


# ---- G. FINAL MODELLING DATASET (LISTWISE DELETION) ----

df_model <- df %>%
  select(-EARNHRC2, -ppp) %>%
  drop_na()

# Remove factor levels with no observations in the final sample
# (e.g. an empty "" level of ISIC1C left over from excluded rows)
df_model <- droplevels(df_model)

# Check: PV2-PV10 and SPFWT0 must not remove any observation: the
# sample must be the same as with the model variables alone
# (N = 26,336).
n_model_vars_only <- df %>%
  select(-EARNHRC2, -ppp, -all_of(setdiff(pv_all, cognitive)), -all_of(weight_var)) %>%
  drop_na() %>%
  nrow()
stopifnot(nrow(df_model) == n_model_vars_only)

# Check: every column of df_model belongs to a block in
# 00_variable_blocks.R (predictors, outcome or auxiliary columns)
check_blocks(df_model)

cat("\nObservations after listwise deletion:", nrow(df_model), "\n")
cat("Variables in modelling dataset:", ncol(df_model), "\n\n")
cat("Variable list:", paste(names(df_model), collapse = ", "), "\n\n")

cat("Sample size by country (final modelling dataset):\n")
df_model %>% count(CNTRYID) %>% arrange(desc(n)) %>% print(n = 15)


# ---- H. TRAIN/TEST SPLIT (80/20, stratified by country) ----

# The split depends only on the seed and on the CNTRYID column, so
# adding PV2-PV10 and SPFWT0 leaves train and test unchanged.
set.seed(42)
train_index <- createDataPartition(df_model$CNTRYID, p = 0.8, list = FALSE)

train <- df_model[ train_index, ]
test  <- df_model[-train_index, ]

cat("\nTraining set:", nrow(train), "observations\n")
cat("Test set:    ", nrow(test),  "observations\n\n")

# Verify country balance across the split
bind_rows(
  train %>% count(CNTRYID) %>% mutate(set = "train"),
  test  %>% count(CNTRYID) %>% mutate(set = "test")
) %>%
  pivot_wider(names_from = set, values_from = n) %>%
  mutate(pct_train = round(train / (train + test) * 100, 1)) %>%
  print(n = 15)


# ---- I. SAVE THE CLEANED OBJECTS FOR THE OTHER SCRIPTS ----

saveRDS(df,       "data/processed/df_clean.rds")   # full cleaned data, pre-listwise-deletion (for descriptives)
saveRDS(df_model, "data/processed/df_model.rds")   # final modelling dataset (post listwise deletion)
saveRDS(train,    "data/processed/train.rds")
saveRDS(test,     "data/processed/test.rds")

cat("\nSaved: df_clean.rds, df_model.rds, train.rds, test.rds\n")


# END OF SCRIPT 01 ----
