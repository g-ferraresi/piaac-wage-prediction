# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 08_rq2_sensitivity.R
# RQ2 — two additions to 07_multilevel.R (Section 5.3)
#
#   PART A — Overall return to cognitive skills by country:
#            wage change for a worker who is 1 SD higher in ALL
#            three skills at once (sum of the three country slopes).
#            With highly correlated skills this is more stable than
#            the three separate returns, and is the number used in
#            the text of Chapter 5; the skill-by-skill returns go to
#            the appendix.
#   PART B — Leave-one-country-out: ML1 re-estimated 12 times,
#            each time without one country. Shows whether the slope
#            variances (the answer to RQ2) depend on a single country,
#            Chile in particular (numeracy return +15.4% in ML1).
#
# Same data, standardization and model as 07 (ML1 = random intercept
# + uncorrelated random slopes on the three skills, REML).
# Runtime: 13 model fits, a few minutes.


# Run from the repository root (see README).
library(tidyverse)
library(lme4)

source("R/00_variable_blocks.R")

train <- readRDS("data/processed/train.rds")
check_blocks(train)

country_labels <- c(
  "152" = "Chile", "191" = "Croatia", "203" = "Czechia",
  "250" = "France", "372" = "Ireland", "376" = "Israel",
  "380" = "Italy", "428" = "Latvia", "440" = "Lithuania",
  "554" = "New Zealand", "616" = "Poland", "724" = "Spain"
)

# Standardize skills with the training mean and SD (as in 06-07)
skill_mean <- setNames(map_dbl(cognitive, ~ mean(train[[.x]])), cognitive)
skill_sd   <- setNames(map_dbl(cognitive, ~ sd(train[[.x]])),   cognitive)
train_s <- train
for (v in cognitive) train_s[[v]] <- (train_s[[v]] - skill_mean[v]) / skill_sd[v]

# ML1 formula (as in 07)
fixed_vars <- c(setdiff(controls, "CNTRYID"), cognitive, bigfive)
f_ml1 <- as.formula(paste(
  "LOG_WAGE ~", paste(fixed_vars, collapse = " + "),
  "+ (1 + ", paste(cognitive, collapse = " + "), " || CNTRYID)"
))
f_ols_s4 <- make_formula(c(controls, cognitive, bigfive, paste0("CNTRYID:", cognitive)))

fit_ml1 <- function(data) {
  lmer(f_ml1, data = data, REML = TRUE,
       control = lmerControl(optimizer = "bobyqa"))
}

ml1    <- fit_ml1(train_s)
ols_s4 <- lm(f_ols_s4, data = train_s)
cat("ML1 refitted. Singular:", isSingular(ml1), "\n\n")


# ---- A. OVERALL RETURN TO COGNITIVE SKILLS BY COUNTRY ----

# ML1: country slopes = common slope + predicted country deviation
ml_overall <- coef(ml1)$CNTRYID[, cognitive] %>%
  as_tibble(rownames = "CNTRYID") %>%
  mutate(ml_sum = PVLIT1 + PVNUM1 + PVAPS1) %>%
  select(CNTRYID, ml_sum)

# OLS S4: main effect + interaction (0 for the reference country)
b <- coef(ols_s4)
ref_country <- levels(train_s$CNTRYID)[1]
ols_overall <- tibble(CNTRYID = ml_overall$CNTRYID) %>%
  mutate(ols_sum = map_dbl(CNTRYID, function(cn) {
    sum(map_dbl(cognitive, function(sk) {
      b[[sk]] + if (cn == ref_country) 0 else b[[paste0("CNTRYID", cn, ":", sk)]]
    }))
  }))

common_sum <- sum(fixef(ml1)[cognitive])

overall <- ml_overall %>%
  left_join(ols_overall, by = "CNTRYID") %>%
  mutate(country    = unname(country_labels[CNTRYID]),
         pct_ml     = 100 * (exp(ml_sum) - 1),
         pct_ols_s4 = 100 * (exp(ols_sum) - 1)) %>%
  arrange(desc(pct_ml))

cat("=== A. OVERALL RETURN: % WAGE CHANGE FOR +1 SD IN ALL THREE SKILLS ===\n")
cat("Common return (ML1 fixed part):", sprintf("%+.1f%%", 100 * (exp(common_sum) - 1)), "\n\n")
overall %>%
  transmute(country,
            ML1    = sprintf("%+.1f%%", pct_ml),
            OLS_S4 = sprintf("%+.1f%%", pct_ols_s4)) %>%
  print(n = Inf)

write_csv(overall, "output/tables/table_3x_rq2_overall_returns.csv")

# Figure: overall return to cognitive skills by country (ML1)
p_returns <- overall %>%
  mutate(country = fct_reorder(country, pct_ml)) %>%
  ggplot(aes(x = pct_ml, y = country)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct_ml)),
            hjust = -0.15, size = 3.2, colour = "grey30") +
  geom_vline(xintercept = 0, colour = "grey40") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Returns to cognitive skills vary widely across countries",
       subtitle = "Multilevel model (ML1) with country-specific slopes",
       x = "% wage change for +1 SD in all three cognitive skills",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/fig_3x_rq2_returns_by_country.png", p_returns,
       width = 7, height = 5, dpi = 300)


# ---- B. LEAVE-ONE-COUNTRY-OUT ----

loo <- map_dfr(names(country_labels), function(cn) {
  d <- train_s %>% filter(CNTRYID != cn) %>% droplevels()
  m <- fit_ml1(d)
  vc <- as.data.frame(VarCorr(m))
  tibble(
    excluded  = country_labels[[cn]],
    singular  = isSingular(m),
    sd_int    = vc$sdcor[vc$var1 == "(Intercept)" & !is.na(vc$var1)],
    sd_lit    = vc$sdcor[vc$var1 == "PVLIT1" & !is.na(vc$var1)],
    sd_num    = vc$sdcor[vc$var1 == "PVNUM1" & !is.na(vc$var1)],
    sd_aps    = vc$sdcor[vc$var1 == "PVAPS1" & !is.na(vc$var1)],
    pct_lit   = 100 * (exp(fixef(m)[["PVLIT1"]]) - 1),
    pct_num   = 100 * (exp(fixef(m)[["PVNUM1"]]) - 1),
    pct_aps   = 100 * (exp(fixef(m)[["PVAPS1"]]) - 1)
  )
})

vc_all <- as.data.frame(VarCorr(ml1))
full_row <- tibble(
  excluded = "(none: full ML1)", singular = isSingular(ml1),
  sd_int = vc_all$sdcor[vc_all$var1 == "(Intercept)" & !is.na(vc_all$var1)],
  sd_lit = vc_all$sdcor[vc_all$var1 == "PVLIT1" & !is.na(vc_all$var1)],
  sd_num = vc_all$sdcor[vc_all$var1 == "PVNUM1" & !is.na(vc_all$var1)],
  sd_aps = vc_all$sdcor[vc_all$var1 == "PVAPS1" & !is.na(vc_all$var1)],
  pct_lit = 100 * (exp(fixef(ml1)[["PVLIT1"]]) - 1),
  pct_num = 100 * (exp(fixef(ml1)[["PVNUM1"]]) - 1),
  pct_aps = 100 * (exp(fixef(ml1)[["PVAPS1"]]) - 1)
)

loo_table <- bind_rows(full_row, loo)

cat("\n=== B. LEAVE-ONE-COUNTRY-OUT: ML1 WITHOUT ONE COUNTRY AT A TIME ===\n")
cat("sd_*: SD across countries of the random effect (log points per 1 SD of skill)\n")
cat("pct_*: common return, % wage change per 1 SD\n\n")
loo_table %>%
  mutate(across(starts_with("sd_"), ~ round(.x, 4)),
         across(starts_with("pct_"), ~ round(.x, 2))) %>%
  print(n = Inf, width = Inf)

write_csv(loo_table, "output/tables/table_3x_rq2_leave_one_out.csv")


# END OF SCRIPT 08 ----
