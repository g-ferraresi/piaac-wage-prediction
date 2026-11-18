# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 02_descriptives.R
#
# Descriptive statistics for Chapter 3 (Data).
# Reads df_model.rds (01_data_cleaning.R): the final analytic
# sample (N = 26,336), the same used by all the estimation scripts.
#   A. Sample composition by country
#   B. LOG_WAGE: overall and by country
#   C. Cognitive skills
#   D. Big Five personality traits

# Run from the repository root (see README).
library(tidyverse)

df_model <- readRDS("data/processed/df_model.rds")

country_labels <- c(
  "152" = "Chile", "191" = "Croatia", "203" = "Czechia",
  "250" = "France", "372" = "Ireland", "376" = "Israel",
  "380" = "Italy", "428" = "Latvia", "440" = "Lithuania",
  "554" = "New Zealand", "616" = "Poland", "724" = "Spain"
)

df_model <- df_model %>%
  mutate(country_name = dplyr::recode(as.character(CNTRYID), !!!country_labels))

df_model <- df_model %>%
  mutate(wage_usd = exp(LOG_WAGE))   # hourly wage in USD PPP

# ---- A. OVERALL SAMPLE SUMMARY ----

cat("=== A. OVERALL SAMPLE SUMMARY ===\n")

n_total <- nrow(df_model)
cat("Total observations (12 countries, final analytic/modelling sample):", n_total, "\n")

sample_summary <- df_model %>%
  group_by(country_name) %>%
  summarise(
    n           = n(),
    pct_sample  = round(n() / n_total * 100, 1)
  ) %>%
  arrange(desc(n))

print(sample_summary, n = 15)

# Sample composition by country
write_csv(sample_summary, "output/tables/table_2x_sample_composition.csv")


# ---- B. DEPENDENT VARIABLE — LOG_WAGE (detailed) ----

cat("\n=== B. LOG_WAGE — OVERALL DESCRIPTIVES ===\n")

y_overall <- df_model %>%
  filter(!is.na(LOG_WAGE)) %>%
  summarise(
    n     = n(),
    mean  = round(mean(LOG_WAGE), 3),
    sd    = round(sd(LOG_WAGE), 3),
    min   = round(min(LOG_WAGE), 3),
    p5    = round(quantile(LOG_WAGE, 0.05), 3),
    p25   = round(quantile(LOG_WAGE, 0.25), 3),
    median= round(median(LOG_WAGE), 3),
    p75   = round(quantile(LOG_WAGE, 0.75), 3),
    p95   = round(quantile(LOG_WAGE, 0.95), 3),
    max   = round(max(LOG_WAGE), 3),
    skewness = round(mean((LOG_WAGE - mean(LOG_WAGE))^3) / sd(LOG_WAGE)^3, 3)
  )
print(y_overall)

cat("\nImplied hourly wage (2022 USD PPP):\n")

wage_usd_summary <- df_model %>%
  filter(!is.na(EARNHRC2_PPP)) %>%
  summarise(
    arith_mean = round(mean(EARNHRC2_PPP), 2),
    geom_mean  = round(exp(mean(log(EARNHRC2_PPP))), 2),
    min        = round(min(EARNHRC2_PPP), 2),
    p5         = round(quantile(EARNHRC2_PPP, 0.05), 2),
    p25        = round(quantile(EARNHRC2_PPP, 0.25), 2),
    median     = round(median(EARNHRC2_PPP), 2),
    p75        = round(quantile(EARNHRC2_PPP, 0.75), 2),
    p95        = round(quantile(EARNHRC2_PPP, 0.95), 2),
    max        = round(max(EARNHRC2_PPP), 2)
  )
print(wage_usd_summary)

# Cross-country heterogeneity (motivates multilevel modelling) 
cat("\n=== LOG_WAGE BY COUNTRY (cross-country heterogeneity — RQ2) ===\n")

y_by_country <- df_model %>%
  filter(!is.na(LOG_WAGE)) %>%
  group_by(country_name) %>%
  summarise(
    n      = n(),
    mean   = round(mean(LOG_WAGE), 3),
    sd     = round(sd(LOG_WAGE), 3),
    median = round(median(LOG_WAGE), 3)
  ) %>%
  arrange(desc(mean))

print(y_by_country, n = 15)

cat("\nRange of country means:",
    round(min(y_by_country$mean), 3), "(", y_by_country$country_name[which.min(y_by_country$mean)], ")",
    "to",
    round(max(y_by_country$mean), 3), "(", y_by_country$country_name[which.max(y_by_country$mean)], ")\n")
cat("Ratio of implied wages (max/min country mean):",
    round(exp(max(y_by_country$mean)) / exp(min(y_by_country$mean)), 2), "\n")

write_csv(y_by_country, "output/tables/table_2x_logwage_by_country.csv")

# Plots 

# Distribution of LOG_WAGE
p_hist <- ggplot(df_model %>% filter(!is.na(LOG_WAGE)), aes(x = LOG_WAGE)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean(df_model$LOG_WAGE, na.rm = TRUE),
             color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Distribution of log hourly wage (2022 USD PPP)",
       subtitle = "Dashed line = sample mean",
       x = "log(hourly wage, PPP-adjusted)", y = "Count") +
  theme_minimal()
print(p_hist)
ggsave("output/figures/fig_2x_logwage_histogram.png", p_hist, width = 7, height = 4.5, dpi = 300)

# Cross-country boxplot
p_box <- df_model %>%
  filter(!is.na(LOG_WAGE)) %>%
  ggplot(aes(x = reorder(country_name, LOG_WAGE, FUN = median), y = LOG_WAGE)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.4) +
  coord_flip() +
  labs(title = "Log hourly wage (PPP) by country",
       subtitle = "Countries ordered by median",
       x = NULL, y = "log(hourly wage, PPP-adjusted)") +
  theme_minimal()
print(p_box)
ggsave("output/figures/fig_2x_logwage_by_country.png", p_box, width = 7, height = 5, dpi = 300)


# ---- C. COGNITIVE SKILLS — compact descriptives ----

cat("\n=== C. COGNITIVE SKILLS — DESCRIPTIVES ===\n")

cognitive_summary <- df_model %>%
  select(PVLIT1, PVNUM1, PVAPS1) %>%
  pivot_longer(everything(), names_to = "skill", values_to = "score") %>%
  filter(!is.na(score)) %>%
  group_by(skill) %>%
  summarise(
    n    = n(),
    mean = round(mean(score), 1),
    sd   = round(sd(score), 1),
    min  = round(min(score), 1),
    max  = round(max(score), 1)
  ) %>%
  mutate(skill = dplyr::recode(skill,
                               PVLIT1 = "Literacy",
                               PVNUM1 = "Numeracy",
                               PVAPS1 = "Adaptive problem solving"))

print(cognitive_summary)
write_csv(cognitive_summary, "output/tables/table_2x_cognitive_skills.csv")

# Brief correlation note (anticipates the VIF/multicollinearity check
# in the OLS assumption checking block)
cog_cor <- df_model %>%
  select(PVLIT1, PVNUM1, PVAPS1) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2)

cat("\nCorrelation among cognitive skill measures:\n")
print(cog_cor)


# ---- D. BIG FIVE PERSONALITY TRAITS — compact descriptives ----

cat("\n=== D. BIG FIVE — DESCRIPTIVES ===\n")

bigfive_summary <- df_model %>%
  select(AGRE, CONS, EMOS, EXTR, OPEM) %>%
  pivot_longer(everything(), names_to = "trait", values_to = "score") %>%
  filter(!is.na(score)) %>%
  group_by(trait) %>%
  summarise(
    n    = n(),
    mean = round(mean(score), 3),
    sd   = round(sd(score), 3),
    min  = round(min(score), 3),
    max  = round(max(score), 3)
  ) %>%
  mutate(trait = dplyr::recode(trait,
                               AGRE = "Agreeableness",
                               CONS = "Conscientiousness",
                               EMOS = "Emotional stability",
                               EXTR = "Extraversion",
                               OPEM = "Openness to experience"))

print(bigfive_summary)
write_csv(bigfive_summary, "output/tables/table_2x_bigfive.csv")

# Check: within-country standardization of the Big Five
# The OECD provides the scores as within-country z-scores computed
# on the whole adult population: in the analytic sample, country
# means close to 0 and SDs close to 1 are expected.

cat("\n=== CHECK: Big Five standardisation — pooled vs. country-specific ===\n")
bigfive_by_country <- df_model %>%
  group_by(country_name) %>%
  summarise(across(c(AGRE, CONS, EMOS, EXTR, OPEM),
                   ~round(mean(., na.rm = TRUE), 3)))
print(bigfive_by_country, n = 15)

# Standard deviation by country
bigfive_sd_by_country <- df_model %>%
  group_by(country_name) %>%
  summarise(across(c(AGRE, CONS, EMOS, EXTR, OPEM),
                   ~round(sd(., na.rm = TRUE), 3)))
cat("\nSD by country (decisive check):\n")
print(bigfive_sd_by_country, n = 15)


write_csv(bigfive_by_country, "output/tables/table_2x_bigfive_mean_by_country.csv")
write_csv(bigfive_sd_by_country, "output/tables/table_2x_bigfive_sd_by_country.csv")


# END OF SCRIPT 02 ----
