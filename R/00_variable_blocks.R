# ==============================================================================
# Predicting wages across OECD countries: a multilevel analysis using 
# cross-national adult skills data
# Script: 00_variable_blocks.R
#
# Single definition of the variable blocks used by every script.
# It is not run on its own: the other scripts load it with
#   source("R/00_variable_blocks.R")
#
# Every model formula is built explicitly from these blocks with
# make_formula(), never with "LOG_WAGE ~ .", because the data also
# contain auxiliary columns (plausible values 2-10 and the sampling
# weight) that must not be used as predictors.


# ---- A. PREDICTOR BLOCKS ----

# Personal and family controls (country dummies included)
personal_controls <- c("CNTRYID", "AGEG10LFS", "GENDER_R", "EDCAT8_TC1",
                       "NFE12C2", "C2_Q10", "NATIVELANG", "IMGENC2",
                       "PAREDC2", "PARTNER")

# Job characteristics (channels through which skills affect wages)
job_controls <- c("YRSGETC2", "ISIC1C", "ISCOSKIL4", "D2_Q07a", "D2_Q11")

controls <- c(personal_controls, job_controls)

# Cognitive skills — plausible value 1 (main estimates, test metrics)
cognitive <- c("PVLIT1", "PVNUM1", "PVAPS1")

# Big Five — within-country z-scores as provided by the OECD in the
# Public Use File (standardisation done by the OECD, not in this
# pipeline). 3-item domain scores common to BFI-2-S and BFI-2-XS.
bigfive <- c("AGRE", "CONS", "EMOS", "EXTR", "OPEM")


# ---- B. PLAUSIBLE VALUES AND WEIGHTS ----
# Names verified in piaac-cy2-international-codebook.xlsx, sheet PUF.

pv_lit <- paste0("PVLIT", 1:10)
pv_num <- paste0("PVNUM", 1:10)
pv_aps <- paste0("PVAPS", 1:10)
pv_all <- c(pv_lit, pv_num, pv_aps)

# The three cognitive-skill variables for plausible value m (1-10).
# Used to re-estimate a model once per PV (Rubin's rules).
cognitive_pv <- function(m) c(paste0("PVLIT", m), paste0("PVNUM", m), paste0("PVAPS", m))

weight_var <- "SPFWT0"   # "Final full sample weight"


# ---- C. AUXILIARY COLUMNS ----
# Kept in df_model/train/test but NEVER used as predictors:
#   - EARNHRC2_PPP: wage in USD PPP (LOG_WAGE = log(EARNHRC2_PPP))
#   - plausible values 2-10: used only one PV at a time, in place of PV1
#   - SPFWT0: used only for the weighted robustness check

aux_vars <- c("EARNHRC2_PPP", setdiff(pv_all, cognitive), weight_var)


# ---- D. HELPERS ----

# Model formula from a vector of predictor names
make_formula <- function(vars, response = "LOG_WAGE") {
  reformulate(vars, response = response)
}

# Stops the script if the data contain a column that is not assigned
# to any block (or a block variable is missing). If a variable is
# added to or removed from 01_data_cleaning.R, update this file too.
check_blocks <- function(data) {
  expected <- c(controls, cognitive, bigfive, "LOG_WAGE", aux_vars)
  missing_cols <- setdiff(expected, names(data))
  extra_cols   <- setdiff(names(data), expected)
  if (length(missing_cols) > 0 || length(extra_cols) > 0) {
    stop("Variable blocks do not match the data.\n",
         "  Missing from data: ", paste(missing_cols, collapse = ", "), "\n",
         "  Not assigned to any block: ", paste(extra_cols, collapse = ", "))
  }
  invisible(TRUE)
}


# END OF SCRIPT 00 ----
