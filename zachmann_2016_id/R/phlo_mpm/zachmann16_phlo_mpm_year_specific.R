# MPM year specific - Zachmann 2016 Idaho - Phlox longifolia

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.13

# Publication : https://doi.org/10.1890/10-0404.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder


# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year <- c('zachmann_2016')
# Define region abbreviation
v_region_abb  <- c('id')
# Define the species variable
v_species     <- c('Phlox longifolia')
# Type of population model
v_mod_type    <- c('mpm')


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
v_years_re <- c()
# Age:
#  Set a age threshold at which all ages get the same survival probability
#  If unspecified it's 'c(1)'
surv_age_threshold <- c(4)


# MPM pipeline -- with model comparison -- -------------------------------------
source('pipeline/mpm_year_specific.R')


# Output -----------------------------------------------------------------------
# Data frames:
# Recruitment: year specific recruitment predictions
skim(recr_pred_df)

# Population count: year specific lambdas and underlying survival 
#  and recruitment predictions 
#  alongside observed population count and growth rate   
skim(pop_counts)

# Models:
# Survival: 
summary(mod_surv)

# Recruit1: recruitment being there or not 
summary(yesno_mod)

# Recruit2: per-capita recruitment (pcr) *conditional* on recruitment happening
summary(pcr_mod)

# Recruit3: per capita recruitment *unconditinoal*
summary(pcr_uc_mod)
