# MPM year specific - Anderson 2016 Montana - Allium textile

# Author: Diana Spurite
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : diana.spurite@posteo.de
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.03.14

# Publication : https://doi.org/10.1890/11-0193.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder


# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'anderson_2016'
# Define region abbreviation
region_abb  <- 'mt'
# Define the species variable
species     <- 'Allium textile'
# Type of population model
mod_type    <- 'mpm'


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
years_re <- c()
# Age:
#  Set a age threshold at which all ages get the same survival probability
#  If unspecified it's 'c(1)'
surv_age_threshold <- c(1)


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

