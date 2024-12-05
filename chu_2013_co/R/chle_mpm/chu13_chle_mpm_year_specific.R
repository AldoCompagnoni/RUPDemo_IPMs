# MPM year specific - chu 2013 colorado - Chenopodium leptophyllum

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# Publication : https://doi.org/10.1890/13-0121.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder
# 3. check the global environment for all the other objects
# 3.1 recr_pred_df: view(recr_pred_df)
#      year specific recruitment predictions
# 3.2 pop_counts:   view(pop_counts)
#      year specific lambdas and underlying survival and recruitment predictions
#      alongside observed population count and growth rate
# 3.3 model survival: summary(mod_surv)
# 3.4 model recruit1: summary(yesno_mod)
#      recruitment being there or not
# 3.4 model recruit2: summary(pcr_mod)
#      per-capita recruitment (pcr) *conditional* on recruitment happening
# 3.5 model recruit3: summary(pcr_uc_mod)
#      per capita recruitment *unconditinoal*


rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'chu_2013'
# Define region abbreviation
region_abb  <- 'co'
# Define the species variable
species     <- 'Chenopodium leptophyllum'
# Type of population model
mod_type    <- 'mpm'


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
years_re <- c()


# MPM pipeline -- with model comparison -- -------------------------------------
source('helper_functions/mpm_year_specific.R')


# Published vs our lambdas -----------------------------------------------------
rmse(compute_df$obs_pgr,
     compute_df$pub_lam)
rmse(compute_df$obs_pgr,
     compute_df$lam)
rmse(compute_df$obs_pgr,
     compute_df$pub_proj_lam)
rmse(compute_df$obs_pgr,
     compute_df$proj_lam)

