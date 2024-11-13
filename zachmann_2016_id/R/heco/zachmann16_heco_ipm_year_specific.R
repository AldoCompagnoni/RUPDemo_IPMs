# IPM year specific - chu 2013 colorado - Hesperostipa comata

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.11.11


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker and IPM mean if the data does not exist
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the data in their respective folder
# 3. check the global environment for all the other objects

# 3.1 df     : view(df)
# 3.2 surv_df: view(surv_df)
# 3.3 grow_df: view(grow_df)
# 3.4 recr_df: view(recr_df)
# pop_counts : view(pop_counts)
# all_pars: view(all_pars)

# 3.5 growth models    : list(gr_mod_yr, gr_mod_yr_2, gr_mod_yr_3)
# 3.6 gorwth var. model: gr_var
# 3.7 survival models  : list(su_mod_yr, su_mod_yr_2, su_mod_yr_3)
# 3.8 recruitment model: rec_mod

# lambda mean: lam_mean_kern


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
rm(list = ls()) 


# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'zachmann_2016'
# Define region abbreviation
region_abb <- 'id'
# Define species 
species <- 'Hesperostipa comata'


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
years_re <- c()

# Models:
#  Going down in complexity of the survival and/or growth model.
# Survival model, 0 means keep the complexity (takes: 0-2)
su_complex_down_by <- c(0)
# Growth model, 0 means keep the complexity (takes: 0-2)
gr_complex_down_by <- c(0)


# Main code --------------------------------------------------------------------
# Run the IPM year specific wrapper function
source('helper_functions/ipm_year_specific.R')
