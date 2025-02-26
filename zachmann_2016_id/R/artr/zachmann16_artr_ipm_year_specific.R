# IPM year specific - Zachmann 2016 Idaho - Artemisia tripartita

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.02.26

# Publication: https://doi.org/10.1890/10-0404.1


# Comments ---------------------------------------------------------------------
# 0. !!! Please define all key variables in the in the corresponding section !!!
# 1. The pipeline runs plant tracker and IPM mean only if needed 
# 2. Find all the graphics in the result folder of the respective species
#     and the growth, survival and recruitment data in the data folder


# Clean up ---------------------------------------------------------------------
# rm(list = ls())


# Key variables ----------------------------------------------------------------
# Define publication 
v_author_year <- c('zachmann_2016')
# Define region abbreviation
v_region_abb  <- c('id')
# Define species 
v_species     <- c('Artemisia tripartita')


# CHECK -- Adaptions -----------------------------------------------------------
# Removal of certain years if unspecified nothing is removed
v_years_re       <- c()
# Define size threshold
v_size_threshold <- c()
# Set a complexity to the growth and survival model 
# (NULL = highest AIC / 0 = intercept / 1 = linear / 2 = quadratic / 3 = cubic)
v_mod_set_gr     <- c()
v_mod_set_su     <- c()


# Main pipeline ----------------------------------------------------------------
# Run the IPM year specific wrapper function
source('pipeline/ipm_year_specific.R')


# Data -------------------------------------------------------------------------

# Dataframe
skim(df)

# Survival
skim(surv_df)

# Grow_df
skim(grow_df)

# Recruitment 
skim(recr_df)

# Population counts
skim(pop_counts)


# Models -----------------------------------------------------------------------
# Survival 
list(su_mod_yr, su_mod_yr_2, su_mod_yr_3)

# Growth
list(gr_mod_yr, gr_mod_yr_2, gr_mod_yr_3)

# Growth variation
gr_var

# Recruitment
rec_mod


# Building the IPM from scratch ------------------------------------------------
# All parameters
skim(all_pars)

# Mean population growth rate
lam_mean_yr

# Observed population growth rate
skim(pop_counts) 

# Geometric mean of yearly population growth rates
lam_mean_count 

# Overall (aggregated) population growth rate
lam_mean_overall

# Mean lambda 
lam_mean_kern


# Building the IPM with ipmr ---------------------------------------------------
ipmr_yr 
