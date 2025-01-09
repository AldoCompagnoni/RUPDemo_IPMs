# IPM mean - Chu 2013 Colorado - Sitanion hystrix

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.12.03

# Publication: https://doi.org/10.1890/13-0121.1

# Read in and clean the data
#  explore the overall-years rates
#  set up the vital rate data-frames for the year specific 
#  build the ipm from scratch
#  build the ipm with `ipmr`


# Comments ---------------------------------------------------------------------
# 0. !!! Please define all key variables in the in the corresponding section !!!
# 1. The pipeline runs plant tracker only if the data does not exist already
# 2. Find all the graphics in the result folder of the respective species
#     and the growth, survival and recruitment data in the data folder


# Clean up ---------------------------------------------------------------------
rm(list = ls())


# Key variables ----------------------------------------------------------------
# Define publication 
author_year <- 'chu_2013'
# Define region abbreviation
region_abb <- 'co'
# Define species 
species <- 'Sitanion hystrix'


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
years_re <- c()

# Models:
#  Changing to the next best complexity of the survival and/or growth model.
# Survival model, 0 means keep the complexity (takes: 0-2)
su_complex <- c(0)
# Growth model, 0 means keep the complexity (takes: 0-2)
gr_complex <- c(0)


# Main pipeline ----------------------------------------------------------------
# Run the ipm mean wraper function
source('pipeline/ipm_mean.R')


# Data 1 -----------------------------------------------------------------------
# Raw
skim(df)

# Survival
skim(surv_df)

# Growth
skim(grow_df)

# Recruitment
skim(recr_df)


# Models -----------------------------------------------------------------------
# Survival 
list(su_mod_mean, su_mod_mean_2, su_mod_mean_3)

# Growth
list(gr_mod_mean, gr_mod_mean_2, gr_mod_mean_3)

# Growth variation
gr_var_m

# Recruitment
rec_mod_mean


# Building the IPM from scratch ------------------------------------------------
# Parameters
tibble(parameter = names(pars), value = unlist(pars))

# Mean population growth rate
lam_mean

# Observed population growth rate
skim(pop_counts) 

# Geometric mean of yearly population growth rates
lam_mean_count 

# Overall (aggregated) population growth rate
lam_mean_overall 


# Building the IPM with ipmr ---------------------------------------------------
ipmr_p 
plot(ipmr_p)
