# IPM mean - Moore 2021 Arizona - Koeleria macrantha

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.02.21

# Publication: https://doi.org/10.1002/ecy.3661

# Read in and clean the data, explore the overall-years rates,
#  set up the vital rate data-frames for the year specific,
#  build the ipm from scratch, build the ipm with `ipmr`


# Clean up ---------------------------------------------------------------------
# rm(list = ls())


# Key variables ----------------------------------------------------------------
# Define publication 
v_author_year <- c('moore_2021')
# Define region abbreviation
v_region_abb  <- c('az')
# Define species 
v_species     <- c('Koeleria macrantha')


# CHECK -- Adaptions -----------------------------------------------------------
# Removal of certain years if unspecified nothing is removed
v_years_re       <- c()
# Define size threshold
v_size_threshold <- c(-12.7)
# Set a complexity to the growth and survival model 
# (NULL = highest AIC / 0 = intercept / 1 = linear / 2 = quadratic / 3 = cubic)
v_mod_set_gr     <- c(2)
v_mod_set_su     <- c()


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
