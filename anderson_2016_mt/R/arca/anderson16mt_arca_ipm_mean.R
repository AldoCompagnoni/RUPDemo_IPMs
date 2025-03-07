# IPM mean - Anderson 2016 Montana - Artemisia cana

# Author: Diana Spurite
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : diana.spurite@posteo.de
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.03.07

# Publication: https://doi.org/10.1890/11-0193.1

# Read in and clean the data, explore the overall-years rates,
#  set up the vital rate data-frames for the year specific,
#  build the ipm from scratch, build the ipm with `ipmr`


# Clean up ---------------------------------------------------------------------
# rm(list = ls())


# Key variables ----------------------------------------------------------------
# Define publication 
v_author_year <- c('anderson_2016')
# Define region abbreviation
v_region_abb  <- c('mt')
# Define species 
v_species     <- c('Artemisia cana')


# CHECK -- Adaptions -----------------------------------------------------------
# Removal of certain years if unspecified nothing is removed
v_years_re       <- c(33-38)
# Define size threshold
v_size_threshold <- c(-10.47)
# Set a complexity to the growth and survival model 
# (NULL = highest AIC / 0 = intercept / 1 = linear / 2 = quadratic / 3 = cubic)
v_mod_set_gr     <- c()
v_mod_set_su     <- c()


# Main pipeline ----------------------------------------------------------------
# Run the ipm mean wraper function
source('pipeline/ipm_mean.R')

traceback()
summary(x)
summary(y)
print(any(is.na(x)))  # TRUE if x has NA
print(any(is.na(y)))  # TRUE if y has NA
print(any(is.infinite(x)))  # TRUE if x has Inf
print(any(is.infinite(y)))  # TRUE if y has Inf
print(range(y, na.rm = TRUE))  # Check for extreme values

data_clean <- na.omit(data.frame(x, y))

summary(recr_nona_nr_quad$nr_quad)
print(sum(is.na(recr_nona_nr_quad$nr_quad)))  # Count NA values
print(sum(is.infinite(recr_nona_nr_quad$nr_quad)))  # Count Inf values

print(any(recr_nona_nr_quad$nr_quad <= 0))
recr_nona_nr_quad$nr_quad <- recr_nona_nr_quad$nr_quad + 0.01

summary(recr_nona_nr_quad$nr_quad)
table(recr_nona_nr_quad$nr_quad)


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

