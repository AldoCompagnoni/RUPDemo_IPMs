# IPM mean - chu 2013 colorado - Hesperostipa comata

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.11.08

# reading in, and cleaning the data
#  exploring the overall-years rates
#  setting up the vital rate data-frames for the year specific 

# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the growth, survival and recruitment data in their respective folder
# 3. check the global environment for all the other objects

# 3.1 df     : view(df)
# 3.2 surv_df: view(surv_df)
# 3.3 grow_df: view(grow_df)
# 3.4 recr_df: view(recr_df)

# 3.5 growth models    : list(gr_mod_mean, gr_mod_mean_2, gr_mod_mean_3)
# 3.6 gorwth var. model: gr_var_m
# 3.7 survival models  : list(su_mod_mean, su_mod_mean_2, su_mod_mean_3)
# 3.8 recruitment model: rec_mod_mean

rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'zachmann_2016'
# Define region abbreviation
region_abb <- 'id'
# Define species 
species <- 'Hesperostipa comata'

# Run the ipm mean wraper function
source('helper_functions/ipm_mean.R')

# Provide a summary of the cleaned data
skim(df)

# Output -----------------------------------------------------------------------
# Reproduction per capita summary
repr_pc_m