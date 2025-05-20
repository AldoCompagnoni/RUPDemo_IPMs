# MPM mean - Anderson 2016 Arizona - Ambrosia artemisiifolia

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.20

# Publication: https://doi.org/10.1890/0012-9658(2007)88[2673:LMQFKP]2.0.CO;2


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder


# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year <- c('anderson_2016')
# Define region abbreviation
v_region_abb  <- c('az')
# Define the species variable
v_species     <- c('Ambrosia artemisiifolia')
# Type of population model
v_mod_type    <- c('mpm')


# CHECK -- Adaptions to the models ---------------------------------------------
# Years:
#  Removal of certain years if unspecified nothing is removed
v_years_re  <- c()
# Age:
#  Set a max age of age classes
#  Undefined it is max age is 1, therefore we have age classes 0 and 1 only
v_age_class <- c()


# MPM mean pipeline ------------------------------------------------------------
source('pipeline/mpm_mean.R')


# Output -----------------------------------------------------------------------
# Models:
# Survival:
summary(mod_su)

# Recruit1: recruitment being there or not
summary(mod_re_yn)

# Recruit2: per-capita recruitment (pcr) *conditional* on recruitment happening
summary(mod_re_pc_c)

# Recruit3: per capita recruitment *unconditinoal*
summary(mod_re_pc_uc)


# Plots: 
fig_age_hist
fig_age_surv
fig_su
fig_re


# Data frames:
# Population count: mean lambda and survival and recruitment predictions
#  alongside observed population count and growth rate
print(pop_counts)