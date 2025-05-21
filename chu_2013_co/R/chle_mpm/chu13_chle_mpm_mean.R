# MPM mean - Chu 2013 Colorado - Chenopodium leptophyllum

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.13

# Publication: https://doi.org/10.1890/13-0121.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder


# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year <- c('chu_2013')
# Define region abbreviation
v_region_abb  <- c('co')
# Define the species variable
v_species     <- c('Chenopodium leptophyllum')
# Type of population model
v_mod_type    <- c('mpm')
