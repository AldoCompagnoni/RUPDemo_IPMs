# plantTracker - Anderson 2016 Montana - Sphaeralcea coccinea

# Author: 
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.12.03

# Publication : 

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'anderson_2016'
# Define region abbreviation
region_abb <- 'mt'
# Define growth form (grass, forb, shrub, c4)
gr_form    <- 'forb'
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c(',')

source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list)
target_spec <- sp_list %>% .[c(1),]  

source('pipeline/plant_tracker_02.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
buff