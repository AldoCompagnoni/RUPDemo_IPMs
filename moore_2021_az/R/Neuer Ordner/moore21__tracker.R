# plantTracker - Moore 2021 Arizona -

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.01.27

# Publication: 

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'moore_21'
# Define region abbreviation
region_abb <- 'az'
# Define growth form (grass, forb, shrub, c4)
gr_form    <- 'grass'
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c(',')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- as.data.frame(c('Muhlenbergia minutissima'))

source('pipeline/plant_tracker_02.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
buff
