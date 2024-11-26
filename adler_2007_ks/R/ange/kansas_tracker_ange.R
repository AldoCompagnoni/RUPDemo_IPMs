# plantTracker - adler 2007 kansas - Andropogon gerardii

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# Publication: https://doi.org/10.1890/10-0404.1

rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'adler_2007'
# Define region abbreviation
region_abb <- 'ks'
# Define growth form (grass, forb, shrub, c4)
gr_form    <- 'grass'
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c()

source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 10)
target_spec <- sp_list %>% .[c(2),]  

source('pipeline/plant_tracker_02.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
