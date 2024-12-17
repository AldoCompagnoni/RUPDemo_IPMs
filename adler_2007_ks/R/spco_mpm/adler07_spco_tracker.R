# plantTracker - Adler 2007 Kansas - Sphaeralcea coccinea

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.12.11

# Publication: https://doi.org/10.1890/0012-9658(2007)88[2673:LMQFKP]2.0.CO;2

rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'adler_2007'
# Define region abbreviation
region_abb <- 'ks'
# Define growth form (grass, forb, shrub, c4)
gr_form    <- 'forb'
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c(',')

source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- sp_list %>% .[c(5),]  

source('pipeline/plant_tracker_02.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
buff
