# plantTracker for chu 2013 colorado Pseudoroegneria spicata

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# 

rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'zachmann_2016'
# Define region abbreviation
region_abb <- 'id'
# Define growth form (grass, forb, shrub)
gr_form    <- 'grass'

source('helper_functions/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list)
target_spec <- sp_list %>% .[c(1),]  

source('helper_functions/plant_tracker_02.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat
# Buffer size - regular and genet
st_bbox(dat_target_spec)
buff
