# plantTracker - Christensen 2021 New Mexico - Sporobolus airoides

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.12.11

# Publication: https://doi.org/10.1002/ecy.3530

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year      <- c('christensen_2021')
# Define region abbreviation
v_region_abb       <- c('nm')
# Define growth form (grass, forb, shrub, c4)
v_gr_form          <- c('grass')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c(',')


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list)
target_spec <- sp_list %>% .[c(9),]  

# Removing problematic quadrats
mod_plot <- c('L5')

source('pipeline/plant_tracker_02.R')


# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
v_buff
