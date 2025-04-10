# plantTracker - Zachmann 2016 Idaho - Arnica fulgens

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# Publication: https://doi.org/10.1890/10-0404.1

# rm(list = ls())

# Specifications ---------------------------------------------------------------
# Define publication 
v_author_year      <- c('zachmann_2016')
# Define region abbreviation
v_region_abb       <- c('id')
# Define growth form (grass, forb, shrub, c4)
v_gr_form          <- c('forb')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list)
target_spec <- sp_list %>% .[c(3),]  

source('pipeline/plant_tracker_02.R')


# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
v_buff
