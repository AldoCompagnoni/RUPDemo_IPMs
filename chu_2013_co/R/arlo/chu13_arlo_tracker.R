# plantTracker - Chu 2013 Colorado - Aristida longiseta

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.12.11

# Publication: https://doi.org/10.1890/13-0121.1

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year      <- c('chu_2013')
# Define region abbreviation
v_region_abb       <- c('co')
# Define growth form (grass, forb, shrub, c4)
v_gr_form          <- c('grass')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c('\t')


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 10)
target_spec <- sp_list %>% .[c(7),]  

source('pipeline/plant_tracker_02.R')


# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
v_buff
