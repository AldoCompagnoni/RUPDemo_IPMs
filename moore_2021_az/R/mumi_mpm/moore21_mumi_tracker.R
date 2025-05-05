# plantTracker - Moore 2021 Arizona - Muhlenbergia minutissima

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.04.30

# Publication: https://doi.org/10.1002/ecy.3661

# rm(list = ls())

# Specifications ---------------------------------------------------------------
# Define publication 
v_author_year      <- c('moore_2021')
# Define region abbreviation
v_region_abb       <- c('az')
# Define growth form (density, cover)
v_gr_form          <- c('Density')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c(',')


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- sp_list %>% .[c(1),]  

source('pipeline/plant_tracker_02.R')


# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
v_buff