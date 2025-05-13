# plantTracker - Anderson 2016 Montana - Bouteloua dactyloides

# Author: Diana Spurite
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : diana.spurite@posteo.de
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.13

# Publication: https://doi.org/10.1890/11-0193.1

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
v_author_year      <- c('anderson_2016')
# Define region abbreviation
v_region_abb       <- c('mt')
# Define growth form (grass - including c4, c3 and short grasses-, forb, shrub)
v_gr_form          <- c('grass')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c(',')


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- sp_list %>% .[c(2),]  


# Modifications to the data structure ------------------------------------------
# Specific plots to exclude, a list of plots
v_mod_plot <- c()

source('pipeline/plant_tracker_02.R')

problem_quadrats <- unique(data$quadrat[!(data$year %in% inv$year)])
print(problem_quadrats)

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
