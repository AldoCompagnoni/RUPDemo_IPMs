# plantTracker - Anderson 2016 Montana - Allium textile


# Author: Diana Spurite
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : diana.spurite@posteo.de
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.03.18

# Publication: https://doi.org/10.1890/11-0193.1

# rm(list = ls())


# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'Anderson_2016'
# Define region abbreviation
region_abb <- 'mt'
# Define growth form (grass - including c4, c3 and short grasses-, forb, shrub)
gr_form    <- 'grass'
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c(',')

source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- sp_list %>% .[c(5),]  


# Modifications to the data structure ------------------------------------------
# Specific plots to exclude, a list of plots
mod_plot <- c()

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

