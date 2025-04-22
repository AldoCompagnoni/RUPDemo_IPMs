# plantTracker - Moore 2021 Arizona - Antennaria parvifolia

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.04.17

# Publication: https://doi.org/10.1002/ecy.3661

# rm(list = ls())

# Data -------------------------------------------------------------------------
# Define publication 
author_year      <- c('moore_2021')
# Define region abbreviation
region_abb       <- c('az')
# Define growth form (density, cover)
m_type           <- c('Cover')
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c(',')

source('pipeline/plant_tracker_01_moore.R')

# Select the x_th species (target species)
head(sp_list, 20)
target_spec <- sp_list[13,]

source('pipeline/plant_tracker_02_moore.R')

# Exploration ------------------------------------------------------------------
# Quadrat inventory
quad_inv
# Polygon data
dat_target_spec
# Buffer size - regular and genet
st_bbox(dat_target_spec)
buff
