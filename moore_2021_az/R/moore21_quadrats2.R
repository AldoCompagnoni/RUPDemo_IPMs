# Quadrat data - Moore 2021 Arizona 

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.01.30

# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs

# Outline
# This script uses the plantTracker package 
# (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data 
# that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) 
# of the perennial grasses.


# Publication: https://doi.org/10.1002/ecy.3530

# rm(list = ls())


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)


# Specifications ---------------------------------------------------------------
# Define publication 
author_year <- 'moore_2021'
# Define region abbreviation
region_abb  <- 'az'

# Directories
dir_publ        <- paste0(author_year, '_', region_abb, '/')
dir_data        <- paste0(dir_publ, "data/")
dir_data_ancill <- paste0(dir_data, 'Ancillary_Data_CSVs/')
dir_data_quad   <- paste0(dir_data, 'quadrat_data/')
dir_shp         <- paste0(dir_publ,  "data/Species_Shapefile_Extractions/")


# Data -------------------------------------------------------------------------

# Species 
sp_list <- read_csv(paste0(dir_data_ancill, "Plant_Species_List.csv"))

# Quadrat inventory
quad_inv <- read_csv(paste0(
  dir_data_ancill,"Summarize_Quadrats_by_Year.csv")) %>%
  select(-c('Years_Surveyed', 'Proportion', 'Comments')) %>%
  column_to_rownames("Quadrat") %>%   # Move the first column to row names
  t() %>%                             # Transpose the dataframe
  as.data.frame() %>%                 # Convert back to dataframe
  rownames_to_column("Year") %>%
  mutate(Year = substr(as.character(Year), 3, 4))

quad_inv[,-1] <- apply(quad_inv[,-1], 2, function(x) ifelse(x == "X", quad_inv$Year, x))


quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_sgs <- lapply(quadInv_list, function(x) as.numeric(x))

# Create a list of all shapefiles in the directory
shpFiles <- list.files(
  paste0(dir_shp, "/"),
  pattern = ".shp$",
  recursive = TRUE
)

shpFiles[1]


# ---------------------------------------------------------
# SF
for (n in 1:length(shpFiles)) {
  shape_x  <- paste0(shp_dir, shpFiles[n])
  shapeNow <- sf::st_read(dsn = shape_x, crs = NA)
  # Add Site, Quad, and Year
  # shapeNow$quad <- strsplit(shape_x, split = "/")[[1]][22]
  # shapeNow$year <- as.numeric(
  #   strsplit(strsplit(shape_x, split = "_")[[1]][23], split = '.')[[1]][1])
  if (n == 1) {
    dat <- shapeNow
  } else {
    dat <- rbind(dat, shapeNow)
  }
}

dat <- dat %>% 
  mutate(year = as.numeric(z_Year))

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, file = paste0(dir_data_quad, "moore21_quadrats_full.rds"))
dat <- readRDS(file = paste0(dir_data_quad, "moore21_quadrats_full.rds"))


# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "species", site = "Site", 
         quad = "Quadrat", year = "year", geometry = "geometry")

