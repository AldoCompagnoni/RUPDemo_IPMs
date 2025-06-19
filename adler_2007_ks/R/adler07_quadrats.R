# Quadrat data - Moore 2021 Arizona 

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.04.15

# Publication: https://doi.org/10.1890/0012-9658(2007)88[2673:LMQFKP]2.0.CO;2

# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)
library(janitor)


# Directories ------------------------------------------------------------------
dir_pub <- file.path('adler_2007_ks')
dir_dat <- file.path(dir_pub, 'data')
dir_qud <- file.path(dir_dat, 'quad_data')
dir_shp <- file.path(dir_qud, 'arcexport')


# Key variables ----------------------------------------------------------------
v_delimiter <- c(',')


# Data -------------------------------------------------------------------------
# Species list
sp_list <- read_delim(
  file.path(dir_qud, "species_list.csv"), 
  delim = v_delimiter, escape_double = FALSE, trim_ws = TRUE) %>% 
  clean_names()
  
# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv      <- read_delim(
  file.path(dir_qud, "quadrat_inventory.csv"),
  delim = v_delimiter, escape_double = FALSE, trim_ws = TRUE)
quadInv_list  <- as.list(quad_inv)
quadInv_list  <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_ks <- quadInv_list
names(inv_ks) <- gsub( '\\.','-',names(inv_ks) )

# Read in shapefiles to create sf dataframe to use as 'dat' in plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)
# Create list of quad names
quadNames <- list.files(dir_shp)

# Use for loop to download data from each quad folder
for(i in 1:length(quadNames)){
  
  quadNow      <- quadNames[i]
  # quad/year name combination
  quadYears    <- file.path(dir_shp, quadNow) %>% 
    list.files( pattern = ".e00$" ) 
  
  # loop over each year  
  for(j in 1:length(quadYears)){
    
    quad_yr_name  <- quadYears[j]
    shapeNow      <- sf::st_read( 
      dsn   = file.path(dir_shp, quadNow, quad_yr_name ),
      layer = 'PAL' ) %>% 
      dplyr::select( SCI_NAME, geometry) %>% 
      clean_names() %>% 
      rename(species = sci_name)
    shapeNow$site <- "ks"
    shapeNow$quad <- quadNow
    shapeNow$year <- quad_yr_name %>% 
      gsub(quadNow,'',.) %>% 
      gsub('.e00','',.) %>%
      unlist %>% 
      as.numeric()
    
    # start final data frame
    if(i == 1 & j == 1) {
      dat <- shapeNow
      # "append" data to the initial data frame  
    } else {
      dat <- rbind(dat, shapeNow)
    }
  }
}

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat,   file.path(dir_qud, 'adler07_quadrats_full.rds'))
dat <- readRDS(file.path(dir_qud, 'adler07_quadrats_full.rds'))


# Clean data -------------------------------------------------------------------
# Check the inv and dat arguments
checkDat(dat, inv_ks, species = "species", site = "site", quad = "quad", 
         year = "year", geometry = "geometry")

dat01 <- dat

# Save the data
saveRDS(dat01, file.path(dir_qud,'adler07_quadrats_filtered.rds'))
