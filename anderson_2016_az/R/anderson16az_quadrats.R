# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# Aspen Workman, fall 2023
#
#
# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.
#
# Publication: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3530
#
# Setup
#
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)

base_dir <- ('anderson_2016_az')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")
shp_dir <- paste(base_dir, "/data/quadrat_data/shapefiles/", sep="")


# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read.csv(paste0(dat_dir, "species_list.csv"))

sp_name_changes <- read.csv(paste0(dat_dir, "species_name_changes.csv")) 
#  will use to check names later on

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read.csv(paste0(dat_dir, "quad_inventory.csv"))
quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_sgs <- quadInv_list

# Read in shapefiles to create sf dataframe to use as 'dat' in plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)
# Create list of quad names
quadNames <- list.files(paste0(dat_dir, 'shapefiles'))
# Use for loop to download data from each quad folder
for(i in 1:length(quadNames)){
  quadNow <- quadNames[i]
  quadYears <- unlist(strsplit(list.files(
    paste0(shp_dir,quadNow,"/"),
    pattern = ".shp$"), split = ".shp"))
  for (j in 1:length(quadYears)) {
    quadYearNow <- quadYears[j]
    shapeNow <- sf::st_read(dsn = paste0(shp_dir,quadNow),
                            layer = quadYearNow)
    shapeNow$Site <- "AZ"
    shapeNow$Quad <- quadNow
    shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
    if (grepl(quadYearNow, pattern = "_D")) {
      # Keep only the relevant columns for point data
      shapeNow <- shapeNow %>%
        select(Species, Site, Quad, Year, geometry) %>%
        mutate(type = "point")
      
      # Apply buffer if necessary for point data
      shapeNow <- st_buffer(shapeNow, dist = .0025)
    } else {
      # Keep only the relevant columns for polygon data
      shapeNow <- shapeNow %>%
        select(Species, Site, Quad, Year, geometry) %>%
        mutate(type = "polygon")
    }
    if (i == 1 & j == 1) {
      dat <- shapeNow
    } else {
      dat <- rbind(dat, shapeNow)
    }
  }
}

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, file = paste0(dat_dir, "anderson16az_quadrats_full.rds"))
dat <- readRDS(file = paste0(dat_dir, "anderson16az_quadrats_full.rds"))


# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  17707, 17838, 18256, 19466, 21327, 35642, 39386, 46376, 60034, 60042, 60364, 
  60374, 97081, 100701, 111066, 113144, 120701, 128870, 128929, 129569, 130959, 
  134738)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  61364, 74426, 74463, 74467)

dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat02, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat02, file = paste0(dat_dir, "anderson16az_quadrats_filtered.rds"))

