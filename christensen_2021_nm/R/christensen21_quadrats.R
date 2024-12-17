# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# Aspen Workman, fall 2023
#
#
# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.
#
# Publication:  https://doi.org/10.1002/ecy.3530
#
# Setup
#
library(sf) #ver 1.0-1.2
# library(plantTracker) #ver 1.1.0
library(tidyverse)

base_dir <- ('christensen_2021_nm')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")
shp_dir <- paste(base_dir, "/data/quadrat_data/Jornada_shapefiles/", sep="")


# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read.csv(paste0(dat_dir, "species_list.csv")) %>% 
  mutate(species_bn2 = species,
         species = paste(genus, species)) %>%
  select(species, everything())

# sp_name_changes <- read.csv(paste0(dat_dir, "species_name_changes.csv")) 
#  will use to check names later on

# Read in quad inventory to use as 'inv' list in plantTracker
# Read the data
quad_inv <- read.csv(paste0(dat_dir, "Jornada_quadrat_sampling_dates.csv"))[c(1:2)] %>%
  mutate(quadrat = as.factor(quadrat))

quadInv_list <- split(quad_inv$year, quad_inv$quadrat)
inv_sgs <- quadInv_list

# Read in shapefiles to create sf dataframe to use as 'dat' in plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)
# Create list of quad names
quadNames <- list.files(paste0(dat_dir, 'Jornada_shapefiles'))
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
    shapeNow$Site <- "NM"
    shapeNow$Quad <- quadNow
    shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
    if (grepl(quadYearNow, pattern = "_pnt")) {
      # Keep only the relevant columns for point data
      shapeNow <- shapeNow %>%
        select(species, Site, Quad, Year, geometry) %>%
        mutate(type = "point")
      
      # Apply buffer if necessary for point data
      shapeNow <- st_buffer(shapeNow, dist = .0025)
    } else {
      # Keep only the relevant columns for polygon data
      shapeNow <- shapeNow %>%
        select(species, Site, Quad, Year, geometry) %>%
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
invalid_geom3 <- c(
  84097, 84276, 84278, 84284, 84287, 84788, 84790, 84797, 84891, 84983, 85272, 
  85308, 85432, 85722, 85915, 86004, 86006, 86136, 86855, 87068, 87679, 87700, 
  87851, 87980, 88221, 88296, 88600, 88704, 88705, 89198, 89209, 89250, 89325, 
  89828, 90245, 90560, 90908, 91531, 93498, 93810, 94125, 94231, 94421, 94516, 
  94537, 94774, 94957, 95176, 95177, 95590, 95593, 95625, 95659, 95681, 95705, 
  95722, 95948, 96222, 96339, 96876, 99074, 99278, 99528, 100124, 100136, 
  100774, 100807, 100833, 101630, 101737, 102260, 102600, 102605, 102630, 
  107938, 108186, 108274, 108644, 109236, 109762, 110054, 110428, 110810, 
  110939, 111117, 111211, 111543, 111666, 111770, 113269, 113853, 114078, 
  114299, 114300, 114343, 114691, 115359, 115535, 115554, 115560, 115689, 
  115904, 116045, 116069, 116103, 116217, 116304, 118893, 119007, 119208, 
  119220, 119562, 119792, 120299, 120302, 120726, 120873, 121031, 121034, 
  121039, 121411, 121461, 121544, 121691, 121805, 121887
)
for(i in 1:length(invalid_geom3)){
  dat01[invalid_geom3[i],6] <- st_make_valid(dat01[invalid_geom3[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom4 <- c(
  123325, 123790, 123984, 124009, 124070, 124334, 124681, 125120, 125129, 
  125839, 126288, 130048, 130659, 130935, 130961, 131130, 131281, 131344, 
  131349, 131406, 131475, 131606, 131674, 131676, 131751, 131762, 131841, 
  132125, 132232, 132253, 132258, 132548, 132609, 132623, 133163, 133274
)

for(i in 1:length(invalid_geom4)){
  dat01[invalid_geom4[i],6] <- st_make_valid(dat01[invalid_geom4[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
dat02 <- dat01[!is.na(dat01$Species), ]
checkDat(dat02, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  401, 1010, 1661, 24487, 25606, 25670, 26255, 32384, 38391, 39569, 55076,
  55089)

dat03 <- dat02[!(row.names(dat02) %in% drop_rows),]
checkDat(dat03, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows2 <- c(
  55112, 58091, 64803, 72238, 75507, 77291, 79196, 96626, 103604, 116003, 
  128419)

dat03 <- dat03[!(row.names(dat03) %in% drop_rows2),]
checkDat(dat03, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat03, file = paste0(dat_dir, "SGS_LTER_plantTracker_all_filtered.rds"))

