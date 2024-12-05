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
# library(plantTracker) #ver 1.1.0
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
    shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][4])
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
  2311, 2347, 2649, 2929, 3229, 3301, 3874, 3877, 3879, 3933, 4008, 4012, 4163, 
  4248, 4275, 4305, 4480, 4509, 4890, 5020, 6461, 6858, 6892, 7025, 7094, 7118, 
  7322, 7634, 7684, 8187, 8315, 8397, 8415, 8521, 8537, 8742, 9350, 9513, 10506, 
  10836, 11043, 11117, 11121, 11153, 11464, 11487, 11822, 11849, 11906, 11926, 
  12160, 12200, 12254, 12704, 13489, 14743, 15185, 15647, 15650, 15892, 15978, 
  16376, 16499, 16906, 17117, 17583, 18024, 18347, 18374, 18406, 18427, 18778, 
  18786, 19427, 20258, 20665, 21173, 21340, 21818, 22169, 22202, 22372, 23128, 
  23155, 23478, 23950, 24536, 25893, 25935, 27171, 27523, 27775, 28132, 28525, 
  28529, 29028, 29045, 29052, 29071, 29102, 29105, 29205, 29539, 29943, 30065, 
  30219, 30635, 31042, 31617, 31634, 34811, 35096, 35388, 35683, 35948, 36047, 
  36139, 36174, 36194, 36211, 36219, 36220, 36552, 36582, 36738, 36897, 37347, 
  37385, 37388, 37471, 39161, 39549, 39733, 39927, 40093, 40315, 40324, 40459, 
  40783, 40834, 41218, 41226, 41252, 41258, 41324, 41330, 41502, 41602, 41645, 
  42148, 42336, 43659, 44177, 44285, 44376, 44593, 44756, 44778, 44847, 44892, 
  44894, 44962, 45351, 45515, 45559, 45652, 45682, 47282, 47496, 47577, 47904, 
  48020, 48378, 48764, 49028, 49057, 49097, 49159, 49341, 49370, 49440, 49603, 
  50608, 50703, 51166, 51175, 51377, 51402, 51557, 51683, 51849, 52111, 52454, 
  52457, 53041, 53448, 53449, 53471, 53523, 54208, 55464, 55764, 55872, 55905, 
  56139, 56314, 56496, 56826, 57085, 57086, 57156, 57183, 57508, 57766, 57850, 
  57975, 58194, 58795, 59964, 60244, 60261, 60906, 61559, 61957, 61962, 62564, 
  63315, 63419, 63736, 63739, 63764, 63808, 63809, 64229, 64320, 65145, 66132, 
  68407, 68586, 68878, 68922, 69033, 69300, 69315, 69401, 69827, 69983, 70051, 
  70375, 70435, 70468, 70470, 70493, 70503, 70552, 70663, 70912, 71118, 74592, 
  75064, 75067, 75433, 75800, 75962, 76001, 76473, 77872, 78722, 79486, 79550, 
  80555, 82966, 83268, 83403, 83447, 83746)
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

