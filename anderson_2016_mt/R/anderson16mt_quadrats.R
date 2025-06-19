# Quadrat data - Anderson 2016 Montana

# Author: Niklas Neisse
# Co    : Diāna Spurīte, Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.06.19

# Publication: https://doi.org/10.1890/11-0193.1


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)
library(janitor)


# Directories ------------------------------------------------------------------
dir_pub <- file.path('anderson_2016_mt')
dir_dat <- file.path(dir_pub, 'data')
dir_qud <- file.path(dir_dat, 'quadrat_data')
dir_shp <- file.path(dir_qud, 'shapefiles')


# Data -------------------------------------------------------------------------
# Species list
sp_list <- read.csv(file.path(dir_qud, "species_list.csv"))

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read.csv(file.path(dir_qud, "quad_inventory.csv"))
quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_sgs <- quadInv_list


# Create a list of all shapefiles in the directory
shpFiles <- list.files(dir_shp)

quadYears <- unlist(strsplit(list.files(
  file.path(dir_shp),
  pattern = ".shp$"), split = ".shp"))

for (j in 1:length(quadYears)) {
  quadYearNow <- quadYears[j]
  
  # Read the shapefile
  shapeNow <- sf::st_read(dsn = file.path(dir_shp), 
                          layer = quadYearNow, 
                          # Because of some custom corrdi ref system,
                          #  might aswell use: st_crs(4326)
                          crs = NA)
  
  # Add Site, Quad, and Year
  shapeNow$Site <- "Mt"
  shapeNow$Quad <- strsplit(quadYearNow, split = "_")[[1]][1]
  shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
  
  # Handle specific columns based on layer type (_D or not)
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
  
  # Combine all data into the `dat` variable
  if (j == 1) {
    dat <- shapeNow
  } else {
    dat <- rbind(dat, shapeNow)
  }
} 

dat <- dat %>% clean_names()

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat,   file = file.path(dir_qud, "anderson16mt_quadrats_full.rds"))
dat <- readRDS(file = file.path(dir_qud, "anderson16mt_quadrats_full.rds"))

# Check the inv and dat arguments ----------------------------------------------
checkDat(dat, inv_sgs, species = "species", site = "site", quad = "quad", 
         year = "year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  452, 481, 686, 713, 2224, 3287, 3736, 18640, 19211, 36792, 42830, 44822, 
  45384, 46597, 46604, 46611, 47825, 49400, 51516, 51518, 58996, 59722, 59999, 
  71625, 76288, 76471, 76694, 76773, 82224, 82507, 82520, 92972, 100890, 103095, 
  103189, 103396, 107938, 107946, 108262, 108275, 108669, 109006, 109482, 
  110340, 110519, 110622, 110643, 110948, 111336, 113642, 114459, 115375, 
  115444, 115883, 115959, 116622, 120225, 120230, 121366, 121895, 122632, 
  124332, 124499, 125595, 126685, 130351, 130352, 130691, 134742, 144658, 
  144661, 144679, 148967, 149133, 150049, 151646, 153203, 154089, 154332, 
  157148, 162078, 163308, 164625, 171411, 175206, 175208, 177238, 177456, 
  180368, 182126, 182928, 185677, 185702, 186382, 188348, 189072, 189088, 
  189117, 189426, 189448, 189473, 189864, 190134, 190152, 190296, 190474, 
  203736, 204777)

dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
 }

checkDat(dat01, inv_sgs, species = "species", site = "site", quad = "quad", 
         year = "year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  20522, 26373, 34704, 67889, 74969, 132988, 135544, 165470, 167341)

dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat02, inv_sgs, species = "species", site = "site", quad = "quad", 
         year = "year", geometry = "geometry")


# Save data
saveRDS(  dat02, file = file.path(dir_qud, "anderson16mt_quadrats_filtered.rds"))
dat02 <- readRDS(file = file.path(dir_qud, "anderson16mt_quadrats_filtered.rds"))