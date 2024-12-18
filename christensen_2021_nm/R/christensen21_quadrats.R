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

# rm(list = ls())

library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)
library(stringr)
library(dplyr)

base_dir <- ('christensen_2021_nm')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")
shp_dir <- paste(base_dir, "/data/quadrat_data/Jornada_shapefiles/", sep="")


# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read.csv(paste0(dat_dir, "Jornada_quadrat_species_list.csv")) %>% 
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


# Set the path to your main folder
main_folder <- "christensen_2021_nm/data/quadrat_data/Jornada_shapefiles//"

# List all subfolders within the main folder
folder_names <- list.dirs(main_folder, full.names = FALSE, recursive = FALSE)

# Create an empty list to store the years for each folder
years_list_per_folder <- list()

# Iterate over each folder
for (folder in folder_names) {
  
  # List all .shp files in the current folder
  folder_path <- file.path(main_folder, folder)
  files <- list.files(folder_path, pattern = "\\.shp$", full.names = TRUE)
  
  # Initialize an empty vector to hold all years for the current folder
  folder_years <- c()
  
  # Iterate over the .shp files in the folder
  for (file in files) {
    # Print the filename to debug
    print(paste("Processing file:", basename(file)))
    
    # Extract the year from the filename (assuming it's a 4-digit year)
    years <- str_extract(basename(file), "\\d{4}")
    
    # If a year is found, add it to the folder_years vector
    if (!is.na(years)) {
      print(paste("Year extracted:", years))
      folder_years <- c(folder_years, years)
    } else {
      print("No year found in filename")
    }
  }
  
  # Store the unique years for this folder in the list
  if (length(folder_years) > 0) {
    years_list_per_folder[[folder]] <- unique(folder_years)
  } else {
    years_list_per_folder[[folder]] <- NA
  }
}

for (folder in names(years_list_per_folder)) {
  # Convert the years from character to integer
  years_list_per_folder[[folder]] <- as.integer(years_list_per_folder[[folder]])
}


# Now we need to create a dataframe with rownames as all unique years across all folders
# Find all unique years across all folders
all_unique_years <- sort(unique(unlist(years_list_per_folder)))

# Create an empty dataframe with all unique years as row names
years_df <- data.frame(matrix(ncol = length(folder_names), nrow = length(all_unique_years)))
rownames(years_df) <- all_unique_years
colnames(years_df) <- folder_names

# Fill the dataframe: For each folder, check if the year is present
for (folder in folder_names) {
  # Get the unique years for the current folder
  folder_years <- years_list_per_folder[[folder]]
  
  # For each unique year, fill the dataframe with the year if it's in the folder's unique years
  years_df[[folder]] <- ifelse(rownames(years_df) %in% folder_years, rownames(years_df), NA)
}

write.csv(row.names = F, years_df, paste0(dat_dir, "/quadrat_inventory.csv"))

inv_sgs <- years_list_per_folder

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
saveRDS(dat, file = paste0(dat_dir, "christensen21_quadrats_full.rds"))
dat <- readRDS(file = paste0(dat_dir, "/christensen21_quadrats_full.rds")) %>% 
  rename(Species = species)

dat$geometry <- st_make_valid(dat$geometry)
invalid_polygons <- st_is_valid(dat$geometry) & !sapply(dat$geometry, function(x) length(x[[1]]) >= 4)
dat <- dat[!invalid_polygons, ]
dat <- dat[!is.na(dat$Species), ]


# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  215388)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  215388)

dat03 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat03, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")



# Perform a left join to match Species in dat03 with species_code in sp_list
dat04 <- dat03 %>%
  left_join(sp_list %>%
              select(species_code, species), by = c("Species" = "species_code")) %>%
  mutate(Species = coalesce(species, Species)) %>%  # Update Species with species from sp_list
  select(-species)  # Optionally drop the extra species column

species_summary <- dat %>%
  select(Species) %>% 
  as.data.frame(.) %>% 
  count(Species) %>%
  rename(count = n) %>%
  arrange(desc(count))

sp_list <- sp_list %>%
  left_join(species_summary, by = c("species" = "Species"))

write.csv(row.names = F, sp_list, paste0(dat_dir, "species_list.csv"))

saveRDS(dat04, file = paste0(dat_dir, "christensen21_quadrats_filtered.rds"))
dat04 <- readRDS(file = paste0(dat_dir, "christensen21_quadrats_filtered.rds"))
