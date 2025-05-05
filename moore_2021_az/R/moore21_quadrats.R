# Quadrat data - Moore 2021 Arizona 

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.04.15

# Publication: https://doi.org/10.1002/ecy.3661

# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)
library(janitor)


# Specifications ---------------------------------------------------------------
# Define publication 
v_author_year <- c('moore_2021')
# Define region abbreviation
v_region_abb  <- c('az')


# Directories ------------------------------------------------------------------
dir_publ        <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_data        <- file.path(dir_publ, 'data')
dir_data_ancill <- file.path(dir_data, 'Ancillary_Data_CSVs')
dir_data_quad   <- file.path(dir_data, 'quadrat_data')
dir_shp         <- file.path(dir_publ, 'data', 'Species_Shapefile_Extractions')
dir_geo         <- file.path(dir_data, 'Quadrat_Spatial_Data', 'Combined_by_Site.gdb')


# Quadrat list -----------------------------------------------------------------
#  List of all quadrats and their respecitve years

# List of all shapefiles in the moore folder
shapefiles <- list.files(
  dir_shp, pattern = '\\.shp$', recursive = TRUE, full.names = TRUE)

# Function to extract year and quadrat information from file names
extract_info <- function(file_path) {
  # Extract the year (assuming it's a 4-digit number in the filename)
  year <- sub('.*_(\\d{4})\\.shp$', '\\1', basename(file_path))
  
  # Extract quadrat identifier (assuming it's before the year in the filename)
  quadrat <- sub('^(.*)_(\\d{4})\\.shp$', '\\1', basename(file_path))
  
  return(c(quadrat = quadrat, year = year))  # Return as a vector
}

# Apply the function to all shapefiles and store results
file_info <- lapply(shapefiles, extract_info)

# Convert to a data frame for easier manipulation
file_info_df <- do.call(rbind, file_info)
colnames(file_info_df) <- c('quadrat', 'year')
file_info_df <- as.data.frame(file_info_df)

file_info_df <- file_info_df %>% 
  mutate(year = str_sub(as.factor(year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(quadrat = gsub('Quadrat_', '', quadrat),
         quadrat = gsub('_bar_', ' / ', quadrat))

unique_quadrats_by_plot <- bind_rows(
  tibble(
    quadrat = 'year',  # New quadrat
    year = list(unique(file_info_df$year))),
  file_info_df %>%
  group_by(quadrat) %>%
  summarise(year = list(unique(year))) %>%
  arrange(quadrat)
  )

unique_year_list <- setNames(
  unique_quadrats_by_plot$year, unique_quadrats_by_plot$quadrat)
inv_sgs <- unique_year_list

saveRDS(
  inv_sgs, file = file.path(dir_data_quad, 'moore21_quadrat_inventory.RData'))


# Data from geo-data -----------------------------------------------------------
cover_all <- sf::st_read( 
  dsn = dir_geo, layer = 'Cover_All') %>% 
  rename(geometry = Shape) %>% 
  clean_names()

density_all <- sf::st_read(
  dsn = dir_geo, layer = 'Density_All' ) %>% 
  rename(geometry = Shape) %>% 
  clean_names()


# Centroids --------------------------------------------------------------------
# Some of the data has area sizes of -14 and lower. We might want to change that
hist(log(cover_all$shape_area))

# Centroids for Cover data
# Subset the data where Shape_Area < 2.5e-6 and Shape_Area > 0
subset_to_update <- cover_all[
  cover_all$shape_area < 2.5e-6 & cover_all$shape_area > 0, ]

# Apply centroid and buffer transformations to the subset
subset_to_update <- st_buffer(st_centroid(subset_to_update), dist = 0.001)

# Remove the updated subset from the original data
cover_all_remaining <- cover_all[
  cover_all$shape_area >= 2.5e-6 | cover_all$shape_area == 0, ]

# Combine the updated subset with the remaining original data
cover_all_updated <- rbind(cover_all_remaining, subset_to_update)


# Centroids for density data
density_all1 <- density_all %>% 
  mutate(area = NA,
         n_flower = NA) %>%
  rename(x = coords_x1,
         y = coords_x2)

density_all2 <- st_centroid(density_all1) %>%  
#  st_buffer(dist = 0.0003120444) # what we decided with aspen
  st_buffer(dist = 0.001) # this gives us the lowest in the graphs 
 # which is 3.14*10^-6 m2 -- log --> ~ -13

# Bind density and cover
df0 <- rbind(density_all2, cover_all_updated) %>% 
  select(-c('area')) %>% 
  mutate(year = as.numeric(year))

saveRDS(df0, file.path(dir_data_quad, 'moore21_quadrats_full.rds'))


# Species list -----------------------------------------------------------------
## For each species, summarize the total number of quads, the number of years it
## was observed, and the total instances of observation

df_species <- st_drop_geometry(df0) %>%
  group_by(species, type) %>% 
  summarise(
    quads  = length(unique(quadrat)),
    years  = length(unique(year)),
    counts = n()) %>%
  arrange(desc(counts))

write.csv(row.names = F, df_species, 
          file.path(dir_data_quad, 'moore21_species_list.csv'))


# Filter data ------------------------------------------------------------------
# Check the inv and dat arguments
checkDat(df0, inv_sgs, species = 'species', site = 'site', quad = 'quadrat', 
         year = 'year', geometry = 'geometry')
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  128861, 134027, 143973, 145470, 153222, 159612, 160762, 175001, 186821, 
  190840, 190874, 190903, 196026, 196117, 196167, 196297, 196388, 198730, 
  199882, 204361, 205862, 206551, 206628, 206659, 206680, 209765, 214580, 
  215726, 217726, 219736, 220399, 222227, 227217, 229959, 229976, 235643, 
  238753, 239018, 241560, 241561, 241566, 241576, 241583, 241586, 241587, 
  241604)

dat01 <- df0

for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],15] <- st_make_valid(dat01[invalid_geom[i],15])
}

checkDat(dat01, inv_sgs, species = 'species', site = 'site', quad = 'quadrat', 
         year = 'year', geometry = 'geometry')
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  8684, 37028, 76227, 76229, 76231, 76233, 76235, 76237, 76239)

dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]

checkDat(dat02, inv_sgs, species = 'species', site = 'site', quad = 'quadrat', 
         year = 'year', geometry = 'geometry')


# Save -------------------------------------------------------------------------
dat03 <- dat02 %>% 
  rename(location = site) %>% 
  mutate(site = 'az')


saveRDS(dat03, file = file.path(dir_data_quad, 'moore21_quadrats_filtered.rds'))
dat03 <- readRDS(file = file.path(dir_data_quad, 'moore21_quadrats_filtered.rds'))
