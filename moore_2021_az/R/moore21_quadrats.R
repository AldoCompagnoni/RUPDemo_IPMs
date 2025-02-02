# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs

# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.

# Publication: https://doi.org/10.1002/ecy.3661

library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)

base_dir <- ('moore_2021_az/')
dat_dir <- paste(base_dir, "data/", sep="")
shp_dir <- paste(base_dir,  "data/Species_Shapefile_Extractions/", sep="")


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


# setwd(dat_dir)

# # Read in species list, species name changes, and subset species list to perennial grasses
# # with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# # species with the lowest cover later.
# sp_list <- read_csv(paste0(dat_dir, "Plant_Species_List.csv"))
# 
# # Read in quad inventory:
# #  a table that has all the quadrats as columns and the respective year 
# #   they've been sampled as rows with NA for years that are missing
# #  to use as 'inv' list in plantTracker
# quad_inv <- read_csv(paste0(dat_dir,"Summarize_Quadrats_by_Year.csv")) %>% 
#   select(-c('Years_Surveyed', 'Proportion', 'Comments')) %>%
#   column_to_rownames("Quadrat") %>%   # Move the first column to row names
#   t() %>%                             # Transpose the dataframe
#   as.data.frame() %>%                 # Convert back to dataframe
#   rownames_to_column("Year") %>%
#   mutate(Year = substr(as.character(Year), 3, 4))
# 
# quad_inv[,-1] <- apply(quad_inv[,-1], 2, function(x) ifelse(x == "X", quad_inv$Year, x))
# 
# 
# quadInv_list <- as.list(quad_inv)
# quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
# inv_sgs <- quadInv_list
# 
# 
# # Create a list of all shapefiles in the directory
# shpFiles <- list.files(shp_dir)
# 
# quadYears <- unlist(strsplit(list.files(
#   paste0(shp_dir,"/"),
#   pattern = ".shp$"), split = ".shp"))
# 
# 
# spec_names <- list.dirs(shp_dir, full.names = F, recursive = FALSE)
# for (j in 1:length(spec_names)) {
#   spec_name_now <- spec_names[j]
# }
# 
# 
# list.files(paste0(
#   shp_dir, 'Achillea_millefolium_Jan_8_2023'))
# 
# shapeNow <- sf::st_read(dsn = paste0(
#   shp_dir, "Achillea_millefolium_Jan_8_2023/Fry_Park_Jan_8_2023/Quadrat_30735_Jan_8_2023"),
#             layer = 'Quadrat_30735_2005')
# shapeNow$Site <- "AZ"
# 
# 
# 
# for (j in 1:length(quadYears)) {
#   quadYearNow <- quadYears[j]
#   shapeNow <- sf::st_read(dsn = paste0(shp_dir),
#                           layer = quadYearNow)
#   shapeNow$Site <- "Id"
#   shapeNow$Quad <- strsplit(quadYearNow, split = "_")[[1]][1]
#   shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
#   if (grepl(quadYearNow, pattern = "_D")) {
#     shapeNow <- shapeNow[,!(names(shapeNow)
#                             %in% c("OBJECTID", "seedling", "stem", "x", "y"))]
#     shapeNow <- sf::st_buffer(x = shapeNow, dist = .0025)
#     shapeNow$type <- "point"
#   } else {
#     shapeNow <- shapeNow[,!(names(shapeNow) %in% c("SP_ID", "stemID", "area", "x", "y"))]
#     
#     shapeNow$type <- "polygon"
#   }
#   if (j == 1) {
#     dat <- shapeNow
#   } else {
#     dat <- rbind(dat, shapeNow)
#   }
# } 


# Quadrat list -----------------------------------------------------------------
shapefiles <- list.files("C:/code/RUPDemo_IPMs/moore_2021_az/data/Species_Shapefile_Extractions/", 
                         pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)


# Function to extract year and quadrat information from file names
extract_info <- function(file_path) {
  # Extract the year (assuming it's a 4-digit number in the filename)
  year <- sub(".*_(\\d{4})\\.shp$", "\\1", basename(file_path))
  
  # Extract quadrat identifier (assuming it's before the year in the filename)
  quadrat <- sub("^(.*)_(\\d{4})\\.shp$", "\\1", basename(file_path))
  
  return(c(quadrat = quadrat, year = year))  # Return as a vector
}

# Apply the function to all shapefiles and store results
file_info <- lapply(shapefiles, extract_info)

# Convert to a data frame for easier manipulation
file_info_df <- do.call(rbind, file_info)
colnames(file_info_df) <- c("quadrat", "year")
file_info_df <- as.data.frame(file_info_df)

file_info_df <- file_info_df %>% 
  mutate(year = str_sub(as.factor(year))) %>%
  rename(Year = year) %>% 
  mutate(Year = as.numeric(Year)) %>%
  mutate(quadrat = gsub("Quadrat_", "", quadrat),
         quadrat = gsub("_bar_", " / ", quadrat))


unique_quadrats_by_plot <- bind_rows(
  tibble(
    quadrat = "year",  # New quadrat
    Year = list(unique(file_info_df$Year))),
  file_info_df %>%
  group_by(quadrat) %>%
  summarise(Year = list(unique(Year))) %>%
  arrange(quadrat)
  )

unique_year_list     <- setNames(unique_quadrats_by_plot$Year, 
                                 unique_quadrats_by_plot$quadrat)
inv_sgs <- unique_year_list


saveRDS(inv_sgs, 
        file = paste0(dir_data_quad, "moore21_quadrat_inventory.RData"))


# Data from geo-data -----------------------------------------------------------
wdName <- paste0(base_dir, "data/Quadrat_Spatial_Data/Combined_by_Site.gdb" )
"moore_2021_az/data/Quadrat_Spatial_Data"
cover_all <- sf::st_read( dsn = wdName,
                          layer = "Cover_All" ) %>% 
  rename(geometry = Shape,
         species  = species,
         Quadrat = Quadrat)

density_all <- sf::st_read( dsn = wdName,
                            layer = "Density_All" ) %>% 
  rename(geometry = Shape,
         species  = species,
         Quadrat = Quadrat) 


# Species list -----------------------------------------------------------------
names(cover_all)
summary(cover_all$area)
hist(log(cover_all$Shape_Area))


names(density_all)
density_all1 <- density_all %>% 
  mutate(area = NA,
         N_Flower = NA) %>%
  rename(x = coords_x1,
         y = coords_x2)


density_all2 <- st_centroid(density_all1) %>%  
  st_buffer(dist = 0.0003120444)

names(density_all2)
names(cover_all)
df0 <- rbind(density_all2, cover_all)

names(df0)
df0 <- df0 %>% 
  select(-c('area')) %>% 
  mutate(Year = as.numeric(Year))
  


# class(dplyr::select(cover_all,geometry))
# df<- rbind(dplyr::select(density_all,geometry), 
#            dplyr::select(cover_all,geometry))
# class(df)
# 
# dens <- st_drop_geometry( density_all )
# cov <- st_drop_geometry( cover_all )
# 
# colnames( dens ) <- c( "id", "species", "seedling", "x", "y", "site", "spcode", 
#                        "quadrat", "year", "type", "is_empty", "shape_length", 
#                        "shape_area" )
# 
# cov <- cov[,c(1,2,4:6,8:15)]
# colnames( cov ) <- c( "id", "species", "seedling", "x", "y", "site", "spcode", 
#                       "quadrat", "year", "type", "is_empty", "shape_length", 
#                       "shape_area" )
# 
# all <- rbind( dens, cov )
# 
# df3 <- st_set_geometry(all, df)
# 
# sf::st_join(df, all)
# class(df1)




  
saveRDS(df0, paste0(dat_dir, "moore21_quadrats_full.rds"))


# # Drop instances which we will not consider for modeling
# # e.g. taxa which were not identified to the species level
# 
# all <- all[-grep( " sp.", all$species ),]
# all <- all[-grep( "Unknown", all$species ),]
# all <- all[-which( all$species %in% c( "No Density Species Observed",
#                                        "No Cover Species Observed",
#                                        "Nama dichotoma" ) ),]


# Create summary table
## For each species, summarize the total number of quads, the number of years it
## was observed, and the total instances of observation

summary <- st_drop_geometry(df0) %>% 
  group_by(species, Type) %>% 
  summarise( quads = length( unique( Quadrat ) ),
                              years = length( unique( Year ) ),
                              counts = n()) %>%
  arrange(desc(counts))

write.csv(row.names = F, summary, 
          paste0(dir_data_quad, "moore21_species_list.csv"))



# # Save the output file so that it doesn't need to be recreated ever again
# saveRDS(dat, paste0(dat_dir, "moore21_quadrats_full.rds"))
# dat <- readRDS(paste0(dat_dir,"moore21_quadrats_full.rds"))

# Check the inv and dat arguments
checkDat(df0, inv_sgs, species = "species", site = "Site", quad = "Quadrat", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  128863, 134036, 144006, 145503, 153280, 159760, 160926, 175316, 187285, 191348, 191382, 191411, 196535, 196626, 196676, 196806, 196897, 199329, 200481, 205058, 206569, 207258, 207335, 207366, 207387, 210523, 215432, 216605, 218610, 220632, 221304, 223156, 228187, 230930, 230947, 236626, 239737, 240002, 242644, 242645, 242650, 242660, 242667, 242670, 242671, 242688)
dat01 <- df0
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],15] <- st_make_valid(dat01[invalid_geom[i],15])
}

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quadrat", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  8684, 37028, 76227, 76229, 76231, 76233, 76235, 76237, 76239)

dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat02, inv_sgs, species = "species", site = "Site", quad = "Quadrat", year = "Year", geometry = "geometry")



# FIX----------------------------------------------------------------------------

saveRDS(dat02, file = paste0(dir_data_quad, "moore21_quadrats_full.rds"))
dat02 <- readRDS(file = paste0(dat_dir, "moore21_quadrats_filtered.rds"))
