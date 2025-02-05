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
  mutate(Year = as.numeric(str_sub(Year, start = -2))) %>%
  mutate(quadrat = gsub("Quadrat_", "", quadrat),
         quadrat = gsub("_bar_", " / ", quadrat))


unique_quadrats_by_plot <- file_info_df %>%
  group_by(quadrat) %>%
  summarise(Year = list(unique(Year))) %>%
  arrange(quadrat)

unique_year_list     <- setNames(unique_quadrats_by_plot$Year, 
                                 unique_quadrats_by_plot$quadrat)
inv_sgs <- unique_year_list



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


dens <- st_drop_geometry( density_all )
cov <- st_drop_geometry( cover_all )

colnames( dens ) <- c( "id", "species", "seedling", "x", "y", "site", "spcode", 
                       "quadrat", "year", "type", "is_empty", "shape_length", 
                       "shape_area" )

cov <- cov[,c(1,2,4:6,8:15)]
colnames( cov ) <- c( "id", "species", "seedling", "x", "y", "site", "spcode", 
                      "quadrat", "year", "type", "is_empty", "shape_length", 
                      "shape_area" )

all <- rbind( dens, cov )


# Drop instances which we will not consider for modeling
# e.g. taxa which were not identified to the species level

all <- all[-grep( " sp.", all$species ),]
all <- all[-grep( "Unknown", all$species ),]
all <- all[-which( all$species %in% c( "No Density Species Observed",
                                       "No Cover Species Observed",
                                       "Nama dichotoma" ) ),]


# Create summary table
## For each species, summarize the total number of quads, the number of years it
## was observed, and the total instances of observation

summary <- all %>% 
  group_by(species, type) %>% 
  summarise( quads = length( unique( quadrat ) ),
                              years = length( unique( year ) ),
                              counts = n()) %>%
  arrange(desc(counts))

write.csv(row.names = F, summary, 
          paste0(dat_dir, "quadrat_data/moore21_species_list.csv"))
saveRDS(inv_sgs, 
        file = paste0(dat_dir, "quadrat_data/moore21_quadrat_inventory.RData"))



# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, paste0(dat_dir, "moore21_quadrats_full.rds"))
dat <- readRDS(paste0(dat_dir,"moore21_quadrats_full.rds"))

# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  3146, 8319, 18289, 19786, 27563, 34043, 35209, 49599, 61568, 65631, 65665, 
  65694, 70818, 70909, 70959, 71089, 71180, 73612, 74764, 79341, 80852, 81541, 
  81618, 81649, 81670, 84806, 89715, 90888, 92893, 94915, 95587, 97439, 102470, 
  105213, 105230, 110909, 114020, 114285, 116927, 116928, 116933, 116943, 
  116950, 116953, 116954, 116971)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  132536, 160880, 200079, 200081, 200083, 200085, 200087, 200089, 200091)

dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat02, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat02, file = paste0(dat_dir, "moore21_quadrats_filtered.rds"))
dat02 <- readRDS(file = paste0(dat_dir, "moore21_quadrats_filtered.rds"))
