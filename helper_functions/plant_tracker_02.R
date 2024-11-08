# Second part of the plant tracker pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)

# Define the species variable and abbreviation
species <- target_spec[1,1]
sp_abb  <- tolower(gsub(' ', '', paste(substr(
  strsplit(species, ' ')[[1]], 1, 2), collapse = '')))

# Directory 2
pub_dir    <- file.path(paste0(author_year, '_', region_abb))
R_dir      <- file.path(pub_dir, 'R',      sp_abb)
data_dir   <- file.path(pub_dir, 'data',   sp_abb)
result_dir <- file.path(pub_dir, 'result', sp_abb)

# Read in quadrat inventory data and prepare for plantTracker
quad_inv        <- 
  as.list(read_csv(paste0(dat_dir, '/quad_inventory.csv')) %>% 
            as.data.frame())

# Remove NAs
inv          <- lapply(X = quad_inv, 
                          FUN = function(x) x[is.na(x) == FALSE])
# Replace dots in names with dashes
names(inv)   <- gsub( '\\.','-',names(inv) )  

# Read spatial data (polygon for each species per quadrat)
dat             <- readRDS(
  file=paste0(dat_dir, '/SGS_LTER_plantTracker_all_filtered.rds'))

# Subset data for the target species
dat_target_spec <- dat[dat$species %in% target_spec$species,] %>%
  select(-c(type, Species)) %>% 
  # Rename columns
  setNames(quote_bare(Species, Site, Quad, Year, geometry))  

buff <- if_else(st_bbox(dat_target_spec)[3] < 1.1, 0.05, 5)

# Prepare data for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec, inv,
                        # Dormancy flag
                        dorm         = 1,
                        # Buffer size
                        buff         = buff,
                        # Allow for clonal tracking
                        clonal       = TRUE,
                        # Buffer for genet
                        buffGenet    = buff,
                        # Aggregate by genet
                        aggByGenet   = TRUE,
                        # Flag potential issues
                        flagSuspects = TRUE)   

# Create output directory if it doesn't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

# Save the tracked data to a CSV file
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select(Site, Quad, Species, trackID, Year, basalArea_genet, recruit,
                survives_tplus1, age, size_tplus1, nearEdge, Suspect) %>%
  write.csv(paste0(data_dir, '/', script_prefix, '_', sp_abb, '.csv'), 
            row.names = FALSE)  
