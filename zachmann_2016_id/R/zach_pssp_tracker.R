# plantTracker for chu 2013 colorado Pseudoroegneria spicata

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com
# Date: 2024.10.24

# Code adapted from: https://github.com/aestears/plantTracker
 # Adapted from plantTracker How to (Stears et al. 2022)

rm(list = ls())

# Packages ---------------------------------------------------------------------
# Load packages, verify, and download if needed
source('helper_functions/load_packages.R')
load_packages(sf, plantTracker, tidyverse)


# Data -------------------------------------------------------------------------
# Directory
base_dir <- ('zachmann_2016_id')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")

# Function to quote bare names for tidy evaluation
quote_bare <- function(...){
  # creates a list of the arguments, but keeps their expressions unevaluated
  # wraps the list, original expressions (the bare names) without evaluation
  substitute(alist(...)) %>% 
    # evaluates the substituted list, which results in a list of the bare name
    eval() %>% 
    # converts the list into their character representations
    # apply to each element, returning a character vector of the names
    sapply(deparse)  
}


# Read species list and filter for target species
sp_list     <- read_csv(paste0(dat_dir, "/species_list.csv")) %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(cover)) 

# Select the x_th species (target species)
target_spec <- sp_list %>% .[c(1),]  

# Define the species variable and abbreviation
species <- target_spec[1,1]
sp_abb  <- tolower(gsub(" ", "", paste(substr(
  strsplit(species, " ")[[1]], 1, 2), collapse = "")))

# Read in quadrat inventory data and prepare for plantTracker
quad_inv        <- 
  as.list(read_csv(paste0(dat_dir, "/quad_inventory.csv")) %>% 
  as.data.frame())

# Remove NAs
inv_ks          <- lapply(X = quad_inv, 
                          FUN = function(x) x[is.na(x) == FALSE])
# Replace dots in names with dashes
names(inv_ks)   <- gsub( '\\.','-',names(inv_ks) )  

# Read spatial data (polygon for each species per quadrat)
dat             <- readRDS(
  file=paste0(dat_dir, "/SGS_LTER_plantTracker_all_filtered.rds"))

# Subset data for the target species
dat_target_spec <- dat[dat$species %in% target_spec$species,] %>%
  select(-c(type, Species)) %>% 
  # Rename columns
  setNames(quote_bare(Species, Site, Quad, Year, geometry))  

# Prepare data for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec, inv_ks,
                        # Dormancy flag
                        dorm         = 1,
                        # Buffer size
                        buff         = 5,
                        # Allow for clonal tracking
                        clonal       = TRUE,
                        # Buffer for genet
                        buffGenet    = 5,
                        # Aggregate by genet
                        aggByGenet   = TRUE,
                        # Flag potential issues
                        flagSuspects = TRUE)   

# Create output directory if it doesn't exist
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))
}

# Save the tracked data to a CSV file
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select(Site, Quad, Species, trackID, Year, basalArea_genet, recruit,
                survives_tplus1, age, size_tplus1, nearEdge, Suspect) %>%
  write.csv(paste0('adler_2007_ks/data/', sp_abb, '/ks_', sp_abb, '.csv'), 
            row.names = FALSE)  
