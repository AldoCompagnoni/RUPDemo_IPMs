# plantTracker for Adler 2007 Kansas Hedyotis nigricans

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com
# Date: 2024.10.28

# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)

# load packages
source( 'helper_functions/load_packages.R' )
load_packages( sf, plantTracker, tidyverse )

# Type of population model
mod_type <- 'mpm'

# Directories 1
data0_dir <- 'adler_2007_ks/data/'
quad_dir  <- paste0( data0_dir, 'quadrat_data/')
sh_dir    <- paste0( quad_dir,  'arcexport/')          

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
forb_list     <- read.csv(paste0(data0_dir,
                                 "quadrat_data/species_list.csv"))  %>% 
  dplyr::arrange(desc(count)) %>% 
  filter(type == 'forb') %>% 
  head(25)

# Select the x_th species (target species)
target_spec <- forb_list %>% .[c(8),]  

# Define the species variable and abbreviation
species <- target_spec$species
sp_abb  <- tolower(gsub(" ", "", paste(substr(
  strsplit(species, " ")[[1]], 1, 2), collapse = "")))

# Directories 2
data_dir   <- paste0('adler_2007_ks/data/',    sp_abb, '_', mod_type)
result_dir <- paste0('adler_2007_ks/results/', sp_abb, '_', mod_type)
R_dir      <- paste0('adler_2007_ks/R/',       sp_abb, '_', mod_type)


# Read in quadrat inventory data and prepare for plantTracker
quad_inv        <- read.csv(paste0(data0_dir,
                                   "quadrat_data/quadrat_inventory.csv"), 
                            sep=',') %>% 
  dplyr::select(-year) %>% 
  as.list

# Remove NAs
inv_ks          <- lapply(X = quad_inv, 
                          FUN = function(x) x[is.na(x) == FALSE])

# Replace dots in names with dashes
names(inv_ks)   <- gsub( '\\.','-',names(inv_ks) )  

# Read spatial data (polygon for each species per quadrat)
dat             <- readRDS(file=paste0(quad_dir,
                                       "KS_polygons_full.rds"))

# Subset data for the target species
dat_target_spec <- dat[dat$SCI_NAME %in% target_spec$species,] %>%
  # Rename columns
  setNames( quote_bare(Species, Site, 
                       Quad, Year, geometry) )  

# Prepare data for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec, 
                        inv_ks,
                        # Number of years in dormancy
                        dorm         = 0,
                        # Buffer size
                        buff         = 5,
                        # Allow for clonal tracking
                        clonal       = FALSE,
                        # Buffer for genet
                        buffGenet    = 5,
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
  write.csv(paste0(data_dir, '/ks_', sp_abb, '.csv'), 
            row.names = FALSE)  
