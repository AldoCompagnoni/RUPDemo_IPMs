# plantTracker for Adler 2007 Kansas Andropogon gerardii

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com
# Date: 2024.10.17

# Code adapted from: https://github.com/aestears/plantTracker
 # Adapted from plantTracker How to (Stears et al. 2022)

# List of required CRAN packages
.cran_packages <- c(
  # Spatial data handling
  'sf',
  # Plant tracking analysis
  'plantTracker') 

# Check if all required packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install any missing packages from CRAN
  install.packages(.cran_packages[!.inst]) 
}
# Load the necessary packages into the R environment
sapply(.cran_packages, require, character.only = TRUE)

# Define directories for data
data_directory   <- 'adler_2007_ks/data/'
quadrat_data_dir <- file.path(data_directory, 'quadrat_data/')
shape_files_dir  <- file.path(quadrat_data_dir, 'arcexport/')          

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
sp_list     <- read.csv(paste0(data_directory,
                               "quadrat_data/species_list.csv"))  %>% 
  dplyr::arrange(desc(count)) %>% head(25)
# Select the x_th species (target species)
target_spec <- sp_list %>% .[c(5),]  

# Define the species variable and abbreviation
species <- target_spec[1,1]
sp_abb  <- tolower(gsub(" ", "", paste(substr(
  strsplit(species, " ")[[1]], 1, 2), collapse = "")))

# Read in quadrat inventory data and prepare for plantTracker
quad_inv        <- 
  as.list(read.csv(paste0(data_directory,
                          "quadrat_data/quadrat_inventory.csv"), sep=',') %>% 
  dplyr::select(-year))
# Remove NAs
inv_ks          <- lapply(X = quad_inv, 
                          FUN = function(x) x[is.na(x) == FALSE])
# Replace dots in names with dashes
names(inv_ks)   <- gsub( '\\.','-',names(inv_ks) )  

# Read spatial data (polygon for each species per quadrat)
dat             <- readRDS(file=paste0(data_directory,
                                       "quadrat_data/KS_polygons_full.rds"))

# Subset data for the target species
dat_target_spec <- dat[dat$SCI_NAME %in% target_spec$species,] %>%
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
