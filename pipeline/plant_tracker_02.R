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
if (gr_form == 'forb') {
  folder_suffix <- paste0(sp_abb, '_mpm')
  R_dir      <- file.path(pub_dir, 'R', folder_suffix)
  data_dir   <- file.path(pub_dir, 'data', folder_suffix)
  result_dir <- file.path(pub_dir, 'results', folder_suffix)
} else {
  R_dir      <- file.path(pub_dir, 'R', sp_abb)
  data_dir   <- file.path(pub_dir, 'data', sp_abb)
  result_dir <- file.path(pub_dir, 'results', sp_abb)
}

# Identify quadrat inventory file 
quad_file      <- list.files(
                    path    = dat_dir, full.names = TRUE,
                    pattern = 'quad_inventory|quadrat_inventory' 
                    )

# Read anf format quadrat inventory file 
quad_inv        <- read_delim( quad_file, delim = delimiter, 
                               escape_double = FALSE, trim_ws = TRUE) %>%
                    as.data.frame %>% 
                    as.list

# Remove NAs
inv          <- lapply(X   = quad_inv, 
                       FUN = function(x) x[is.na(x) == FALSE])

# Replace dots in names with dashes
names(inv)   <- gsub( '\\.','-',names(inv) )  

# Read spatial data (polygon for each species per quadrat)
data <- {
  
    # Check for files ending with '_quadrats_filtered.rds' in the specified directory
    rds_files <- list.files( dat_dir, 
                             pattern = "_quadrats_filtered\\.rds$", 
                             full.names = TRUE)
    
    # If files matching the pattern are found, read the first one
    if (length(rds_files) > 0) {
      readRDS(file = rds_files[1])
    } else {
      # Fallback to the default file if no files matching the pattern are found
      readRDS(file = list.files(dat_dir, pattern = "all_filtered\\.rds$", full.names = TRUE)[1])
    }
    
  } %>%
  # Remove the 'type' column
  select(-any_of(c("type"))) %>%  
  # Rename columns to a standardized format
  setNames(quote_bare(Species, Site, Quad, Year, geometry))


# Subset data for the target species
dat_target_spec <- data %>% 
  subset( Species %in% target_spec$species ) %>% 
  # Exclude certain samples
  filter(
    if (exists("mod_plot") && length(mod_plot) > 0) !(Quad %in% mod_plot) else TRUE
  )

# set the buffer to 0.05 or 5 depending on the unit of measure used in GIS file
buff    <- if_else(st_bbox(dat_target_spec)[3] < 1.1, 0.05, 5)

# Forbs are not "clonal"; we assume non-forbs are clonal.
clonal  <- if_else(gr_form == 'forb', F, T)

# Prepare data for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec, inv,
                        # Dormancy flag
                        dorm         = 1,
                        # Buffer size
                        buff         = buff,
                        # Allow for clonal tracking
                        clonal       = clonal,
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
