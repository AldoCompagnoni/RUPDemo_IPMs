# Second part of the plant tracker pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.29

# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)


# Specifications ----------------------------------------------------------------
# Define the species variable and abbreviation
v_species <- target_spec[1,1]
v_sp_abb  <- v_species %>% 
              strsplit(' ') %>% 
              pluck(1) %>% 
              substr( 1, 2 ) %>% 
              paste(collapse = '') %>% 
              str_replace_all(' ', '') %>% 
              tolower


# Directory --------------------------------------------------------------------
if ((tolower(v_gr_form) == 'forb' | 
     tolower(v_gr_form) == 'shrub' | 
     tolower(v_gr_form) == 'density')) {
  v_folder_suffix <- paste0(v_sp_abb, '_mpm')
  dir_R      <- file.path(dir_pub, 'R', v_folder_suffix)
  dir_data   <- file.path(dir_pub, 'data', v_folder_suffix)
  dir_result <- file.path(dir_pub, 'results', v_folder_suffix)
} else {
  dir_R      <- file.path(dir_pub, 'R', v_sp_abb)
  dir_data   <- file.path(dir_pub, 'data', v_sp_abb)
  dir_result <- file.path(dir_pub, 'results', v_sp_abb)
}


# Data -------------------------------------------------------------------------
# Read anf format quadrat inventory file 
# Identify quadrat inventory file 
quad_file <- list.files(
  path = dir_qud, full.names = TRUE,
  pattern = 'quad_inventory|quadrat_inventory'
)[1]

if (grepl("\\.RData$", quad_file)) {
  quad_inv <- readRDS(quad_file)
} else {
  quad_inv <- read_delim(
    quad_file, delim = v_delimiter, escape_double = FALSE, trim_ws = TRUE) %>%
    as.data.frame %>% 
    as.list()
}

# Remove NAs
inv        <- lapply(
  X = quad_inv, FUN = function(x) x[is.na(x) == FALSE])

# Replace dots in names with dashes
names(inv) <- gsub('\\.','-',names(inv))  

# Read spatial data (polygon for each species per quadrat)
data <- {
  # Check for files ending with '_quadrats_filtered.rds' in the specified directory
  rds_files <- list.files( 
    dir_qud, pattern = '_quadrats_filtered\\.rds$', full.names = TRUE)
  # If files matching the pattern are found, read the first one
  if (length(rds_files) > 0) {
    readRDS(file = rds_files[1])
  } else {
    # If not fallback to the default file
    readRDS(file = list.files(
      dir_qud, pattern = 'all_filtered\\.rds$', full.names = TRUE)[1])
  }
  
} %>%
  # Remove the 'type' column
  select(-any_of(c('type'))) %>% 
  clean_names() #%>%   
# This is to accommodate older versions of 'quadrat' scripts
#set_names(c('species', 'site', 'quad', 'year', 'geometry'))


# Subset data for the target species
dat_target_spec <- data %>%
  subset(species %in% target_spec$species) %>%
  filter(
    if (exists("mod_plot") && length(mod_plot) > 0) !(quad %in% mod_plot)
    else TRUE
  ) %>%
  select(
    species, site, 
    matches("quad"),  # selects quad or quadrat
    year, geometry,
    everything()      # all remaining columns
  ) %>% 
  rename(Species = species,
         Site    = site,
         Quad    = matches("quad"),
         Year    = year)


# set the buffer to 0.05 or 5 depending on the unit of measure used in GIS file
buff <- v_buff <- if_else(
  st_bbox(dat_target_spec)[3] < 1.1 | st_bbox(dat_target_spec)[3] > 100, 
  0.05, 5)

# Forbs are not "clonal"; we assume non-forbs are clonal.
v_clonal  <- if_else(v_gr_form == 'forb', F, T)

# Prepare data for the trackSpp function
datTrackSpp <- trackSpp(
  dat_target_spec,
  inv,
  # Dormancy flag
  dorm         = 1,
  # Buffer size
  buff         = v_buff,
  # Allow for clonal tracking
  clonal       = v_clonal,
  # Buffer for genet
  buffGenet    = v_buff,
  # Aggregate by genet
  aggByGenet   = TRUE,
  # Flag potential issues
  flagSuspects = TRUE)   

# Create output directory if it doesn't exist
if (!dir.exists(dir_data)) {
  dir.create(dir_data)
}

# Save the tracked data to a CSV file
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select(
    Site, Quad, Species, trackID, Year, basalArea_genet, recruit, 
    survives_tplus1, age, size_tplus1, nearEdge, Suspect) %>%
  clean_names() %>% 
  write.csv(file.path(
    dir_data, paste0(
      v_script_prefix, '_', v_sp_abb, '.csv')), row.names = FALSE)  


