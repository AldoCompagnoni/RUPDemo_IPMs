# First part of the plant tracker pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.24

# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)

# Packages ---------------------------------------------------------------------
# Load packages, verify, and download if needed
source('helper_functions/load_packages.R')
load_packages(sf, plantTracker, tidyverse)

# Data -------------------------------------------------------------------------
# Directory 1
pub_dir <- file.path(paste0(author_year, '_', region_abb))
dat_dir <- file.path(pub_dir, 'data', '//quadrat_data/', sep='')
# Prefix for the script name
script_prefix <- str_c(str_extract(author_year, '^[^_]+'), 
                       str_sub(str_extract(author_year, '_\\d+$'), -2, -1))


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
sp_list <- read_csv(paste0(dat_dir, '/species_list.csv')) %>% 
  as.data.frame() %>% 
  filter(growthForm == gr_form) %>% 
  {
    if(gr_form == 'forb') {
      filter(., is.na(cover)) %>%
        arrange(desc(density))
    } else {
      arrange(., desc(cover))
    }
  }