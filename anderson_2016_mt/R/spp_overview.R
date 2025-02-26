# Counting number of individuals per species 
# This helps choosing which species to focus on

# Load relevant packages
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse) #ver 1.1.0

# set up working directories
dir     <- 'anderson_2016_mt/data/quadrat_data/'
dat_dir <- dir

# Read in the data
dat <- readRDS(file=paste0(dat_dir,"anderson16mt_plantTracker_all_filtered.rds"))

# species list from the original data
spp_list  <- read.csv( paste0(dat_dir, 'species_list.csv') ) %>% 
               mutate( species = paste0(species, ' ', X) ) %>% 
               dplyr::select(-X) 

# count the number of feature-per-species
spp_count <- dat %>% 
              count(species) %>% 
              as.data.frame %>% 
              dplyr::select(species,n) %>% 
              left_join( spp_list ) %>% 
              arrange(desc(n) )

# Check which species are most abundant
spp_count %>% View