# Combine species binomial in a single column 
library(tidyverse)

# read in species list
spp_l <- read.csv("anderson_2016_mt/data/quadrat_data/species_list.csv")

# check whether the species list has already been changed.
#   in that case, column "X" will not be present
if( "X" %in% names(spp_l) ){
  
  # create species column containing entire species binomial
  spp_l %>%
    mutate(species = paste(species, X, sep = " ") ) %>% 
    # remove column X (now contained in "species")
    dplyr::select(-X) %>% 
    # replace original file
    write.csv("anderson_2016_mt/data/quadrat_data/species_list.csv", 
              row.names = F)  
}

