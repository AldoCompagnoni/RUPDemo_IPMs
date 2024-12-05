# run all plant tracker files, all at once
library(tidyverse)

# Function to quote bare names for tidy evaluation
quote_bare <- function(...){
  substitute(alist(...)) %>% 
    eval() %>% 
    sapply(deparse)  
}

# Set up working directory and species codes
dir         <- 'adler_2007_ks/R/'
spp_codes_v <- quote_bare( ange, arlo, bocu, bogr, 
                           buda, pavi, scsc, spas, spcr )


# function to run mean models (models with all data comprised)
run_mean_models <- function( spp_code_x ){
  print( spp_code_x )
  source( paste0(dir, spp_code_x, '/',
                 'ipm_mean_',spp_code_x,'.R') )
}

# Run plant tracker for all species
lapply( spp_codes_v, run_mean_models )
