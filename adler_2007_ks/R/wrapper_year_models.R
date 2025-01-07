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
                           bohi, scsc, spcr )

# function to run year-specific IPMs
run_yearly_ipms <- function( spp_code_x ){
  print( spp_code_x )
  source( paste0(dir, spp_code_x, '/',
                 'adler07_',spp_code_x,'_ipm_year_specific.R') )
}

# Run plant tracker for all species
lapply( spp_codes_v, run_yearly_ipms )
