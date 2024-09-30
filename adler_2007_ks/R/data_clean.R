# IPM data clean for alder 2007 & chu 2016; kansas; bocu

# Niklas Neisse
# 2024.09.30

# reading in, explore, and cleaning the data
# setting up the vital rate data-frames for the year specific 

# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
rm(list = ls()) 
# Set seed for reproducibility
set.seed(100) 
rm( list = ls() )
options( stringsAsFactors = F )
# Set working directory
setwd("C:/code/RUPDemo_IPMs")

# Packages ---------------------------------------------------------------------

# Define CRAN packages
.cran_packages <- c('tidyverse','patchwork','skimr') 
# Check if CRAN packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}
# Load required packages
sapply(.cran_packages, require, character.only = TRUE) 


# Data -------------------------------------------------------------------------
## 

df <- read_csv("adler_2007_ks/data/KS_grasses_all.csv") %>% 
  filter(Site == 'KS')

write.csv(df, "adler_2007_ks/data/ks_grasses.csv", row.names = F)


