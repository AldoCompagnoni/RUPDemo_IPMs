library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse) 

dat_dir <- 'adler_2007_ks/data/quadrat_data/'
shp_dir <- paste0(dat_dir, "arcexport/")

# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list         <- read.csv( paste0( dat_dir, "species_list.csv") )
grasses         <- sp_list %>% 
                    dplyr::arrange( desc(count) ) %>% 
                    .[1:3,]
# grasses         <- subset(sp_list, growthForm=="grass" & longevity=="P" & cover>100 & species!="Carex spp.")

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv        <- read.csv( paste0(dat_dir, "quadrat_inventory.csv"),
                              sep=',') %>% 
                    dplyr::select(-year)
quadInv_list    <- as.list(quad_inv)
quadInv_list    <- lapply(X = quadInv_list, 
                          FUN = function(x) x[is.na(x) == FALSE])
inv_ks          <- quadInv_list 
names(inv_ks)   <- gsub( '\\.','-',names(inv_ks) )

# Read in shapefiles to create sf dataframe to use as 'dat' in plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)
# Create list of quad names
quadNames <- list.files(shp_dir)

# Use for loop to download data from each quad folder
for(i in 1:length(quadNames)){
  
  quadNow      <- quadNames[i]
  # quad/year name combination
  quadYears    <- paste0(shp_dir,quadNow,"/") %>% 
                    list.files( pattern = ".e00$" ) 
  
  # loop over each year  
  for(j in 1:length(quadYears)){
    
    quad_yr_name  <- quadYears[j]
    shapeNow      <- sf::st_read( dsn   = paste0(shp_dir,quadNow,"/",
                                                 quad_yr_name ),
                                  layer = 'PAL' ) %>% 
                                  dplyr::select( SCI_NAME, geometry)
    shapeNow$Site <- "KS"
    shapeNow$Quad <- quadNow
    shapeNow$Year <- quad_yr_name %>% 
                      gsub(quadNow,'',.) %>% 
                      gsub('.e00','',.) %>%
                      unlist
    
    # start final data frame
    if(i == 1 & j == 1) {
      
      dat <- shapeNow
    
    # "append" data to the initial data frame  
    } else {
      
      dat <- rbind(dat, shapeNow)
      
    }
  }
}

# 
dat$Year <- as.numeric(dat$Year)

# Save the output file so that it doesn't need to be recreated ever again
# saveRDS(dat,file="KS_polygons_full.rds")
dat <- readRDS(file="KS_polygons_full.rds")

