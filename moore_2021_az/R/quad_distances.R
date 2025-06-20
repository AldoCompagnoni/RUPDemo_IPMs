# calculate quatrat distances
library(sf)
library(tidyverse)

# re-read in quadrat data
dir       <- 'moore_2021_az/data/Ancillary_Data_CSVs/'
quad_df   <- read.csv( paste0(dir,'Quadrat_Locations_and_Data.csv') )

# vector of "unique" distances
dist_v    <- quad_df %>% 
              select( Latitude_NAD_1983, Longitude_NAD_1983 ) %>% 
              # mutate( id = 1:nrow(.) ) %>% 
              rename( lat = Latitude_NAD_1983, 
                      lon = Longitude_NAD_1983 ) %>% 
              st_as_sf( coords = c("lon", "lat"), crs = 4326 ) %>% 
              st_distance %>% 
              as.numeric %>% 
              unique  


# ok wow, 20Km distance is relatively common...
hist(dist_v / 1000)
mean(dist_v / 1000)
  