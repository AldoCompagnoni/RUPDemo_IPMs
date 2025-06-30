library(tidyverse)
library(leaflet)
library(sf)

# Load up data from the species of interest (this will vary!)
source('moore_2021_az/R/popr/moore21_popr_ipm_mean.R')

# re-read in quadrat data
dir       <- 'moore_2021_az/data/Ancillary_Data_CSVs/'
quad_df   <- read.csv( paste0(dir,'Quadrat_Locations_and_Data.csv') )

# show how many SITES are represented
quad_df %>% 
  subset( Quadrat %in% unique(df$quad) ) %>% 
  select(Quadrat,Site)

# Do we have enough individuals? (at least 50 a year on average?)
count(df, year) 
  

# Where are the quadrats? Are the Sites very far from each other?
# NOTE: Now I play around with "Site". 
leaflet(quad_df) %>%
  addTiles() %>%
  addLabelOnlyMarkers(lng = ~Longitude_NAD_1983,
                      lat = ~Latitude_NAD_1983,
                      label = ~Site,
                      labelOptions = labelOptions(
                        noHide = TRUE,
                        direction = "top",
                        textOnly = TRUE,
                        style = list(
                          "color" = "red",
                          "font-family" = "serif",
                          "font-style" = "bold",
                          "box-shadow" = "3px 3px rgba(0,0,0,0.25)",
                          "font-size" = "14px",
                          "background-color" = "rgba(255,255,255,0.8)",
                          "padding" = "2px"
                        )
                      )
  )

# in this example, let's just keep "Big Fill" site
# count how many individuals are there

# what are the plot ids?
plot_ids <- quad_df %>% 
              subset( Site == 'Big Fill' ) %>% 
              .$Quadrat

# pretty damn good to keep the Big Fill data!
big_fill_counts <- df %>% 
                     subset( quad %in% plot_ids ) %>% 
                     count(year) %>% 
                     rename( big_fill_n = n )

# now, how many 
all_counts      <- count(df, year) %>% 
                     rename( all_n = n )

# how much data is lost? Increases through time, ~50% by end of series.
left_join( all_counts, big_fill_counts ) %>% 
  mutate( diff_n    = all_n - big_fill_n ) %>% 
  mutate( prop_lost = diff_n / all_n )
