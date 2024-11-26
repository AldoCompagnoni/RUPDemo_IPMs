# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# Aspen Workman, fall 2023

# Data sets were provided by the Shortgrass Steppe Long Term Ecological Research group,
# a partnership between Colorado State University, United States Department of Agriculture,
# Agricultural Research Service, and the U.S. Forest Service Pawnee National Grassland.
# Significant funding for these data was provided by the National Science Foundation Long
# Term Ecological Research program (NSF Grant Number DEB-1027319 and 0823405).

# Data source
# Chengjin Chu, John Norman, Robert Flynn, Nicole Kaplan, William K. Lauenroth,
# Peter B. Adler. 2013. Cover, density, and demographics of shortgrass steppe plants
# mapped 1997–2010 in permanent grazed and ungrazed quadrats. Ecology 94:1435.
# http://dx.doi.org/10.1890/13-0121.1

# Brief data description
# 14 years (1997-2010) of chart quadrat data from the Shortgrass Steppe Long Term Ecological
# Research project. 24 quadrats (1mx1m) were mapped annually. 6 pastures (19, 11, 24, 7, 5a, 5b)
# each contained 4 quadrats. Each pasture had two historically ungrazed and two historically
# grazed quadrats. One ungrazed quadrat was opened to grazing and one grazed quadrat was closed
# to grazing beginning in 1991. Each pasture contained 4 quadrats, ex. gzgz_5a, unun_5a,
# gzun_5a, and ungz_5a. Density-type species (forbs) were mapped as single points. Cover-type
# species (perennial grasses) were mapped as polygons (when mapped as points, converted to
# arbitrarily small square polygons with area of 0.25 cm sq). Some records were lost from 2000.

# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.



# Setup
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)

rm(list = ls())

base_dir <- ('adler_2007_ks')
dat_dir <- paste(base_dir, "/data/quad_data/", sep="")
shp_dir <- paste(base_dir, "/data/quad_data/arcexport/", sep="")

delimiter <- ','

# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read_delim(paste0(dat_dir, "species_list.csv"),
                      delim = delimiter, escape_double = FALSE,
                      trim_ws = TRUE)

# grasses <- subset(sp_list, growthForm=="grass" & longevity=="P" & cover>100 & species!="Carex spp.")

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read_delim(paste0(dat_dir, "quadrat_inventory.csv"),
                       delim = delimiter, escape_double = FALSE,
                       trim_ws = TRUE)
quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_ks <- quadInv_list
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

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, file = paste0(dat_dir, "SGS_LTER_plantTracker_all.rds"))
dat <- readRDS(file = paste0(dat_dir, "SGS_LTER_plantTracker_all.rds"))
dat <- dat %>% mutate(Year = as.numeric(Year))


# # Subset to the species of interest
# dat2 <- dat[dat$Species %in% grasses$species,]
#
# # And save the subsetted file, too
# saveRDS(dat2, file = paste0(dat_dir, "SGS_LTER_plantTracker_grasses.rds"))

# Check the inv and dat arguments
checkDat(dat, inv_ks, species = "SCI_NAME", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat, file = paste0(dat_dir, "SGS_LTER_plantTracker_all_filtered.rds"))

# # Now the data are ready for the trackSpp function
# datTrackSpp <- trackSpp(dat02,
#                         inv_sgs,
#                         dorm=1,
#                         buff=0.05,
#                         clonal=TRUE,
#                         buffGenet = 0.05,
#                         aggByGenet = TRUE,
#                         flagSuspects = TRUE)
#
#
# saveRDS(datTrackSpp,file="SGS_LTER_plantTracker_tracked.rds")
#
# # Playing around with the data
#
# all_spp <- datTrackSpp
# for(i in 1:7552){
#   all_spp$Treatment[i] <- str_split(all_spp$Quad[i],"_")[[1]][1]
#   all_spp$Pasture[i] <- str_split(all_spp$Quad[i],"_")[[1]][2]
# }
#
# Ari_lon <- all_spp[all_spp$Species == "Aristida longiseta",] # 337 occurrences
# Bou_gra <- all_spp[all_spp$Species == "Bouteloua gracilis",] # 3547 occurrences
# Buc_dac <- all_spp[all_spp$Species == "Buchloe dactyloides",] # 999 occurrences
# Sit_hys <- all_spp[all_spp$Species == "Sitanion hystrix",] # 471 occurrences
# Spo_cry <- all_spp[all_spp$Species == "Sporobolus cryptandrus",] # 1256 occurrences
# Sti_com <- all_spp[all_spp$Species == "Stipa comata",] # 880 occurrences
#
#
# ggplot(Bou_gra) +
#   geom_point(aes( x = log(basalArea_genet),y = survives_tplus1,group=Treatment, color=Treatment))
#
# do.call( cbind, list( Bou_gra$Treatment,
#                       Bou_gra$Year
#                       ) ) %>%
#   as.data.frame %>%
#   count(V1,V2) %>%
#   arrange(V2,V1)
#
# Bou_gra_df <- Bou_gra %>% st_drop_geometry()
# Bou_gra_df$logsize <- log(Bou_gra_df$basalArea_genet)
# write.csv(Bou_gra_df,"Bou_gra.csv")
#
# Bou_gra_unun <- Bou_gra_df[Bou_gra_df$Treatment=="unun",]
# Bou_gra_ungz <- Bou_gra_df[Bou_gra_df$Treatment=="ungz",]
# Bou_gra_gzgz <- Bou_gra_df[Bou_gra_df$Treatment=="gzgz",]
# Bou_gra_gzun <- Bou_gra_df[Bou_gra_df$Treatment=="gzun",]
#
# Bou_gra_binned_unun <- plot_binned_prop(Bou_gra_unun,20,logsize,survives_tplus1)
# Bou_gra_binned_ungz <- plot_binned_prop(Bou_gra_ungz,20,logsize,survives_tplus1)
# Bou_gra_binned_gzgz <- plot_binned_prop(Bou_gra_gzgz,20,logsize,survives_tplus1)
# Bou_gra_binned_gzun <- plot_binned_prop(Bou_gra_gzun,20,logsize,survives_tplus1)
#
# Bou_gra_binned_trt <- data.frame(matrix(vector(), 80, 4,
#                        dimnames=list(c(), c("binned_logsize", "binned_survival_t1", "treatment","n"))),
#                 stringsAsFactors=F)
#
# Bou_gra_binned_trt$binned_logsize[1:20] <- Bou_gra_binned_gzgz$x_binned
# Bou_gra_binned_trt$binned_survival_t1[1:20] <- Bou_gra_binned_gzgz$y_binned
# Bou_gra_binned_trt$treatment[1:20] <- "gzgz"
# Bou_gra_binned_trt$n[1:20] <- Bou_gra_binned_gzgz$n_s
#
# Bou_gra_binned_trt$binned_logsize[21:40] <- Bou_gra_binned_gzun$x_binned
# Bou_gra_binned_trt$binned_survival_t1[21:40] <- Bou_gra_binned_gzun$y_binned
# Bou_gra_binned_trt$treatment[21:40] <- "gzun"
# Bou_gra_binned_trt$n[21:40] <- Bou_gra_binned_gzun$n_s
#
# Bou_gra_binned_trt$binned_logsize[41:60] <- Bou_gra_binned_unun$x_binned
# Bou_gra_binned_trt$binned_survival_t1[41:60] <- Bou_gra_binned_unun$y_binned
# Bou_gra_binned_trt$treatment[41:60] <- "unun"
# Bou_gra_binned_trt$n[41:60] <- Bou_gra_binned_unun$n_s
#
# Bou_gra_binned_trt$binned_logsize[61:80] <- Bou_gra_binned_ungz$x_binned
# Bou_gra_binned_trt$binned_survival_t1[61:80] <- Bou_gra_binned_ungz$y_binned
# Bou_gra_binned_trt$treatment[61:80] <- "ungz"
# Bou_gra_binned_trt$n[61:80] <- Bou_gra_binned_ungz$n_s
#
#
# ggplot(Bou_gra_binned_trt) +
#   geom_point(aes(x=binned_logsize,y=binned_survival_t1,group=treatment,color=treatment))
#
#
# Spo_cry_df <- Spo_cry %>% st_drop_geometry()
# Spo_cry_df$logsize <- log(Spo_cry_df$basalArea_genet)
# write.csv(Spo_cry_df,"Spo_cry.csv")
#
#
# Buc_dac_df <- Buc_dac %>% st_drop_geometry()
# Buc_dac_df$logsize <- log(Buc_dac_df$basalArea_genet)
# write.csv(Buc_dac_df,"Buc_dac.csv")
