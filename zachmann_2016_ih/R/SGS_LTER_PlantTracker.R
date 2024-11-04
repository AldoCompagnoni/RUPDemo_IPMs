# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# Aspen Workman, fall 2023
#
# Data sets were provided by the Shortgrass Steppe Long Term Ecological Research group,
# a partnership between Colorado State University, United States Department of Agriculture,
# Agricultural Research Service, and the U.S. Forest Service Pawnee National Grassland.
# Significant funding for these data was provided by the National Science Foundation Long
# Term Ecological Research program (NSF Grant Number DEB-1027319 and 0823405).
#
# Data source
# Chengjin Chu, John Norman, Robert Flynn, Nicole Kaplan, William K. Lauenroth,
# Peter B. Adler. 2013. Cover, density, and demographics of shortgrass steppe plants
# mapped 1997–2010 in permanent grazed and ungrazed quadrats. Ecology 94:1435.
# http://dx.doi.org/10.1890/13-0121.1
#
# Brief data description
# 14 years (1997-2010) of chart quadrat data from the Shortgrass Steppe Long Term Ecological
# Research project. 24 quadrats (1mx1m) were mapped annually. 6 pastures (19, 11, 24, 7, 5a, 5b)
# each contained 4 quadrats. Each pasture had two historically ungrazed and two historically
# grazed quadrats. One ungrazed quadrat was opened to grazing and one grazed quadrat was closed
# to grazing beginning in 1991. Each pasture contained 4 quadrats, ex. gzgz_5a, unun_5a,
# gzun_5a, and ungz_5a. Density-type species (forbs) were mapped as single points. Cover-type
# species (perennial grasses) were mapped as polygons (when mapped as points, converted to
# arbitrarily small square polygons with area of 0.25 cm sq). Some records were lost from 2000.
#
# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.
#
#
#
# Setup
#


library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)

base_dir <- ('zachmann_2016_ih')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")
shp_dir <- paste(base_dir, "/data/quadrat_data/shapefiles/msData/shapefiles/", sep="")

# setwd(dat_dir)

# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read_csv(paste0(dat_dir, "species_list.csv"))
# sp_name_changes <- read.csv("species_name_changes.csv") #will use to check names later on
grasses <- subset(sp_list, growthForm=="grass" & cover>100 & species!="Carex spp.")

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read_csv(paste0(dat_dir,"quad_inventory.csv")) %>% 
  select(-year)
quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_sgs <- quadInv_list


# Create a list of all shapefiles in the directory
shpFiles <- list.files(shp_dir)

quadYears <- unlist(strsplit(list.files(
  paste0(shp_dir,"/"),
  pattern = ".shp$"), split = ".shp"))

for (j in 1:length(quadYears)) {
  quadYearNow <- quadYears[j]
  shapeNow <- sf::st_read(dsn = paste0(shp_dir),
                          layer = quadYearNow)
  shapeNow$Site <- "Ih"
  shapeNow$Quad <- strsplit(quadYearNow, split = "_")[[1]][1]
  shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
  if (grepl(quadYearNow, pattern = "_D")) {
    shapeNow <- shapeNow[,!(names(shapeNow)
                            %in% c("OBJECTID", "seedling", "stem", "x", "y"))]
    shapeNow <- sf::st_buffer(x = shapeNow, dist = .0025)
    shapeNow$type <- "point"
  } else {
    shapeNow <- shapeNow[,!(names(shapeNow) %in% c("SP_ID", "stemID", "area", "x", "y"))]
    
    shapeNow$type <- "polygon"
  }
  if (j == 1) {
    dat <- shapeNow
  } else {
    dat <- rbind(dat, shapeNow)
  }
}


# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, paste0(dat_dir,"SGS_LTER_plantTracker_full.rds"))

# # Subset to the species of interest
# dat2 <- dat[dat$Species %in% grasses$species,]
# # And save the subsetted file, too
# saveRDS(dat2,file="SGS_LTER_plantTracker_grasses.rds")

# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  55, 343, 531, 691, 698, 718, 721, 738, 940, 956, 1123, 1130, 1161, 1330, 1332, 
  1688, 1724, 1927, 3486, 3612, 4404, 4765, 5339, 5724, 6072, 6370, 6374, 7196, 
  8772, 9675, 9733, 9789, 9792, 9800, 9914, 9915, 9927, 9958, 9988, 9998, 10000, 
  10033, 10043, 10175, 10187, 10219, 10233, 10238, 10240, 10430, 10614, 10618, 
  10653, 10952, 11180, 11430, 11494, 11495, 11522, 11530, 12023, 12200, 12218, 
  12265, 12553, 12565, 12571, 12577, 12603, 12607, 12610, 12629, 12886, 12908, 
  12930, 13126, 13153, 13158, 13203, 13420, 13582, 13667, 13675, 13789, 13955, 
  13980, 14233, 14259, 14276, 14386, 14396, 14398, 14406, 14422, 14895, 15992, 
  16192, 16571, 17199, 17594, 21621, 23199, 23212, 24037, 24039, 24049, 24060, 
  24066, 24203, 24205, 24208, 24251, 24342, 24347, 24399, 24403, 24412, 24451, 
  24528, 24793, 24796, 24912, 24916, 24934, 25198, 25245, 25263, 25359, 25403, 
  25533, 25559, 25589, 25682, 25897, 25930, 25941, 26136, 26145, 26162, 26543, 
  26863, 26891, 2690)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
 }
checkDat(dat3, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(25670,58091,75507,116003)
dat4 <- dat3[!(row.names(dat3) %in% drop_rows),]
checkDat(dat4, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat4,file="SGS_LTER_plantTracker_grasses_filtered.rds")

# Now the data are ready for the trackSpp function
datTrackSpp <- trackSpp(dat4,
                        inv_sgs,
                        dorm=1,
                        buff=0.05,
                        clonal=TRUE,
                        buffGenet = 0.05,
                        aggByGenet = TRUE,
                        flagSuspects = TRUE)


saveRDS(datTrackSpp,file="SGS_LTER_plantTracker_tracked.rds")

# Playing around with the data

all_spp <- datTrackSpp
for(i in 1:7552){
  all_spp$Treatment[i] <- str_split(all_spp$Quad[i],"_")[[1]][1]
  all_spp$Pasture[i] <- str_split(all_spp$Quad[i],"_")[[1]][2]
}

Ari_lon <- all_spp[all_spp$Species == "Aristida longiseta",] # 337 occurrences
Bou_gra <- all_spp[all_spp$Species == "Bouteloua gracilis",] # 3547 occurrences
Buc_dac <- all_spp[all_spp$Species == "Buchloe dactyloides",] # 999 occurrences
Sit_hys <- all_spp[all_spp$Species == "Sitanion hystrix",] # 471 occurrences
Spo_cry <- all_spp[all_spp$Species == "Sporobolus cryptandrus",] # 1256 occurrences
Sti_com <- all_spp[all_spp$Species == "Stipa comata",] # 880 occurrences


ggplot(Bou_gra) +
  geom_point(aes( x = log(basalArea_genet),y = survives_tplus1,group=Treatment, color=Treatment))

do.call( cbind, list( Bou_gra$Treatment,
                      Bou_gra$Year
                      ) ) %>% 
  as.data.frame %>% 
  count(V1,V2) %>% 
  arrange(V2,V1)

Bou_gra_df <- Bou_gra %>% st_drop_geometry()
Bou_gra_df$logsize <- log(Bou_gra_df$basalArea_genet)
write.csv(Bou_gra_df,"Bou_gra.csv")

Bou_gra_unun <- Bou_gra_df[Bou_gra_df$Treatment=="unun",]
Bou_gra_ungz <- Bou_gra_df[Bou_gra_df$Treatment=="ungz",]
Bou_gra_gzgz <- Bou_gra_df[Bou_gra_df$Treatment=="gzgz",]
Bou_gra_gzun <- Bou_gra_df[Bou_gra_df$Treatment=="gzun",]

Bou_gra_binned_unun <- plot_binned_prop(Bou_gra_unun,20,logsize,survives_tplus1)
Bou_gra_binned_ungz <- plot_binned_prop(Bou_gra_ungz,20,logsize,survives_tplus1)
Bou_gra_binned_gzgz <- plot_binned_prop(Bou_gra_gzgz,20,logsize,survives_tplus1)
Bou_gra_binned_gzun <- plot_binned_prop(Bou_gra_gzun,20,logsize,survives_tplus1)

Bou_gra_binned_trt <- data.frame(matrix(vector(), 80, 4,
                       dimnames=list(c(), c("binned_logsize", "binned_survival_t1", "treatment","n"))),
                stringsAsFactors=F)

Bou_gra_binned_trt$binned_logsize[1:20] <- Bou_gra_binned_gzgz$x_binned
Bou_gra_binned_trt$binned_survival_t1[1:20] <- Bou_gra_binned_gzgz$y_binned
Bou_gra_binned_trt$treatment[1:20] <- "gzgz"
Bou_gra_binned_trt$n[1:20] <- Bou_gra_binned_gzgz$n_s

Bou_gra_binned_trt$binned_logsize[21:40] <- Bou_gra_binned_gzun$x_binned
Bou_gra_binned_trt$binned_survival_t1[21:40] <- Bou_gra_binned_gzun$y_binned
Bou_gra_binned_trt$treatment[21:40] <- "gzun"
Bou_gra_binned_trt$n[21:40] <- Bou_gra_binned_gzun$n_s

Bou_gra_binned_trt$binned_logsize[41:60] <- Bou_gra_binned_unun$x_binned
Bou_gra_binned_trt$binned_survival_t1[41:60] <- Bou_gra_binned_unun$y_binned
Bou_gra_binned_trt$treatment[41:60] <- "unun"
Bou_gra_binned_trt$n[41:60] <- Bou_gra_binned_unun$n_s

Bou_gra_binned_trt$binned_logsize[61:80] <- Bou_gra_binned_ungz$x_binned
Bou_gra_binned_trt$binned_survival_t1[61:80] <- Bou_gra_binned_ungz$y_binned
Bou_gra_binned_trt$treatment[61:80] <- "ungz"
Bou_gra_binned_trt$n[61:80] <- Bou_gra_binned_ungz$n_s


ggplot(Bou_gra_binned_trt) +
  geom_point(aes(x=binned_logsize,y=binned_survival_t1,group=treatment,color=treatment))


Spo_cry_df <- Spo_cry %>% st_drop_geometry()
Spo_cry_df$logsize <- log(Spo_cry_df$basalArea_genet)
write.csv(Spo_cry_df,"Spo_cry.csv")
              

Buc_dac_df <- Buc_dac %>% st_drop_geometry()
Buc_dac_df$logsize <- log(Buc_dac_df$basalArea_genet)
write.csv(Buc_dac_df,"Buc_dac.csv")
