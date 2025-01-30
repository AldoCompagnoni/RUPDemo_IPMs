# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs

# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)


# Specifications ---------------------------------------------------------------
# Directories
base_dir <- ('anderson_2016_mt')
dat_dir <- paste(base_dir, "/data/quadrat_data/", sep="")
shp_dir <- paste(base_dir, 
                 "/data/quadrat_data/shapefiles/", sep="")


# Data -------------------------------------------------------------------------
# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read.delim(sep = ';', paste0(dat_dir, "species_list.csv"))

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read.csv(paste0(dat_dir,"quad_inventory.csv"))
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
  
  # Read the shapefile
  shapeNow <- sf::st_read(dsn = paste0(shp_dir), 
                          layer = quadYearNow, 
                          # Because of some custom corrdi ref system,
                          #  might aswell use: st_crs(4326)
                          crs = NA)
  
  # Add Site, Quad, and Year
  shapeNow$Site <- "Mt"
  shapeNow$Quad <- strsplit(quadYearNow, split = "_")[[1]][1]
  shapeNow$Year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][2])
  
  # Handle specific columns based on layer type (_D or not)
  if (grepl(quadYearNow, pattern = "_D")) {
    # Keep only the relevant columns for point data
    shapeNow <- shapeNow %>%
      select(Species, Site, Quad, Year, geometry) %>%
      mutate(type = "point")
    
    # Apply buffer if necessary for point data
    shapeNow <- st_buffer(shapeNow, dist = .0025)
  } else {
    # Keep only the relevant columns for polygon data
    shapeNow <- shapeNow %>%
      select(Species, Site, Quad, Year, geometry) %>%
      mutate(type = "polygon")
  }
  
  # Combine all data into the `dat` variable
  if (j == 1) {
    dat <- shapeNow
  } else {
    dat <- rbind(dat, shapeNow)
  }
} 

dat <- dat %>% rename(species = Species)

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, paste0(dat_dir,"anderson16mt_SGS_LTER_plantTracker_all.rds"))
dat <- readRDS(paste0(dat_dir,"anderson16mt_SGS_LTER_plantTracker_all.rds"))

# # Subset to the species of interest
# dat2 <- dat[dat$Species %in% grasses$species,]
# # And save the subsetted file, too
# saveRDS(dat2,file="SGS_LTER_plantTracker_grasses.rds")

# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  452, 481, 686, 713, 2224, 3287, 3736, 18640, 19211, 36792, 42830, 44822, 
  45384, 46597, 46604, 46611, 47825, 49400, 51516, 51518, 58996, 59722, 59999, 
  71625, 76288, 76471, 76694, 76773, 82224, 82507, 82520, 92972, 100890, 103095, 
  103189, 103396, 107938, 107946, 108262, 108275, 108669, 109006, 109482, 
  110340, 110519, 110622, 110643, 110948, 111336, 113642, 114459, 115375, 
  115444, 115883, 115959, 116622, 120225, 120230, 121366, 121895, 122632, 
  124332, 124499, 125595, 126685, 130351, 130352, 130691, 134742, 144658, 
  144661, 144679, 148967, 149133, 150049, 151646, 153203, 154089, 154332, 
  157148, 162078, 163308, 164625, 171411, 175206, 175208, 177238, 177456, 
  180368, 182126, 182928, 185677, 185702, 186382, 188348, 189072, 189088, 
  189117, 189426, 189448, 189473, 189864, 190134, 190152, 190296, 190474, 
  203736, 204777)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
 }

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  20522, 26373, 34704, 67889, 74969, 132988, 135544, 165470, 167341)
dat02 <- dat01[!(row.names(dat01) %in% drop_rows),]
checkDat(dat02, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")



saveRDS(dat02, paste0(dat_dir,"anderson16mt_SGS_LTER_plantTracker_all_filtered.rds"))

# # Now the data are ready for the trackSpp function
# datTrackSpp <- trackSpp(dat4,
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
