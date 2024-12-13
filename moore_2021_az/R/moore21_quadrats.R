# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs

# Outline
# This script uses the plantTracker package (Stears et al. 2022, Methods in Ecology & Evolution)
# to convert this chart quadrat data to demographic data that will be used to parameterize vital
# rate models for the construction of integral projection models (IPMs) of the perennial grasses.

# Publication: https://doi.org/10.1002/ecy.3661

library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)

base_dir <- ('moore_2021_az')
dat_dir <- paste(base_dir, "/data/Ancillary_Data_CSVs", sep="")
shp_dir <- paste(dat_dir, 
                 "", sep="")

# setwd(dat_dir)

# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read_csv(paste0(dat_dir, "/Plant_Species_List.csv"))

# Read in quad inventory:
#  a table that has all the quadrats as columns and the respective year 
#   they've been sampled as rows with NA for years that are missing
#  to use as 'inv' list in plantTracker
quad_inv <- read_csv(paste0(dat_dir,"/Summarize_Quadrats_by_Year.csv")) %>% 
  select(-c('Years_Surveyed', 'Proportion', 'Comments')) %>%
  column_to_rownames("Quadrat") %>%   # Move the first column to row names
  t() %>%                             # Transpose the dataframe
  as.data.frame() %>%                 # Convert back to dataframe
  rownames_to_column("Year") %>%
  mutate(Year = substr(as.character(Year), 3, 4))

quad_inv[,-1] <- apply(quad_inv[,-1], 2, function(x) ifelse(x == "X", quad_inv$Year, x))


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
  shapeNow$Site <- "Id"
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
dat <- readRDS(paste0(dat_dir,"SGS_LTER_plantTracker_full.rds"))

# # Subset to the species of interest
# dat2 <- dat[dat$Species %in% grasses$species,]
# # And save the subsetted file, too
# saveRDS(dat2,file="SGS_LTER_plantTracker_grasses.rds")

# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  255, 343, 531, 691, 698, 718, 721, 738, 940, 956, 1123, 1130, 1161, 1330, 1332, 
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
  26863, 26891)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
 }

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom2 <- c(
  26902, 26911, 26929, 27207, 27214, 27223, 27376, 27405, 27438, 27605, 29292, 
  32885, 36773, 36810, 36854, 36865, 36879, 36942, 36982, 36991, 37003, 37006, 
  37120, 37227, 37818, 38062, 38069, 38071, 38076, 38082, 38092, 38106, 38350, 
  38359, 38366, 38378, 38570, 38619, 38857, 38868, 38883, 38886, 39150, 39153, 
  39160, 39165, 39176, 39179, 39500, 39993, 39996, 40032, 40257, 40259, 40406, 
  40418, 40430, 40618, 40934, 40943, 40965, 41268, 41278, 41491, 41523, 41572, 
  41661, 41662, 41695, 42100, 42117, 42149, 42177, 42248, 42295, 42448, 42459, 
  42495, 42511, 47870, 48337, 49282, 49285, 49333, 49552, 50796, 51316, 57159, 
  57160, 57184, 57186, 57190, 57193, 57229, 57328, 57332, 57334, 57390, 57558, 
  57562, 57609, 57611, 57841, 57842, 57843, 57844, 58081, 58087, 58113, 58125, 
  58426, 58484, 58485, 58488, 58824, 58825, 58826, 59131, 59407, 59411, 59814, 
  60338, 60339, 60340, 60936, 60966, 60995, 60997, 61457, 61488, 61509, 61546, 
  61903, 62024, 62048)
for(i in 1:length(invalid_geom2)){
  dat01[invalid_geom2[i],6] <- st_make_valid(dat01[invalid_geom2[i],6])
}

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom3 <- c(
  62053, 62055, 62337, 62394, 62411, 62442, 62444, 62446, 62680, 62718, 62723, 
  62791, 62792, 62793, 62794, 62796, 63095, 63099, 63119, 63120, 63121, 63122, 
  63402, 63500, 63502, 63867, 65053, 66083, 69089, 70091, 71068, 72536, 72586, 
  72635, 72895, 72935, 72985, 72989, 73239, 73245, 73264, 73270, 73280, 73508, 
  73531, 80562, 80617, 81052, 81061, 81100, 81101, 81103, 81619, 81670, 81963, 
  82079, 82101, 82158, 82176, 82361, 82430, 82648, 82654, 82658, 82865, 82876, 
  82888, 82902, 83250, 83262, 83267, 83525, 83564, 83571, 83839, 83872, 84151, 
  84175, 84462, 84468, 84483, 84487, 84507, 84512, 84525, 84526, 84528, 84529, 
  84530, 84532, 84534, 84784, 84805, 84818, 84820, 84830, 84886, 84889, 84891, 
  85263, 85292, 85294, 85298, 85302, 85323, 85326, 85331, 85332, 85333, 86129, 
  86135, 86136, 86155, 86184, 86519, 87149, 87206, 87399, 87724, 88443, 88483, 
  88484, 88491, 88658, 88674, 88693, 89244, 90461, 90768, 91250, 91831, 91917, 
  92426, 92428, 93055)
for(i in 1:length(invalid_geom3)){
  dat01[invalid_geom3[i],6] <- st_make_valid(dat01[invalid_geom3[i],6])
}

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom4 <- c(
  93090, 93100, 93520, 93820, 93847, 94129, 94462, 94755, 94979, 95287, 95333, 
  97101, 98687, 101685, 101757, 101864, 101992, 102231, 102258, 103113, 103119, 
  103124, 103143, 103479, 103480, 103511, 103517, 103822, 103870, 103875, 
  103882, 104161, 104165, 104182, 104187, 104216, 104406, 104408, 104423, 
  104424, 104630, 104658, 104687, 104704, 104966, 105219, 105501, 105550, 
  105810, 106033, 106040, 106231, 106243, 106266, 106268, 106270, 106498, 
  106925, 107994, 108122, 108437, 109119, 109120, 109321, 109356, 109357, 
  109664, 110134, 111643, 111731, 112878, 113623, 114117, 115306, 117155, 
  117162, 117163, 117685, 120601)
for(i in 1:length(invalid_geom4)){
  dat01[invalid_geom4[i],6] <- st_make_valid(dat01[invalid_geom4[i],6])
}

checkDat(dat01, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

# # Still have a couple of repeated rows, somehow, so we will drop those
# drop_rows <- c(25670,58091,75507,116003)
# dat4 <- dat3[!(row.names(dat3) %in% drop_rows),]
# checkDat(dat4, inv_sgs, species = "species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")


dat01 <- dat01 %>% 
  rename(Species = species)

saveRDS(dat01, paste0(dat_dir,"SGS_LTER_plantTracker_all_filtered.rds"))

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
