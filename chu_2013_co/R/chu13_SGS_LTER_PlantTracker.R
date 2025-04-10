# Using plantTracker to convert chart quadrat data from SGS LTER to demographic data for IPMs
# Aspen Workman, fall 2023 - Modified by Niklas Neiße, spring 2025
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


# Packages ---------------------------------------------------------------------
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse)
library(janitor)


# Directories ------------------------------------------------------------------
dir_pub <- file.path('chu_2013_co')
dir_dat <- file.path(dir_pub, 'data')
dir_qud <- file.path(dir_dat, 'quadrat_data')
dir_shp <- file.path(dir_qud, 'shapefiles')


# Data -------------------------------------------------------------------------
# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list <- read.delim(file.path(dir_qud, "species_list.csv"))
sp_name_changes <- read.csv(file.path(dir_qud, "species_name_changes.csv")) 



# #  will use to check names later on
# grasses <- subset(sp_list, growthForm=="grass" & longevity=="P" & cover>100 & species!="Carex spp.")



# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv <- read.delim(file.path(dir_qud, "quad_inventory.csv"))
quadInv_list <- as.list(quad_inv)
quadInv_list <- lapply(X = quadInv_list, FUN = function(x) x[is.na(x) == FALSE])
inv_sgs <- quadInv_list

# Read in shapefiles to create sf dataframe to use as 'dat' in plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)
# Create list of quad names
quadNames <- list.files(file.path(dir_qud, 'shapefiles'))
# Use for loop to download data from each quad folder
for(i in 1:length(quadNames)){
  quadNow <- quadNames[i]
  quadYears <- unlist(strsplit(list.files(
    file.path(dir_shp, quadNow),
    pattern = ".shp$"), split = ".shp"))
  for (j in 1:length(quadYears)) {
    quadYearNow <- quadYears[j]
    shapeNow <- sf::st_read(dsn = paste0(dir_shp,quadNow),
                            layer = quadYearNow)
    shapeNow$site <- "co"
    shapeNow$quad <- quadNow
    shapeNow$year <- as.numeric(strsplit(quadYearNow, split = "_")[[1]][4])
    if (grepl(quadYearNow, pattern = "pnt")) {
      shapeNow <- shapeNow[,!(names(shapeNow)
                              %in% c("coords_x1", "coords_x2", "coords_x1_", "coords_x2_", "coords_x1.1", "coords_x2.1"))]
      shapeNow <- sf::st_buffer(x = shapeNow, dist = .0025)
      shapeNow$type <- "point"
    } else {
      shapeNow <- shapeNow[,!(names(shapeNow) %in% c("SP_ID", "SP_ID_1", "area", "x", "y"))]
      
      shapeNow$type <- "polygon"
    }
    if (i == 1 & j == 1) {
      dat <- shapeNow
    } else {
      dat <- rbind(dat, shapeNow)
    }
  }
}

# Save the output file so that it doesn't need to be recreated ever again
saveRDS(dat, file = paste0(dir_qud, "SGS_LTER_plantTracker_full.rds"))
dat <- readRDS(file = paste0(dir_qud, "SGS_LTER_plantTracker_full.rds"))

# # Subset to the species of interest
# dat2 <- dat[dat$Species %in% grasses$species,]
# 
# # And save the subsetted file, too
# saveRDS(dat2, file = paste0(dir_qud, "SGS_LTER_plantTracker_grasses.rds"))

# Check the inv and dat arguments
checkDat(dat, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
# Some rows had invalid geometry, so we fix the geometries
invalid_geom <- c(
  2311, 2347, 2649, 2929, 3229, 3301, 3874, 3877, 3879, 3933, 4008, 4012, 4163, 
  4248, 4275, 4305, 4480, 4509, 4890, 5020, 6461, 6858, 6892, 7025, 7094, 7118, 
  7322, 7634, 7684, 8187, 8315, 8397, 8415, 8521, 8537, 8742, 9350, 9513, 10506, 
  10836, 11043, 11117, 11121, 11153, 11464, 11487, 11822, 11849, 11906, 11926, 
  12160, 12200, 12254, 12704, 13489, 14743, 15185, 15647, 15650, 15892, 15978, 
  16376, 16499, 16906, 17117, 17583, 18024, 18347, 18374, 18406, 18427, 18778, 
  18786, 19427, 20258, 20665, 21173, 21340, 21818, 22169, 22202, 22372, 23128, 
  23155, 23478, 23950, 24536, 25893, 25935, 27171, 27523, 27775, 28132, 28525, 
  28529, 29028, 29045, 29052, 29071, 29102, 29105, 29205, 29539, 29943, 30065, 
  30219, 30635, 31042, 31617, 31634, 34811, 35096, 35388, 35683, 35948, 36047, 
  36139, 36174, 36194, 36211, 36219, 36220, 36552, 36582, 36738, 36897, 37347, 
  37385, 37388, 37471, 39161, 39549, 39733, 39927, 40093, 40315, 40324, 40459, 
  40783, 40834, 41218, 41226, 41252, 41258, 41324, 41330, 41502, 41602, 41645, 
  42148, 42336, 43659, 44177, 44285, 44376, 44593, 44756, 44778, 44847, 44892, 
  44894, 44962, 45351, 45515, 45559, 45652, 45682, 47282, 47496, 47577, 47904, 
  48020, 48378, 48764, 49028, 49057, 49097, 49159, 49341, 49370, 49440, 49603, 
  50608, 50703, 51166, 51175, 51377, 51402, 51557, 51683, 51849, 52111, 52454, 
  52457, 53041, 53448, 53449, 53471, 53523, 54208, 55464, 55764, 55872, 55905, 
  56139, 56314, 56496, 56826, 57085, 57086, 57156, 57183, 57508, 57766, 57850, 
  57975, 58194, 58795, 59964, 60244, 60261, 60906, 61559, 61957, 61962, 62564, 
  63315, 63419, 63736, 63739, 63764, 63808, 63809, 64229, 64320, 65145, 66132, 
  68407, 68586, 68878, 68922, 69033, 69300, 69315, 69401, 69827, 69983, 70051, 
  70375, 70435, 70468, 70470, 70493, 70503, 70552, 70663, 70912, 71118, 74592, 
  75064, 75067, 75433, 75800, 75962, 76001, 76473, 77872, 78722, 79486, 79550, 
  80555, 82966, 83268, 83403, 83447, 83746)
dat01 <- dat
for(i in 1:length(invalid_geom)){
  dat01[invalid_geom[i],6] <- st_make_valid(dat01[invalid_geom[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom3 <- c(
  84097, 84276, 84278, 84284, 84287, 84788, 84790, 84797, 84891, 84983, 85272, 
  85308, 85432, 85722, 85915, 86004, 86006, 86136, 86855, 87068, 87679, 87700, 
  87851, 87980, 88221, 88296, 88600, 88704, 88705, 89198, 89209, 89250, 89325, 
  89828, 90245, 90560, 90908, 91531, 93498, 93810, 94125, 94231, 94421, 94516, 
  94537, 94774, 94957, 95176, 95177, 95590, 95593, 95625, 95659, 95681, 95705, 
  95722, 95948, 96222, 96339, 96876, 99074, 99278, 99528, 100124, 100136, 
  100774, 100807, 100833, 101630, 101737, 102260, 102600, 102605, 102630, 
  107938, 108186, 108274, 108644, 109236, 109762, 110054, 110428, 110810, 
  110939, 111117, 111211, 111543, 111666, 111770, 113269, 113853, 114078, 
  114299, 114300, 114343, 114691, 115359, 115535, 115554, 115560, 115689, 
  115904, 116045, 116069, 116103, 116217, 116304, 118893, 119007, 119208, 
  119220, 119562, 119792, 120299, 120302, 120726, 120873, 121031, 121034, 
  121039, 121411, 121461, 121544, 121691, 121805, 121887
)
for(i in 1:length(invalid_geom3)){
  dat01[invalid_geom3[i],6] <- st_make_valid(dat01[invalid_geom3[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
invalid_geom4 <- c(
  123325, 123790, 123984, 124009, 124070, 124334, 124681, 125120, 125129, 
  125839, 126288, 130048, 130659, 130935, 130961, 131130, 131281, 131344, 
  131349, 131406, 131475, 131606, 131674, 131676, 131751, 131762, 131841, 
  132125, 132232, 132253, 132258, 132548, 132609, 132623, 133163, 133274
)

for(i in 1:length(invalid_geom4)){
  dat01[invalid_geom4[i],6] <- st_make_valid(dat01[invalid_geom4[i],6])
}

checkDat(dat01, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")
dat02 <- dat01[!is.na(dat01$Species), ]
checkDat(dat02, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows <- c(
  401, 1010, 1661, 24487, 25606, 25670, 26255, 32384, 38391, 39569, 55076,
  55089)

dat03 <- dat02[!(row.names(dat02) %in% drop_rows),]
checkDat(dat03, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

# Still have a couple of repeated rows, somehow, so we will drop those
drop_rows2 <- c(
  55112, 58091, 64803, 72238, 75507, 77291, 79196, 96626, 103604, 116003, 
  128419)

dat03 <- dat03[!(row.names(dat03) %in% drop_rows2),]
checkDat(dat03, inv_sgs, species = "Species", site = "Site", quad = "Quad", year = "Year", geometry = "geometry")

saveRDS(dat03, file = paste0(dir_qud, "SGS_LTER_plantTracker_all_filtered.rds"))

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
