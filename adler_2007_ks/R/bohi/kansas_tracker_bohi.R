# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)

library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0

dir     <- 'adler_2007_ks/data/quadrat_data/'
dat_dir <- dir
shp_dir <- paste0(dir, "arcexport/")

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# Read in species list, species name changes, and subset species list to perennial target_spec
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list         <- read.csv(paste0(dat_dir,"species_list.csv")) 
sp_list %>% dplyr::arrange( desc(count) ) %>% head(20)
target_spec     <- sp_list %>% 
                    dplyr::arrange( desc(count) ) %>% 
                    .[c(10),]
# Define the species variable
species <- 'Bouteloua hirsuta'
sp_abb  <- tolower(gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                                       collapse = "")))

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv        <- read.csv(paste0(dat_dir,"quadrat_inventory.csv"),
                            sep=',') %>% 
                    dplyr::select(-year)
quadInv_list    <- as.list(quad_inv)
quadInv_list    <- lapply(X = quadInv_list, 
                          FUN = function(x) x[is.na(x) == FALSE])
inv_ks          <- quadInv_list 
names(inv_ks)   <- gsub( '\\.','-',names(inv_ks) )

# Read in the data
dat <- readRDS(file=paste0(dat_dir,"KS_polygons_full.rds"))

# Subset to the species of interest
dat_target_spec <- dat[dat$SCI_NAME %in% target_spec$species,] %>% 
  setNames( quote_bare(Species, Site, Quad, Year, geometry) )

# Now the data are ready for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec,
                        inv_ks,
                        dorm=1,
                        buff=0.05,
                        clonal=TRUE,
                        buffGenet = 0.05,
                        aggByGenet = TRUE,
                        flagSuspects = TRUE)

# create folder
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))
}

# save data
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select( Site, Quad, Species, trackID,
                 Year, basalArea_genet, recruit,
                 survives_tplus1, age, size_tplus1,
                 nearEdge, Suspect) %>% 
  write.csv('adler_2007_ks/data/bohi/ks_bohi.csv', row.names = F )
