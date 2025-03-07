# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)

library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0
library(tidyverse) #ver 1.1.0

dir     <- 'anderson_2016_mt/data/quadrat_data/'
dat_dir <- dir

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# Read in species list, species name changes, and subset species list to perennial target_spec
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list         <- read.csv( paste0(dat_dir,"species_list.csv") ) %>% 
  mutate( species = paste(species,X) ) %>% 
  dplyr::select(-X)
#arrange species by density and cover (deeded in anderson_mt for precision)
sp_list %>% dplyr::arrange(desc(cover), desc(density)) %>% head(50)
target_spec     <- sp_list %>% 
  dplyr::arrange(desc(cover), desc(density)) %>% 
  .[c(32),]
# Define the species variable
species <- 'Carex duriuscula'
sp_abb  <- tolower(gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                                       collapse = "")))

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv        <- read.csv(paste0(dat_dir,"quad_inventory.csv"),
                            sep=',') 
#%>% 
#dplyr::select(-Year)


quad_inv <- read.csv(paste0(dat_dir,"quad_inventory.csv"), sep=',', stringsAsFactors = FALSE)
names(quad_inv)

quadInv_list    <- as.list(quad_inv)
quadInv_list    <- lapply(X = quadInv_list, 
                          FUN = function(x) x[is.na(x) == FALSE])
inv_mt          <- quadInv_list 
names(inv_mt)   <- gsub( '\\.','-',names(inv_mt) )

# Read in the data
dat <- readRDS(file=paste0(dat_dir,"anderson16mt_plantTracker_all_filtered.rds"))

# Subset to the species of interest
dat_target_spec <- dat[dat$species %in% target_spec$species,] %>% 
  setNames( quote_bare(Species, Site, Quad, Year, type, geometry) )

# Now the data are ready for the trackSpp function
datTrackSpp <- trackSpp(dat_target_spec,
                        inv_mt,
                        dorm=1,
                        buff= 0.05,
                        clonal=TRUE,
                        buffGenet = 0.05,
                        aggByGenet = TRUE,
                        flagSuspects = TRUE)



# create folder
if (!dir.exists(paste0("anderson_2016_mt/data/", sp_abb))) {
  dir.create(paste0("anderson_2016_mt/data/", sp_abb))
}

# save data
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select( Site, Quad, Species, trackID,
                 Year, basalArea_genet, recruit,
                 survives_tplus1, age, size_tplus1,
                 nearEdge, Suspect) %>% 
  write.csv( paste0('anderson_2016_mt/data/',
                    sp_abb,'/anderson16mt_',sp_abb,'.csv'), 
             row.names = F )
