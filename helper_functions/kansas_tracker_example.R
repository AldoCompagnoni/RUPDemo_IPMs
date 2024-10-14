# Code adapted from: https://github.com/aestears/plantTracker
library(sf) #ver 1.0-1.2
library(plantTracker) #ver 1.1.0

dir     <- 'C:/Users/ac22qawo/Downloads/3528368/'
dat_dir <- dir
shp_dir <- paste0(dir, "arcexport/")

setwd(dat_dir)

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# Read in species list, species name changes, and subset species list to perennial grasses
# with minimum cover of 100. Also taking out Carex spp.; 8 species total, might exclude some
# species with the lowest cover later.
sp_list         <- read.csv("species_list.csv") 
grasses         <- sp_list %>% 
                    dplyr::arrange( desc(count) ) %>% 
                    .[1:3,]
grasses         <- sp_list %>% 
                    dplyr::arrange( desc(count) ) %>% 
                    .[c(5,8),]
# grasses         <- subset(sp_list, growthForm=="grass" & longevity=="P" & cover>100 & species!="Carex spp.")

# Read in quad inventory to use as 'inv' list in plantTracker
quad_inv        <- read.csv("quadrat_inventory.csv",
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
# quadNames <- quadNames[2]

# read in the GIS files
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
    if(j == 1) {
      
      dat <- shapeNow
      
      # "append" data to the initial data frame  
    } else {
      
      dat <- rbind(dat, shapeNow)
      
    }
  }

  # 
  dat$Year <- as.numeric(dat$Year)
  
  # Subset to the species of interest
  dat_3grasses <- dat[dat$SCI_NAME %in% grasses$species,]
  names(dat_3grasses) <- quote_bare(Species, Site, Quad, Year, geometry )
  
  # Check the inv and dat arguments
  checkDat(dat_3grasses, 
           inv_ks[quadNow], 
           species  = "Species", 
           site     = "Site", 
           quad     = "Quad", 
           year     = "Year", 
           geometry = "geometry")
  
  
  # Now the data are ready for the trackSpp function
  datTrackSpp <- trackSpp(dat_3grasses,
                          inv_ks[quadNow], 
                          dorm=1,
                          buff=0.05,
                          clonal=TRUE,
                          buffGenet = 0.05,
                          aggByGenet = TRUE,
                          flagSuspects = TRUE)
  
  saveRDS(datTrackSpp,
          file= paste0("KS_grasses_",quadNow,".rds") )
  
  rm(dat)
  rm(dat_3grasses)
  
}


# Save the output file so that it doesn't need to be recreated ever again
# saveRDS(dat,file="KS_polygons_full.rds")
# dat <- readRDS(file="KS_polygons_full.rds")

# Subset to the species of interest
dat_3grasses <- dat[dat$SCI_NAME %in% grasses$species,] %>% 
                  setNames( quote_bare(Species, Site, Quad, 
                                       Year, geometry) )
# dat_3grasses$Type <- rep('polygon',93640 )


# And save the subsetted file, too
# saveRDS(dat_3grasses,file="KS_polygons_3grasses.rds")
# dat_3grasses <- readRDS( file="KS_polygons_3grasses.rds" )

# Check the inv and dat arguments
checkDat(dat_3grasses, 
         inv_ks, 
         species  = "Species", 
         site     = "Site", 
         quad     = "Quad", 
         year     = "Year", 
         geometry = "geometry")


# Now the data are ready for the trackSpp function
datTrackSpp <- trackSpp(dat_3grasses,
                        inv_ks,
                        dorm=1,
                        buff=0.05,
                        clonal=TRUE,
                        buffGenet = 0.05,
                        aggByGenet = TRUE,
                        flagSuspects = TRUE)
datTrackSpp %>% 
  as.data.frame %>% 
  dplyr::select( Site, Quad, Species, trackID,
                 Year, basalArea_genet, recruit,
                 survives_tplus1, age, size_tplus1,
                 nearEdge, Suspect) %>% 
  write.csv( 'KS_SCSC_ANGE.csv', row.names = F )
