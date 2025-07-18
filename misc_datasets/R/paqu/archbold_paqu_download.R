# Package ID: edi.9.4 Cataloging System:https://pasta.edirepository.org.
# Data set title: Long Term Research in Environmental Biology: Demographic census data for thirty natural populations of American Ginseng: 1998-2016.
# Data set creator:  James McGraw - West Virginia University 
# Data set creator:  Martha Van der Voort - Adirondack Center for Loon Conservation 
# Data set creator:  Mary Ann Furedi - Western Pennsylvania Conservancy 
# Data set creator:  Anne Lubbers - Centre College 
# Data set creator:  Emily Mooney - University of Colorado - Colorado Springs 
# Data set creator:  Sara Souther - West Virginia University 
# Data set creator:  Jessica Turner-Skoff - Morton Arboretum 
# Data set creator:  Jennifer Chandler - Appalachian State University 
# Data set creator:  Emily Thyroff - Purdue University 
# Contact:  James McGraw -  West Virginia University  - jmcgraw56@gmail.com
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       



options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/9/4/480547e45c7b8beb3f70b104d8de0bd6" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "population",     
                 "year",     
                 "id",     
                 "cluster",     
                 "age",     
                 "obs_leaf_num",     
                 "inf_leaf_num",     
                 "lflt_arr",     
                 "lflt_tot",     
                 "stalk_ht",     
                 "lll1",     
                 "wll1",     
                 "lll2",     
                 "wll2",     
                 "lll3",     
                 "wll3",     
                 "lll4",     
                 "wll4",     
                 "obs_la",     
                 "inf_la",     
                 "f_buds",     
                 "seeds",     
                 "red_frts",     
                 "grn_frts",     
                 "tot_frts",     
                 "loc",     
                 "harvest",     
                 "harv_time",     
                 "browse",     
                 "browse_prcnt",     
                 "browse_time",     
                 "persistence",     
                 "insect",     
                 "thrips",     
                 "fungal",     
                 "fungal_cat"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$population)=="factor") dt1$population <-as.numeric(levels(dt1$population))[as.integer(dt1$population) ]               
if (class(dt1$population)=="character") dt1$population <-as.numeric(dt1$population)
if (class(dt1$id)=="factor") dt1$id <-as.numeric(levels(dt1$id))[as.integer(dt1$id) ]               
if (class(dt1$id)=="character") dt1$id <-as.numeric(dt1$id)
if (class(dt1$cluster)=="factor") dt1$cluster <-as.numeric(levels(dt1$cluster))[as.integer(dt1$cluster) ]               
if (class(dt1$cluster)=="character") dt1$cluster <-as.numeric(dt1$cluster)
if (class(dt1$age)=="factor") dt1$age <-as.numeric(levels(dt1$age))[as.integer(dt1$age) ]               
if (class(dt1$age)=="character") dt1$age <-as.numeric(dt1$age)
if (class(dt1$obs_leaf_num)=="factor") dt1$obs_leaf_num <-as.numeric(levels(dt1$obs_leaf_num))[as.integer(dt1$obs_leaf_num) ]               
if (class(dt1$obs_leaf_num)=="character") dt1$obs_leaf_num <-as.numeric(dt1$obs_leaf_num)
if (class(dt1$inf_leaf_num)=="factor") dt1$inf_leaf_num <-as.numeric(levels(dt1$inf_leaf_num))[as.integer(dt1$inf_leaf_num) ]               
if (class(dt1$inf_leaf_num)=="character") dt1$inf_leaf_num <-as.numeric(dt1$inf_leaf_num)
if (class(dt1$lflt_arr)!="factor") dt1$lflt_arr<- as.factor(dt1$lflt_arr)
if (class(dt1$lflt_tot)=="factor") dt1$lflt_tot <-as.numeric(levels(dt1$lflt_tot))[as.integer(dt1$lflt_tot) ]               
if (class(dt1$lflt_tot)=="character") dt1$lflt_tot <-as.numeric(dt1$lflt_tot)
if (class(dt1$stalk_ht)=="factor") dt1$stalk_ht <-as.numeric(levels(dt1$stalk_ht))[as.integer(dt1$stalk_ht) ]               
if (class(dt1$stalk_ht)=="character") dt1$stalk_ht <-as.numeric(dt1$stalk_ht)
if (class(dt1$lll1)=="factor") dt1$lll1 <-as.numeric(levels(dt1$lll1))[as.integer(dt1$lll1) ]               
if (class(dt1$lll1)=="character") dt1$lll1 <-as.numeric(dt1$lll1)
if (class(dt1$wll1)=="factor") dt1$wll1 <-as.numeric(levels(dt1$wll1))[as.integer(dt1$wll1) ]               
if (class(dt1$wll1)=="character") dt1$wll1 <-as.numeric(dt1$wll1)
if (class(dt1$lll2)=="factor") dt1$lll2 <-as.numeric(levels(dt1$lll2))[as.integer(dt1$lll2) ]               
if (class(dt1$lll2)=="character") dt1$lll2 <-as.numeric(dt1$lll2)
if (class(dt1$wll2)=="factor") dt1$wll2 <-as.numeric(levels(dt1$wll2))[as.integer(dt1$wll2) ]               
if (class(dt1$wll2)=="character") dt1$wll2 <-as.numeric(dt1$wll2)
if (class(dt1$lll3)=="factor") dt1$lll3 <-as.numeric(levels(dt1$lll3))[as.integer(dt1$lll3) ]               
if (class(dt1$lll3)=="character") dt1$lll3 <-as.numeric(dt1$lll3)
if (class(dt1$wll3)=="factor") dt1$wll3 <-as.numeric(levels(dt1$wll3))[as.integer(dt1$wll3) ]               
if (class(dt1$wll3)=="character") dt1$wll3 <-as.numeric(dt1$wll3)
if (class(dt1$lll4)=="factor") dt1$lll4 <-as.numeric(levels(dt1$lll4))[as.integer(dt1$lll4) ]               
if (class(dt1$lll4)=="character") dt1$lll4 <-as.numeric(dt1$lll4)
if (class(dt1$wll4)=="factor") dt1$wll4 <-as.numeric(levels(dt1$wll4))[as.integer(dt1$wll4) ]               
if (class(dt1$wll4)=="character") dt1$wll4 <-as.numeric(dt1$wll4)
if (class(dt1$obs_la)=="factor") dt1$obs_la <-as.numeric(levels(dt1$obs_la))[as.integer(dt1$obs_la) ]               
if (class(dt1$obs_la)=="character") dt1$obs_la <-as.numeric(dt1$obs_la)
if (class(dt1$inf_la)=="factor") dt1$inf_la <-as.numeric(levels(dt1$inf_la))[as.integer(dt1$inf_la) ]               
if (class(dt1$inf_la)=="character") dt1$inf_la <-as.numeric(dt1$inf_la)
if (class(dt1$f_buds)!="factor") dt1$f_buds<- as.factor(dt1$f_buds)
if (class(dt1$seeds)=="factor") dt1$seeds <-as.numeric(levels(dt1$seeds))[as.integer(dt1$seeds) ]               
if (class(dt1$seeds)=="character") dt1$seeds <-as.numeric(dt1$seeds)
if (class(dt1$red_frts)=="factor") dt1$red_frts <-as.numeric(levels(dt1$red_frts))[as.integer(dt1$red_frts) ]               
if (class(dt1$red_frts)=="character") dt1$red_frts <-as.numeric(dt1$red_frts)
if (class(dt1$grn_frts)=="factor") dt1$grn_frts <-as.numeric(levels(dt1$grn_frts))[as.integer(dt1$grn_frts) ]               
if (class(dt1$grn_frts)=="character") dt1$grn_frts <-as.numeric(dt1$grn_frts)
if (class(dt1$tot_frts)=="factor") dt1$tot_frts <-as.numeric(levels(dt1$tot_frts))[as.integer(dt1$tot_frts) ]               
if (class(dt1$tot_frts)=="character") dt1$tot_frts <-as.numeric(dt1$tot_frts)
if (class(dt1$loc)!="factor") dt1$loc<- as.factor(dt1$loc)
if (class(dt1$harvest)!="factor") dt1$harvest<- as.factor(dt1$harvest)
if (class(dt1$harv_time)!="factor") dt1$harv_time<- as.factor(dt1$harv_time)
if (class(dt1$browse)!="factor") dt1$browse<- as.factor(dt1$browse)
if (class(dt1$browse_prcnt)=="factor") dt1$browse_prcnt <-as.numeric(levels(dt1$browse_prcnt))[as.integer(dt1$browse_prcnt) ]               
if (class(dt1$browse_prcnt)=="character") dt1$browse_prcnt <-as.numeric(dt1$browse_prcnt)
if (class(dt1$browse_time)!="factor") dt1$browse_time<- as.factor(dt1$browse_time)
if (class(dt1$persistence)!="factor") dt1$persistence<- as.factor(dt1$persistence)
if (class(dt1$insect)!="factor") dt1$insect<- as.factor(dt1$insect)
if (class(dt1$thrips)!="factor") dt1$thrips<- as.factor(dt1$thrips)
if (class(dt1$fungal)!="factor") dt1$fungal<- as.factor(dt1$fungal)
if (class(dt1$fungal_cat)!="factor") dt1$fungal_cat<- as.factor(dt1$fungal_cat)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(population)
summary(year)
summary(id)
summary(cluster)
summary(age)
summary(obs_leaf_num)
summary(inf_leaf_num)
summary(lflt_arr)
summary(lflt_tot)
summary(stalk_ht)
summary(lll1)
summary(wll1)
summary(lll2)
summary(wll2)
summary(lll3)
summary(wll3)
summary(lll4)
summary(wll4)
summary(obs_la)
summary(inf_la)
summary(f_buds)
summary(seeds)
summary(red_frts)
summary(grn_frts)
summary(tot_frts)
summary(loc)
summary(harvest)
summary(harv_time)
summary(browse)
summary(browse_prcnt)
summary(browse_time)
summary(persistence)
summary(insect)
summary(thrips)
summary(fungal)
summary(fungal_cat) 
# Get more details on character variables

summary(as.factor(dt1$lflt_arr)) 
summary(as.factor(dt1$f_buds)) 
summary(as.factor(dt1$loc)) 
summary(as.factor(dt1$harvest)) 
summary(as.factor(dt1$harv_time)) 
summary(as.factor(dt1$browse)) 
summary(as.factor(dt1$browse_time)) 
summary(as.factor(dt1$persistence)) 
summary(as.factor(dt1$insect)) 
summary(as.factor(dt1$thrips)) 
summary(as.factor(dt1$fungal)) 
summary(as.factor(dt1$fungal_cat))
detach(dt1)               



inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/9/4/2d1981951436bfa5be2354f3d31906f1" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "population",     
                 "state",     
                 "north_bounding_coordinate",     
                 "south_bounding_coordinate",     
                 "east_bounding_coordinate",     
                 "west_bounding_coordinate",     
                 "elevation"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$population)=="factor") dt2$population <-as.numeric(levels(dt2$population))[as.integer(dt2$population) ]               
if (class(dt2$population)=="character") dt2$population <-as.numeric(dt2$population)
if (class(dt2$state)!="factor") dt2$state<- as.factor(dt2$state)
if (class(dt2$north_bounding_coordinate)=="factor") dt2$north_bounding_coordinate <-as.numeric(levels(dt2$north_bounding_coordinate))[as.integer(dt2$north_bounding_coordinate) ]               
if (class(dt2$north_bounding_coordinate)=="character") dt2$north_bounding_coordinate <-as.numeric(dt2$north_bounding_coordinate)
if (class(dt2$south_bounding_coordinate)=="factor") dt2$south_bounding_coordinate <-as.numeric(levels(dt2$south_bounding_coordinate))[as.integer(dt2$south_bounding_coordinate) ]               
if (class(dt2$south_bounding_coordinate)=="character") dt2$south_bounding_coordinate <-as.numeric(dt2$south_bounding_coordinate)
if (class(dt2$east_bounding_coordinate)=="factor") dt2$east_bounding_coordinate <-as.numeric(levels(dt2$east_bounding_coordinate))[as.integer(dt2$east_bounding_coordinate) ]               
if (class(dt2$east_bounding_coordinate)=="character") dt2$east_bounding_coordinate <-as.numeric(dt2$east_bounding_coordinate)
if (class(dt2$west_bounding_coordinate)=="factor") dt2$west_bounding_coordinate <-as.numeric(levels(dt2$west_bounding_coordinate))[as.integer(dt2$west_bounding_coordinate) ]               
if (class(dt2$west_bounding_coordinate)=="character") dt2$west_bounding_coordinate <-as.numeric(dt2$west_bounding_coordinate)
if (class(dt2$elevation)=="factor") dt2$elevation <-as.numeric(levels(dt2$elevation))[as.integer(dt2$elevation) ]               
if (class(dt2$elevation)=="character") dt2$elevation <-as.numeric(dt2$elevation)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt2)                            
attach(dt2)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(population)
summary(state)
summary(north_bounding_coordinate)
summary(south_bounding_coordinate)
summary(east_bounding_coordinate)
summary(west_bounding_coordinate)
summary(elevation) 
# Get more details on character variables

summary(as.factor(dt2$state))
detach(dt2)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('misc_datasets')
# Define species
v_species <- c('Panax quinquefolius')
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c()


# Create a unique species abbreviation for file naming
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Define script prefix
v_script_prefix <- str_c(v_head)


# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_head))
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)

# test if the directory exists
if(!dir.exists(dir_data) ) dir.create(dir_data)

# Store Crotalaria avonensis 
write.csv( dt1, 
           paste0(dir_data, '/panax_quinquefolius_data.csv'), 
           row.names = F )
