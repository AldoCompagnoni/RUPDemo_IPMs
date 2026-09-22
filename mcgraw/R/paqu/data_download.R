# IMPORTANT: Due to new authentication requirements for direct access to EDI data you must provide an 
# authentication key. You can find instructions for obtaining a group-based access key at:
# https://youtu.be/fieZSmHk2H4?si=Wo9a5GsAOYp3dnWS 
# and the IAM web site at: 
# https://auth.edirepository.org
#
# Once you have the access key, just paste it between the quotes below and the code should run correctly and
# automatically download data from EDI for the code to use. 

readRenviron(path.expand("~/.Renviron"))

myEDIAccessKey <- Sys.getenv("EDI_ACCESS_KEY", unset = "")   # put your access key between the quotes 

#
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
# Stylesheet v2.17 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       

# check for required access key
if (myEDIAccessKey == ""){
  stop(paste0("A VALID EDI ACCESS KEY IS REQUIRED FOR THIS CODE TO DOWNLOAD THE NEEDED DATA \n",
              "Please set the myEDIAccessKey variable to a valid key and rerun this script"))
}

options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/9/4/480547e45c7b8beb3f70b104d8de0bd6","?key=",myEDIAccessKey) 
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
print("dt1) Structure")		    
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 


print(" ")
print("Summary of population")
print(summary(population))
print(" ")
print("Summary of year")
print(summary(year))
print(" ")
print("Summary of id")
print(summary(id))
print(" ")
print("Summary of cluster")
print(summary(cluster))
print(" ")
print("Summary of age")
print(summary(age))
print(" ")
print("Summary of obs_leaf_num")
print(summary(obs_leaf_num))
print(" ")
print("Summary of inf_leaf_num")
print(summary(inf_leaf_num))
print(" ")
print("Summary of lflt_arr")
print(summary(lflt_arr))
print(" ")
print("Summary of lflt_tot")
print(summary(lflt_tot))
print(" ")
print("Summary of stalk_ht")
print(summary(stalk_ht))
print(" ")
print("Summary of lll1")
print(summary(lll1))
print(" ")
print("Summary of wll1")
print(summary(wll1))
print(" ")
print("Summary of lll2")
print(summary(lll2))
print(" ")
print("Summary of wll2")
print(summary(wll2))
print(" ")
print("Summary of lll3")
print(summary(lll3))
print(" ")
print("Summary of wll3")
print(summary(wll3))
print(" ")
print("Summary of lll4")
print(summary(lll4))
print(" ")
print("Summary of wll4")
print(summary(wll4))
print(" ")
print("Summary of obs_la")
print(summary(obs_la))
print(" ")
print("Summary of inf_la")
print(summary(inf_la))
print(" ")
print("Summary of f_buds")
print(summary(f_buds))
print(" ")
print("Summary of seeds")
print(summary(seeds))
print(" ")
print("Summary of red_frts")
print(summary(red_frts))
print(" ")
print("Summary of grn_frts")
print(summary(grn_frts))
print(" ")
print("Summary of tot_frts")
print(summary(tot_frts))
print(" ")
print("Summary of loc")
print(summary(loc))
print(" ")
print("Summary of harvest")
print(summary(harvest))
print(" ")
print("Summary of harv_time")
print(summary(harv_time))
print(" ")
print("Summary of browse")
print(summary(browse))
print(" ")
print("Summary of browse_prcnt")
print(summary(browse_prcnt))
print(" ")
print("Summary of browse_time")
print(summary(browse_time))
print(" ")
print("Summary of persistence")
print(summary(persistence))
print(" ")
print("Summary of insect")
print(summary(insect))
print(" ")
print("Summary of thrips")
print(summary(thrips))
print(" ")
print("Summary of fungal")
print(summary(fungal))
print(" ")
print("Summary of fungal_cat")
print(summary(fungal_cat)) 
# Get more details on character variables


print(" ")
print("Summary of lflt_arr")
print(summary(as.factor(dt1$lflt_arr))) 

print(" ")
print("Summary of f_buds")
print(summary(as.factor(dt1$f_buds))) 

print(" ")
print("Summary of loc")
print(summary(as.factor(dt1$loc))) 

print(" ")
print("Summary of harvest")
print(summary(as.factor(dt1$harvest))) 

print(" ")
print("Summary of harv_time")
print(summary(as.factor(dt1$harv_time))) 

print(" ")
print("Summary of browse")
print(summary(as.factor(dt1$browse))) 

print(" ")
print("Summary of browse_time")
print(summary(as.factor(dt1$browse_time))) 

print(" ")
print("Summary of persistence")
print(summary(as.factor(dt1$persistence))) 

print(" ")
print("Summary of insect")
print(summary(as.factor(dt1$insect))) 

print(" ")
print("Summary of thrips")
print(summary(as.factor(dt1$thrips))) 

print(" ")
print("Summary of fungal")
print(summary(as.factor(dt1$fungal))) 

print(" ")
print("Summary of fungal_cat")
print(summary(as.factor(dt1$fungal_cat)))
detach(dt1)               



inUrl2  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/9/4/2d1981951436bfa5be2354f3d31906f1","?key=",myEDIAccessKey) 
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
print("dt2) Structure")		    
str(dt2)                            
attach(dt2)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 


print(" ")
print("Summary of population")
print(summary(population))
print(" ")
print("Summary of state")
print(summary(state))
print(" ")
print("Summary of north_bounding_coordinate")
print(summary(north_bounding_coordinate))
print(" ")
print("Summary of south_bounding_coordinate")
print(summary(south_bounding_coordinate))
print(" ")
print("Summary of east_bounding_coordinate")
print(summary(east_bounding_coordinate))
print(" ")
print("Summary of west_bounding_coordinate")
print(summary(west_bounding_coordinate))
print(" ")
print("Summary of elevation")
print(summary(elevation)) 
# Get more details on character variables


print(" ")
print("Summary of state")
print(summary(as.factor(dt2$state)))
detach(dt2)


# Saving data ------------------------------------------------------------------
write.csv(dt1, file.path("C:/code/RUPDemo_IPMs/mcgraw/paqu/data","df_orig.csv"), row.names = F)
write.csv(dt2, file.path("C:/code/RUPDemo_IPMs/mcgraw/paqu/data","df_meta.csv"), row.names = F)
