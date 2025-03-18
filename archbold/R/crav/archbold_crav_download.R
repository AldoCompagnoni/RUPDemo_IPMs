# EDI code to download the Crotalaria avonensis dataset

# Package ID: edi.219.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Long-term demographic data on the endangered legume Crotalaria avonensis, from Carter Creek Tract, Lake Wales Ridge Wildlife Environmental Area, 1998-2017 (ongoing).
# Data set creator:  Eric Menges - Archbold Biological Station 
# Data set creator:  Stacy Smith - Archbold Biological Station 
# Data set creator:  Sarah Crate - Archbold Biological Station 
# Contact:    - Data Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.14 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       
options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/219/1/b33733fbf967c66161d4472e538cc036" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "site",     
                 "quad",     
                 "mp",     
                 "plant",     
                 "bearing",     
                 "distance",     
                 "caged",     
                 "veg",     
                 "quadyr",     
                 "date",     
                 "burnA",     
                 "burnB",     
                 "burnC",     
                 "burnD",     
                 "burnE",     
                 "burnF",     
                 "s",     
                 "st",     
                 "br",     
                 "fl",     
                 "dv",     
                 "fr",     
                 "h"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$site)=="factor") dt1$site <-as.numeric(levels(dt1$site))[as.integer(dt1$site) ]               
if (class(dt1$site)=="character") dt1$site <-as.numeric(dt1$site)
if (class(dt1$quad)=="factor") dt1$quad <-as.numeric(levels(dt1$quad))[as.integer(dt1$quad) ]               
if (class(dt1$quad)=="character") dt1$quad <-as.numeric(dt1$quad)
if (class(dt1$mp)=="factor") dt1$mp <-as.numeric(levels(dt1$mp))[as.integer(dt1$mp) ]               
if (class(dt1$mp)=="character") dt1$mp <-as.numeric(dt1$mp)
if (class(dt1$plant)=="factor") dt1$plant <-as.numeric(levels(dt1$plant))[as.integer(dt1$plant) ]               
if (class(dt1$plant)=="character") dt1$plant <-as.numeric(dt1$plant)
if (class(dt1$bearing)=="factor") dt1$bearing <-as.numeric(levels(dt1$bearing))[as.integer(dt1$bearing) ]               
if (class(dt1$bearing)=="character") dt1$bearing <-as.numeric(dt1$bearing)
if (class(dt1$distance)=="factor") dt1$distance <-as.numeric(levels(dt1$distance))[as.integer(dt1$distance) ]               
if (class(dt1$distance)=="character") dt1$distance <-as.numeric(dt1$distance)
if (class(dt1$caged)!="factor") dt1$caged<- as.factor(dt1$caged)
if (class(dt1$veg)!="factor") dt1$veg<- as.factor(dt1$veg)
if (class(dt1$burnA)!="factor") dt1$burnA<- as.factor(dt1$burnA)
if (class(dt1$burnB)!="factor") dt1$burnB<- as.factor(dt1$burnB)
if (class(dt1$burnC)!="factor") dt1$burnC<- as.factor(dt1$burnC)
if (class(dt1$burnD)!="factor") dt1$burnD<- as.factor(dt1$burnD)
if (class(dt1$burnE)!="factor") dt1$burnE<- as.factor(dt1$burnE)
if (class(dt1$burnF)!="factor") dt1$burnF<- as.factor(dt1$burnF)
if (class(dt1$s)!="factor") dt1$s<- as.factor(dt1$s)
if (class(dt1$st)=="factor") dt1$st <-as.numeric(levels(dt1$st))[as.integer(dt1$st) ]               
if (class(dt1$st)=="character") dt1$st <-as.numeric(dt1$st)
if (class(dt1$br)=="factor") dt1$br <-as.numeric(levels(dt1$br))[as.integer(dt1$br) ]               
if (class(dt1$br)=="character") dt1$br <-as.numeric(dt1$br)
if (class(dt1$fl)=="factor") dt1$fl <-as.numeric(levels(dt1$fl))[as.integer(dt1$fl) ]               
if (class(dt1$fl)=="character") dt1$fl <-as.numeric(dt1$fl)
if (class(dt1$dv)=="factor") dt1$dv <-as.numeric(levels(dt1$dv))[as.integer(dt1$dv) ]               
if (class(dt1$dv)=="character") dt1$dv <-as.numeric(dt1$dv)
if (class(dt1$fr)=="factor") dt1$fr <-as.numeric(levels(dt1$fr))[as.integer(dt1$fr) ]               
if (class(dt1$fr)=="character") dt1$fr <-as.numeric(dt1$fr)
if (class(dt1$h)!="factor") dt1$h<- as.factor(dt1$h)

# Convert Missing Values to NA for non-dates

dt1$bearing <- ifelse((trimws(as.character(dt1$bearing))==trimws("9999")),NA,dt1$bearing)               
suppressWarnings(dt1$bearing <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$bearing))==as.character(as.numeric("9999"))),NA,dt1$bearing))
dt1$distance <- ifelse((trimws(as.character(dt1$distance))==trimws("9999")),NA,dt1$distance)               
suppressWarnings(dt1$distance <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$distance))==as.character(as.numeric("9999"))),NA,dt1$distance))
dt1$burnA <- as.factor(ifelse((trimws(as.character(dt1$burnA))==trimws("9999")),NA,as.character(dt1$burnA)))
dt1$burnB <- as.factor(ifelse((trimws(as.character(dt1$burnB))==trimws("9999")),NA,as.character(dt1$burnB)))
dt1$burnC <- as.factor(ifelse((trimws(as.character(dt1$burnC))==trimws("9999")),NA,as.character(dt1$burnC)))
dt1$burnD <- as.factor(ifelse((trimws(as.character(dt1$burnD))==trimws("9999")),NA,as.character(dt1$burnD)))
dt1$burnE <- as.factor(ifelse((trimws(as.character(dt1$burnE))==trimws("9999")),NA,as.character(dt1$burnE)))
dt1$burnF <- as.factor(ifelse((trimws(as.character(dt1$burnF))==trimws("9999")),NA,as.character(dt1$burnF)))
dt1$s <- as.factor(ifelse((trimws(as.character(dt1$s))==trimws("9999")),NA,as.character(dt1$s)))
dt1$st <- ifelse((trimws(as.character(dt1$st))==trimws("9999")),NA,dt1$st)               
suppressWarnings(dt1$st <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$st))==as.character(as.numeric("9999"))),NA,dt1$st))
dt1$br <- ifelse((trimws(as.character(dt1$br))==trimws("9999")),NA,dt1$br)               
suppressWarnings(dt1$br <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$br))==as.character(as.numeric("9999"))),NA,dt1$br))
dt1$fl <- ifelse((trimws(as.character(dt1$fl))==trimws("9999")),NA,dt1$fl)               
suppressWarnings(dt1$fl <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$fl))==as.character(as.numeric("9999"))),NA,dt1$fl))
dt1$dv <- ifelse((trimws(as.character(dt1$dv))==trimws("9999")),NA,dt1$dv)               
suppressWarnings(dt1$dv <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$dv))==as.character(as.numeric("9999"))),NA,dt1$dv))
dt1$fr <- ifelse((trimws(as.character(dt1$fr))==trimws("9999")),NA,dt1$fr)               
suppressWarnings(dt1$fr <- ifelse(!is.na(as.numeric("9999")) & (trimws(as.character(dt1$fr))==as.character(as.numeric("9999"))),NA,dt1$fr))
dt1$h <- as.factor(ifelse((trimws(as.character(dt1$h))==trimws("9999")),NA,as.character(dt1$h)))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(site)
summary(quad)
summary(mp)
summary(plant)
summary(bearing)
summary(distance)
summary(caged)
summary(veg)
summary(quadyr)
summary(date)
summary(burnA)
summary(burnB)
summary(burnC)
summary(burnD)
summary(burnE)
summary(burnF)
summary(s)
summary(st)
summary(br)
summary(fl)
summary(dv)
summary(fr)
summary(h) 
# Get more details on character variables

summary(as.factor(dt1$caged)) 
summary(as.factor(dt1$veg)) 
summary(as.factor(dt1$burnA)) 
summary(as.factor(dt1$burnB)) 
summary(as.factor(dt1$burnC)) 
summary(as.factor(dt1$burnD)) 
summary(as.factor(dt1$burnE)) 
summary(as.factor(dt1$burnF)) 
summary(as.factor(dt1$s)) 
summary(as.factor(dt1$h))
detach(dt1)               

# Store Crotalaria avonensis
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Crotalaria avonensis')
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
           paste0(dir_data, '/crotalaria_avonensis_data.csv'), 
           row.names = F )

