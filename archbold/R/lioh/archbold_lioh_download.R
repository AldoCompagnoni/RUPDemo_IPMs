# Package ID: edi.240.2 Cataloging System:https://pasta.edirepository.org.
# Data set title: Demographic measures of Liatris ohlingerae (Asteraceae) in 20 populations          across multiple habitats and time-since-fire intervals in south central Florida from 1997-2017.
# Data set creator:  Eric Menges - Archbold Biological Station 
# Data set creator:  Carl Weekley - Archbold Biological Station 
# Contact:    - Data Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       



options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/240/2/63da0f716a599f3e84516e87f51d589a" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "site",     
                 "population",     
                 "toothpick",     
                 "tag",     
                 "year",     
                 "hurricane",     
                 "habitat",     
                 "burn",     
                 "survival",     
                 "life_stage",     
                 "reproductive_stems",     
                 "topped_stems",     
                 "total_stem_height",     
                 "total_reproductive_heads",     
                 "damaged_heads",     
                 "rosettes",     
                 "leaves",     
                 "comments",     
                 "QAQCcomments"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$site)!="factor") dt1$site<- as.factor(dt1$site)
if (class(dt1$population)!="factor") dt1$population<- as.factor(dt1$population)
if (class(dt1$toothpick)!="factor") dt1$toothpick<- as.factor(dt1$toothpick)
if (class(dt1$tag)=="factor") dt1$tag <-as.numeric(levels(dt1$tag))[as.integer(dt1$tag) ]               
if (class(dt1$tag)=="character") dt1$tag <-as.numeric(dt1$tag)
if (class(dt1$hurricane)!="factor") dt1$hurricane<- as.factor(dt1$hurricane)
if (class(dt1$habitat)!="factor") dt1$habitat<- as.factor(dt1$habitat)
if (class(dt1$burn)!="factor") dt1$burn<- as.factor(dt1$burn)
if (class(dt1$survival)!="factor") dt1$survival<- as.factor(dt1$survival)
if (class(dt1$life_stage)!="factor") dt1$life_stage<- as.factor(dt1$life_stage)
if (class(dt1$reproductive_stems)=="factor") dt1$reproductive_stems <-as.numeric(levels(dt1$reproductive_stems))[as.integer(dt1$reproductive_stems) ]               
if (class(dt1$reproductive_stems)=="character") dt1$reproductive_stems <-as.numeric(dt1$reproductive_stems)
if (class(dt1$topped_stems)=="factor") dt1$topped_stems <-as.numeric(levels(dt1$topped_stems))[as.integer(dt1$topped_stems) ]               
if (class(dt1$topped_stems)=="character") dt1$topped_stems <-as.numeric(dt1$topped_stems)
if (class(dt1$total_stem_height)=="factor") dt1$total_stem_height <-as.numeric(levels(dt1$total_stem_height))[as.integer(dt1$total_stem_height) ]               
if (class(dt1$total_stem_height)=="character") dt1$total_stem_height <-as.numeric(dt1$total_stem_height)
if (class(dt1$total_reproductive_heads)=="factor") dt1$total_reproductive_heads <-as.numeric(levels(dt1$total_reproductive_heads))[as.integer(dt1$total_reproductive_heads) ]               
if (class(dt1$total_reproductive_heads)=="character") dt1$total_reproductive_heads <-as.numeric(dt1$total_reproductive_heads)
if (class(dt1$damaged_heads)=="factor") dt1$damaged_heads <-as.numeric(levels(dt1$damaged_heads))[as.integer(dt1$damaged_heads) ]               
if (class(dt1$damaged_heads)=="character") dt1$damaged_heads <-as.numeric(dt1$damaged_heads)
if (class(dt1$rosettes)=="factor") dt1$rosettes <-as.numeric(levels(dt1$rosettes))[as.integer(dt1$rosettes) ]               
if (class(dt1$rosettes)=="character") dt1$rosettes <-as.numeric(dt1$rosettes)
if (class(dt1$leaves)=="factor") dt1$leaves <-as.numeric(levels(dt1$leaves))[as.integer(dt1$leaves) ]               
if (class(dt1$leaves)=="character") dt1$leaves <-as.numeric(dt1$leaves)
if (class(dt1$comments)!="factor") dt1$comments<- as.factor(dt1$comments)
if (class(dt1$QAQCcomments)!="factor") dt1$QAQCcomments<- as.factor(dt1$QAQCcomments)

# Convert Missing Values to NA for non-dates

dt1$toothpick <- as.factor(ifelse((trimws(as.character(dt1$toothpick))==trimws("NA")),NA,as.character(dt1$toothpick)))
dt1$hurricane <- as.factor(ifelse((trimws(as.character(dt1$hurricane))==trimws("NA")),NA,as.character(dt1$hurricane)))
dt1$burn <- as.factor(ifelse((trimws(as.character(dt1$burn))==trimws("NA")),NA,as.character(dt1$burn)))
dt1$life_stage <- as.factor(ifelse((trimws(as.character(dt1$life_stage))==trimws("NA")),NA,as.character(dt1$life_stage)))
dt1$reproductive_stems <- ifelse((trimws(as.character(dt1$reproductive_stems))==trimws("NA")),NA,dt1$reproductive_stems)               
suppressWarnings(dt1$reproductive_stems <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$reproductive_stems))==as.character(as.numeric("NA"))),NA,dt1$reproductive_stems))
dt1$topped_stems <- ifelse((trimws(as.character(dt1$topped_stems))==trimws("NA")),NA,dt1$topped_stems)               
suppressWarnings(dt1$topped_stems <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$topped_stems))==as.character(as.numeric("NA"))),NA,dt1$topped_stems))
dt1$total_stem_height <- ifelse((trimws(as.character(dt1$total_stem_height))==trimws("NA")),NA,dt1$total_stem_height)               
suppressWarnings(dt1$total_stem_height <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$total_stem_height))==as.character(as.numeric("NA"))),NA,dt1$total_stem_height))
dt1$total_reproductive_heads <- ifelse((trimws(as.character(dt1$total_reproductive_heads))==trimws("NA")),NA,dt1$total_reproductive_heads)               
suppressWarnings(dt1$total_reproductive_heads <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$total_reproductive_heads))==as.character(as.numeric("NA"))),NA,dt1$total_reproductive_heads))
dt1$damaged_heads <- ifelse((trimws(as.character(dt1$damaged_heads))==trimws("NA")),NA,dt1$damaged_heads)               
suppressWarnings(dt1$damaged_heads <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$damaged_heads))==as.character(as.numeric("NA"))),NA,dt1$damaged_heads))
dt1$rosettes <- ifelse((trimws(as.character(dt1$rosettes))==trimws("NA")),NA,dt1$rosettes)               
suppressWarnings(dt1$rosettes <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$rosettes))==as.character(as.numeric("NA"))),NA,dt1$rosettes))
dt1$leaves <- ifelse((trimws(as.character(dt1$leaves))==trimws("NA")),NA,dt1$leaves)               
suppressWarnings(dt1$leaves <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$leaves))==as.character(as.numeric("NA"))),NA,dt1$leaves))
dt1$comments <- as.factor(ifelse((trimws(as.character(dt1$comments))==trimws("NA")),NA,as.character(dt1$comments)))
dt1$QAQCcomments <- as.factor(ifelse((trimws(as.character(dt1$QAQCcomments))==trimws("NA")),NA,as.character(dt1$QAQCcomments)))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(site)
summary(population)
summary(toothpick)
summary(tag)
summary(year)
summary(hurricane)
summary(habitat)
summary(burn)
summary(survival)
summary(life_stage)
summary(reproductive_stems)
summary(topped_stems)
summary(total_stem_height)
summary(total_reproductive_heads)
summary(damaged_heads)
summary(rosettes)
summary(leaves)
summary(comments)
summary(QAQCcomments) 
# Get more details on character variables

summary(as.factor(dt1$site)) 
summary(as.factor(dt1$population)) 
summary(as.factor(dt1$toothpick)) 
summary(as.factor(dt1$hurricane)) 
summary(as.factor(dt1$habitat)) 
summary(as.factor(dt1$burn)) 
summary(as.factor(dt1$survival)) 
summary(as.factor(dt1$life_stage)) 
summary(as.factor(dt1$comments)) 
summary(as.factor(dt1$QAQCcomments))
detach(dt1)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('archbold')
# Define species
v_species <- c('Liatris ohlingerae')
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
           paste0(dir_data, '/liatris_ohlingerae_data.csv'), 
           row.names = F )
