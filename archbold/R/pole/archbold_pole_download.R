# Package ID: edi.318.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Demography of the rare herb, Polygala lewtonii (Polygalaceae), with multiple fires on the Lake Wales Ridge in south-central Florida from 2001-2017 (ongoing).
# Data set creator:  Eric Menges - Archbold Biological Station 
# Data set creator:  Carl Weekley - Archbold Biological Station 
# Data set creator:  Stephanie Koontz - Archbold Biological Station 
# Contact:    - Data Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       
options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/318/1/7e1203f9394667d14c3bf2c4a26910f3" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "date",     
                 "unit",     
                 "quad",     
                 "cohort",     
                 "id",     
                 "angle",     
                 "dist",     
                 "survival",     
                 "stage",     
                 "height",     
                 "crown_diameter",     
                 "stems",     
                 "flowering_stems",     
                 "postburn_plant"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$unit)=="factor") dt1$unit <-as.numeric(levels(dt1$unit))[as.integer(dt1$unit) ]               
if (class(dt1$unit)=="character") dt1$unit <-as.numeric(dt1$unit)
if (class(dt1$quad)=="factor") dt1$quad <-as.numeric(levels(dt1$quad))[as.integer(dt1$quad) ]               
if (class(dt1$quad)=="character") dt1$quad <-as.numeric(dt1$quad)
if (class(dt1$cohort)!="factor") dt1$cohort<- as.factor(dt1$cohort)
if (class(dt1$id)=="factor") dt1$id <-as.numeric(levels(dt1$id))[as.integer(dt1$id) ]               
if (class(dt1$id)=="character") dt1$id <-as.numeric(dt1$id)
if (class(dt1$angle)=="factor") dt1$angle <-as.numeric(levels(dt1$angle))[as.integer(dt1$angle) ]               
if (class(dt1$angle)=="character") dt1$angle <-as.numeric(dt1$angle)
if (class(dt1$dist)=="factor") dt1$dist <-as.numeric(levels(dt1$dist))[as.integer(dt1$dist) ]               
if (class(dt1$dist)=="character") dt1$dist <-as.numeric(dt1$dist)
if (class(dt1$survival)!="factor") dt1$survival<- as.factor(dt1$survival)
if (class(dt1$stage)!="factor") dt1$stage<- as.factor(dt1$stage)
if (class(dt1$height)=="factor") dt1$height <-as.numeric(levels(dt1$height))[as.integer(dt1$height) ]               
if (class(dt1$height)=="character") dt1$height <-as.numeric(dt1$height)
if (class(dt1$crown_diameter)=="factor") dt1$crown_diameter <-as.numeric(levels(dt1$crown_diameter))[as.integer(dt1$crown_diameter) ]               
if (class(dt1$crown_diameter)=="character") dt1$crown_diameter <-as.numeric(dt1$crown_diameter)
if (class(dt1$stems)=="factor") dt1$stems <-as.numeric(levels(dt1$stems))[as.integer(dt1$stems) ]               
if (class(dt1$stems)=="character") dt1$stems <-as.numeric(dt1$stems)
if (class(dt1$flowering_stems)=="factor") dt1$flowering_stems <-as.numeric(levels(dt1$flowering_stems))[as.integer(dt1$flowering_stems) ]               
if (class(dt1$flowering_stems)=="character") dt1$flowering_stems <-as.numeric(dt1$flowering_stems)
if (class(dt1$postburn_plant)!="factor") dt1$postburn_plant<- as.factor(dt1$postburn_plant)

# Convert Missing Values to NA for non-dates

dt1$angle <- ifelse((trimws(as.character(dt1$angle))==trimws("NA")),NA,dt1$angle)               
suppressWarnings(dt1$angle <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$angle))==as.character(as.numeric("NA"))),NA,dt1$angle))
dt1$dist <- ifelse((trimws(as.character(dt1$dist))==trimws("NA")),NA,dt1$dist)               
suppressWarnings(dt1$dist <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$dist))==as.character(as.numeric("NA"))),NA,dt1$dist))
dt1$stage <- as.factor(ifelse((trimws(as.character(dt1$stage))==trimws("NA")),NA,as.character(dt1$stage)))
dt1$height <- ifelse((trimws(as.character(dt1$height))==trimws("NA")),NA,dt1$height)               
suppressWarnings(dt1$height <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$height))==as.character(as.numeric("NA"))),NA,dt1$height))
dt1$crown_diameter <- ifelse((trimws(as.character(dt1$crown_diameter))==trimws("NA")),NA,dt1$crown_diameter)               
suppressWarnings(dt1$crown_diameter <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$crown_diameter))==as.character(as.numeric("NA"))),NA,dt1$crown_diameter))
dt1$stems <- ifelse((trimws(as.character(dt1$stems))==trimws("NA")),NA,dt1$stems)               
suppressWarnings(dt1$stems <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$stems))==as.character(as.numeric("NA"))),NA,dt1$stems))
dt1$flowering_stems <- ifelse((trimws(as.character(dt1$flowering_stems))==trimws("NA")),NA,dt1$flowering_stems)               
suppressWarnings(dt1$flowering_stems <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$flowering_stems))==as.character(as.numeric("NA"))),NA,dt1$flowering_stems))
dt1$postburn_plant <- as.factor(ifelse((trimws(as.character(dt1$postburn_plant))==trimws("NA")),NA,as.character(dt1$postburn_plant)))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(date)
summary(unit)
summary(quad)
summary(cohort)
summary(id)
summary(angle)
summary(dist)
summary(survival)
summary(stage)
summary(height)
summary(crown_diameter)
summary(stems)
summary(flowering_stems)
summary(postburn_plant) 
# Get more details on character variables

summary(as.factor(dt1$cohort)) 
summary(as.factor(dt1$survival)) 
summary(as.factor(dt1$stage)) 
summary(as.factor(dt1$postburn_plant))
detach(dt1)               



inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/318/1/e15db9d220153075e7a565f6344b19b7" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "date",     
                 "unit",     
                 "quad",     
                 "pct_burn",     
                 "burn_status"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$unit)=="factor") dt2$unit <-as.numeric(levels(dt2$unit))[as.integer(dt2$unit) ]               
if (class(dt2$unit)=="character") dt2$unit <-as.numeric(dt2$unit)
if (class(dt2$quad)=="factor") dt2$quad <-as.numeric(levels(dt2$quad))[as.integer(dt2$quad) ]               
if (class(dt2$quad)=="character") dt2$quad <-as.numeric(dt2$quad)
if (class(dt2$pct_burn)!="factor") dt2$pct_burn<- as.factor(dt2$pct_burn)
if (class(dt2$burn_status)!="factor") dt2$burn_status<- as.factor(dt2$burn_status)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt2)                            
attach(dt2)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(date)
summary(unit)
summary(quad)
summary(pct_burn)
summary(burn_status) 
# Get more details on character variables

summary(as.factor(dt2$pct_burn)) 
summary(as.factor(dt2$burn_status))
detach(dt2)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('archbold')
# Define species
v_species <- c('Polygala lewtonii')
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
           paste0(dir_data, '/polygala_lewtonii_data.csv'), 
           row.names = F )
