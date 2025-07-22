# EDI code to download the Eriogonum longifolium var. gnaphalifolium dataset

# Package ID: edi.226.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Demography of the Florida endemic Eriogonum longifolium var. gnaphalifolium (Polygonaceae) at Archbold Biological Station and the Lake Wales Ridge National Wildlife Refuge Carter Creek from 1990 to 2013.
# Data set creator:  Eric Menges - Archbold Biological Station 
# Contact:    - Data Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.14 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.226.1
# Study organism: scrub buckwheat (Eriogonum longifolium var. gnaphalifolium)
# Time period: 1990 - 2013 (Archbold); 2001 - 2010 (Lake Wales Ridge National Wildlife Refuge Carter Creek)
# Ecology: long-lived; dormancy; 
#  resprouts following fire and also recruits many seedlings following post-fire stimulated flowering (minimum of every five years or every 20 years are necessary fpr viablility);


# setwd("C:/users/my_name/my_dir")       

options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/226/1/17019e4dd573d9b8253bb7ce3a6b5e3a" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "date",     
                 "site",     
                 "pop",     
                 "qu",     
                 "plant",     
                 "s",     
                 "stage",     
                 "dia",     
                 "scape",     
                 "invol",     
                 "nb",     
                 "nr",     
                 "herb",     
                 "burn",     
                 "comment"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$site)!="factor") dt1$site<- as.factor(dt1$site)
if (class(dt1$pop)!="factor") dt1$pop<- as.factor(dt1$pop)
if (class(dt1$qu)=="factor") dt1$qu <-as.numeric(levels(dt1$qu))[as.integer(dt1$qu) ]               
if (class(dt1$qu)=="character") dt1$qu <-as.numeric(dt1$qu)
if (class(dt1$plant)=="factor") dt1$plant <-as.numeric(levels(dt1$plant))[as.integer(dt1$plant) ]               
if (class(dt1$plant)=="character") dt1$plant <-as.numeric(dt1$plant)
if (class(dt1$s)!="factor") dt1$s<- as.factor(dt1$s)
if (class(dt1$stage)!="factor") dt1$stage<- as.factor(dt1$stage)
if (class(dt1$dia)=="factor") dt1$dia <-as.numeric(levels(dt1$dia))[as.integer(dt1$dia) ]               
if (class(dt1$dia)=="character") dt1$dia <-as.numeric(dt1$dia)
if (class(dt1$scape)=="factor") dt1$scape <-as.numeric(levels(dt1$scape))[as.integer(dt1$scape) ]               
if (class(dt1$scape)=="character") dt1$scape <-as.numeric(dt1$scape)
if (class(dt1$invol)=="factor") dt1$invol <-as.numeric(levels(dt1$invol))[as.integer(dt1$invol) ]               
if (class(dt1$invol)=="character") dt1$invol <-as.numeric(dt1$invol)
if (class(dt1$nb)=="factor") dt1$nb <-as.numeric(levels(dt1$nb))[as.integer(dt1$nb) ]               
if (class(dt1$nb)=="character") dt1$nb <-as.numeric(dt1$nb)
if (class(dt1$nr)=="factor") dt1$nr <-as.numeric(levels(dt1$nr))[as.integer(dt1$nr) ]               
if (class(dt1$nr)=="character") dt1$nr <-as.numeric(dt1$nr)
if (class(dt1$herb)!="factor") dt1$herb<- as.factor(dt1$herb)
if (class(dt1$burn)!="factor") dt1$burn<- as.factor(dt1$burn)
if (class(dt1$comment)!="factor") dt1$comment<- as.factor(dt1$comment)

# Convert Missing Values to NA for non-dates

dt1$qu <- ifelse((trimws(as.character(dt1$qu))==trimws("NA")),NA,dt1$qu)               
suppressWarnings(dt1$qu <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$qu))==as.character(as.numeric("NA"))),NA,dt1$qu))
dt1$s <- as.factor(ifelse((trimws(as.character(dt1$s))==trimws("NA")),NA,as.character(dt1$s)))
dt1$stage <- as.factor(ifelse((trimws(as.character(dt1$stage))==trimws("NA")),NA,as.character(dt1$stage)))
dt1$dia <- ifelse((trimws(as.character(dt1$dia))==trimws("NA")),NA,dt1$dia)               
suppressWarnings(dt1$dia <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$dia))==as.character(as.numeric("NA"))),NA,dt1$dia))
dt1$scape <- ifelse((trimws(as.character(dt1$scape))==trimws("NA")),NA,dt1$scape)               
suppressWarnings(dt1$scape <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$scape))==as.character(as.numeric("NA"))),NA,dt1$scape))
dt1$invol <- ifelse((trimws(as.character(dt1$invol))==trimws("NA")),NA,dt1$invol)               
suppressWarnings(dt1$invol <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$invol))==as.character(as.numeric("NA"))),NA,dt1$invol))
dt1$nb <- ifelse((trimws(as.character(dt1$nb))==trimws("NA")),NA,dt1$nb)               
suppressWarnings(dt1$nb <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$nb))==as.character(as.numeric("NA"))),NA,dt1$nb))
dt1$nr <- ifelse((trimws(as.character(dt1$nr))==trimws("NA")),NA,dt1$nr)               
suppressWarnings(dt1$nr <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$nr))==as.character(as.numeric("NA"))),NA,dt1$nr))
dt1$herb <- as.factor(ifelse((trimws(as.character(dt1$herb))==trimws("NA")),NA,as.character(dt1$herb)))
dt1$burn <- as.factor(ifelse((trimws(as.character(dt1$burn))==trimws("NA")),NA,as.character(dt1$burn)))
dt1$comment <- as.factor(ifelse((trimws(as.character(dt1$comment))==trimws("NA")),NA,as.character(dt1$comment)))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(date)
summary(site)
summary(pop)
summary(qu)
summary(plant)
summary(s)
summary(stage)
summary(dia)
summary(scape)
summary(invol)
summary(nb)
summary(nr)
summary(herb)
summary(burn)
summary(comment) 
# Get more details on character variables

summary(as.factor(dt1$site)) 
summary(as.factor(dt1$pop)) 
summary(as.factor(dt1$s)) 
summary(as.factor(dt1$stage)) 
summary(as.factor(dt1$herb)) 
summary(as.factor(dt1$burn)) 
summary(as.factor(dt1$comment))
detach(dt1)               


# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('archbold')
# Define species
v_species <- c('Eriogonum longifolium')
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
           paste0(dir_data, '/eriogonum_longifolium_data.csv'), 
           row.names = F )
