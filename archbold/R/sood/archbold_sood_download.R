# Package ID: edi.246.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Dynamics of Chapman's goldenrod (Solidago odora var. chapmanii) with fire at Archbold Biological Station, 1991-2000.
# Data set creator:  Eric Menges - Archbold Biological Station 
# Data set creator:  Richard Root - Cornell University 
# Contact:    - Data Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       
options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/246/1/fef4b6d96f26a2a91735533a9aabf1d2" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "quad",     
                 "fire90",     
                 "fire91",     
                 "fire92",     
                 "fire95",     
                 "fire98",     
                 "nfires",     
                 "year",     
                 "st",     
                 "fl",     
                 "mht"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$quad)=="factor") dt1$quad <-as.numeric(levels(dt1$quad))[as.integer(dt1$quad) ]               
if (class(dt1$quad)=="character") dt1$quad <-as.numeric(dt1$quad)
if (class(dt1$fire90)!="factor") dt1$fire90<- as.factor(dt1$fire90)
if (class(dt1$fire91)!="factor") dt1$fire91<- as.factor(dt1$fire91)
if (class(dt1$fire92)!="factor") dt1$fire92<- as.factor(dt1$fire92)
if (class(dt1$fire95)!="factor") dt1$fire95<- as.factor(dt1$fire95)
if (class(dt1$fire98)!="factor") dt1$fire98<- as.factor(dt1$fire98)
if (class(dt1$nfires)=="factor") dt1$nfires <-as.numeric(levels(dt1$nfires))[as.integer(dt1$nfires) ]               
if (class(dt1$nfires)=="character") dt1$nfires <-as.numeric(dt1$nfires)
if (class(dt1$st)=="factor") dt1$st <-as.numeric(levels(dt1$st))[as.integer(dt1$st) ]               
if (class(dt1$st)=="character") dt1$st <-as.numeric(dt1$st)
if (class(dt1$fl)=="factor") dt1$fl <-as.numeric(levels(dt1$fl))[as.integer(dt1$fl) ]               
if (class(dt1$fl)=="character") dt1$fl <-as.numeric(dt1$fl)
if (class(dt1$mht)=="factor") dt1$mht <-as.numeric(levels(dt1$mht))[as.integer(dt1$mht) ]               
if (class(dt1$mht)=="character") dt1$mht <-as.numeric(dt1$mht)

# Convert Missing Values to NA for non-dates

dt1$st <- ifelse((trimws(as.character(dt1$st))==trimws("NA")),NA,dt1$st)               
suppressWarnings(dt1$st <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$st))==as.character(as.numeric("NA"))),NA,dt1$st))
dt1$fl <- ifelse((trimws(as.character(dt1$fl))==trimws("NA")),NA,dt1$fl)               
suppressWarnings(dt1$fl <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$fl))==as.character(as.numeric("NA"))),NA,dt1$fl))
dt1$mht <- ifelse((trimws(as.character(dt1$mht))==trimws("NA")),NA,dt1$mht)               
suppressWarnings(dt1$mht <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$mht))==as.character(as.numeric("NA"))),NA,dt1$mht))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(quad)
summary(fire90)
summary(fire91)
summary(fire92)
summary(fire95)
summary(fire98)
summary(nfires)
summary(year)
summary(st)
summary(fl)
summary(mht) 
# Get more details on character variables

summary(as.factor(dt1$fire90)) 
summary(as.factor(dt1$fire91)) 
summary(as.factor(dt1$fire92)) 
summary(as.factor(dt1$fire95)) 
summary(as.factor(dt1$fire98))
detach(dt1)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('archbold')
# Define species
v_species <- c('Solidago odora')
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
           paste0(dir_data, '/solidago_odora_data.csv'), 
           row.names = F )
