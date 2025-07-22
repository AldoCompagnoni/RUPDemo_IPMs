# Package ID: edi.181.2 Cataloging System:https://pasta.edirepository.org.
# Data set title: Demographic measures of Hypericum cumulicola (Hypericaceae) in 15 populations in Florida Rosemary Scrub patches with different time-since-fire, at Archbold Biological Station, Highlands County, Florida from 1994-2015.
# Data set creator:  Pedro Quintana-Ascencio - University of Central Florida 
# Data set creator:  Eric Menges - Archbold Biological Station 
# Contact:    - Database Manager Archbold Biological Station  - datamanager@archbold-station.org
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       
options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/181/2/adf4ed268bdbae50f4045549a688ef09" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <- read.csv(infile1, header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Year",     
                 "Site",     
                 "Gap",     
                 "Tag",     
                 "Burn_yr",     
                 "Time",     
                 "Surv_init",     
                 "init_height",     
                 "Rep_init",     
                 "Stems",     
                 "Surv_fin",     
                 "Fin_height",     
                 "Rep_fin"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$Site)=="factor") dt1$Site <-as.numeric(levels(dt1$Site))[as.integer(dt1$Site) ]               
if (class(dt1$Site)=="character") dt1$Site <-as.numeric(dt1$Site)
if (class(dt1$Gap)=="factor") dt1$Gap <-as.numeric(levels(dt1$Gap))[as.integer(dt1$Gap) ]               
if (class(dt1$Gap)=="character") dt1$Gap <-as.numeric(dt1$Gap)
if (class(dt1$Tag)=="factor") dt1$Tag <-as.numeric(levels(dt1$Tag))[as.integer(dt1$Tag) ]               
if (class(dt1$Tag)=="character") dt1$Tag <-as.numeric(dt1$Tag)
if (class(dt1$Surv_init)!="factor") dt1$Surv_init<- as.factor(dt1$Surv_init)
if (class(dt1$init_height)=="factor") dt1$init_height <-as.numeric(levels(dt1$init_height))[as.integer(dt1$init_height) ]               
if (class(dt1$init_height)=="character") dt1$init_height <-as.numeric(dt1$init_height)
if (class(dt1$Rep_init)=="factor") dt1$Rep_init <-as.numeric(levels(dt1$Rep_init))[as.integer(dt1$Rep_init) ]               
if (class(dt1$Rep_init)=="character") dt1$Rep_init <-as.numeric(dt1$Rep_init)
if (class(dt1$Stems)=="factor") dt1$Stems <-as.numeric(levels(dt1$Stems))[as.integer(dt1$Stems) ]               
if (class(dt1$Stems)=="character") dt1$Stems <-as.numeric(dt1$Stems)
if (class(dt1$Surv_fin)!="factor") dt1$Surv_fin<- as.factor(dt1$Surv_fin)
if (class(dt1$Fin_height)=="factor") dt1$Fin_height <-as.numeric(levels(dt1$Fin_height))[as.integer(dt1$Fin_height) ]               
if (class(dt1$Fin_height)=="character") dt1$Fin_height <-as.numeric(dt1$Fin_height)
if (class(dt1$Rep_fin)=="factor") dt1$Rep_fin <-as.numeric(levels(dt1$Rep_fin))[as.integer(dt1$Rep_fin) ]               
if (class(dt1$Rep_fin)=="character") dt1$Rep_fin <-as.numeric(dt1$Rep_fin)

# Convert Missing Values to NA for non-dates

dt1$Site <- ifelse((trimws(as.character(dt1$Site))==trimws("NA")),NA,dt1$Site)               
suppressWarnings(dt1$Site <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Site))==as.character(as.numeric("NA"))),NA,dt1$Site))
dt1$Gap <- ifelse((trimws(as.character(dt1$Gap))==trimws("NA")),NA,dt1$Gap)               
suppressWarnings(dt1$Gap <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Gap))==as.character(as.numeric("NA"))),NA,dt1$Gap))
dt1$Tag <- ifelse((trimws(as.character(dt1$Tag))==trimws("NA")),NA,dt1$Tag)               
suppressWarnings(dt1$Tag <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Tag))==as.character(as.numeric("NA"))),NA,dt1$Tag))
dt1$Surv_init <- as.factor(ifelse((trimws(as.character(dt1$Surv_init))==trimws("NA")),NA,as.character(dt1$Surv_init)))
dt1$init_height <- ifelse((trimws(as.character(dt1$init_height))==trimws("NA")),NA,dt1$init_height)               
suppressWarnings(dt1$init_height <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$init_height))==as.character(as.numeric("NA"))),NA,dt1$init_height))
dt1$Rep_init <- ifelse((trimws(as.character(dt1$Rep_init))==trimws("NA")),NA,dt1$Rep_init)               
suppressWarnings(dt1$Rep_init <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Rep_init))==as.character(as.numeric("NA"))),NA,dt1$Rep_init))
dt1$Stems <- ifelse((trimws(as.character(dt1$Stems))==trimws("NA")),NA,dt1$Stems)               
suppressWarnings(dt1$Stems <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Stems))==as.character(as.numeric("NA"))),NA,dt1$Stems))
dt1$Surv_fin <- as.factor(ifelse((trimws(as.character(dt1$Surv_fin))==trimws("NA")),NA,as.character(dt1$Surv_fin)))
dt1$Fin_height <- ifelse((trimws(as.character(dt1$Fin_height))==trimws("NA")),NA,dt1$Fin_height)               
suppressWarnings(dt1$Fin_height <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Fin_height))==as.character(as.numeric("NA"))),NA,dt1$Fin_height))
dt1$Rep_fin <- ifelse((trimws(as.character(dt1$Rep_fin))==trimws("NA")),NA,dt1$Rep_fin)               
suppressWarnings(dt1$Rep_fin <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Rep_fin))==as.character(as.numeric("NA"))),NA,dt1$Rep_fin))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(Year)
summary(Site)
summary(Gap)
summary(Tag)
summary(Burn_yr)
summary(Time)
summary(Surv_init)
summary(init_height)
summary(Rep_init)
summary(Stems)
summary(Surv_fin)
summary(Fin_height)
summary(Rep_fin) 
# Get more details on character variables

summary(as.factor(dt1$Surv_init)) 
summary(as.factor(dt1$Surv_fin))
detach(dt1)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('archbold')
# Define species
v_species <- c('Hypericum cumulicola')
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
           paste0(dir_data, '/hypericum_cumulicola_data.csv'), 
           row.names = F )
