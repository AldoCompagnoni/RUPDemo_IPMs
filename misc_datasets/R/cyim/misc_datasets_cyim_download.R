# Package ID: knb-lter-sev.323.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Long-term study of tree cholla demography in the Los Pinos mountains, Sevilleta National Wildlife Refuge.
# Data set creator:  Tom Miller - Rice University 
# Contact:    - Information Manager University of New Mexico  - khall001@unm.edu
# Stylesheet v2.15 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       
options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sev/323/1/e3a463a3bc1c2dc1b84f9581ecccd6fb" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Plot",     
                 "TagID",     
                 "Transplant",     
                 "Year_t",     
                 "Height_t",     
                 "Width_t",     
                 "Perp_t",     
                 "NS_t",     
                 "Goodbuds_t",     
                 "TotFlowerbuds_t",     
                 "ABFlowerbuds_t",     
                 "Antcount_t",     
                 "Ant_sp_t",     
                 "Year_t1",     
                 "Newplant",     
                 "Survival_t1",     
                 "Height_t1",     
                 "Width_t1",     
                 "Perp_t1",     
                 "NS_t1",     
                 "Goodbuds_t1",     
                 "TotFlowerbuds_t1",     
                 "ABFlowerbuds_t1",     
                 "Antcount_t1",     
                 "Ant_sp_t1",     
                 "comments"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$Plot)!="factor") dt1$Plot<- as.factor(dt1$Plot)
if (class(dt1$TagID)!="factor") dt1$TagID<- as.factor(dt1$TagID)
if (class(dt1$Transplant)!="factor") dt1$Transplant<- as.factor(dt1$Transplant)
if (class(dt1$Height_t)=="factor") dt1$Height_t <-as.numeric(levels(dt1$Height_t))[as.integer(dt1$Height_t) ]               
if (class(dt1$Height_t)=="character") dt1$Height_t <-as.numeric(dt1$Height_t)
if (class(dt1$Width_t)=="factor") dt1$Width_t <-as.numeric(levels(dt1$Width_t))[as.integer(dt1$Width_t) ]               
if (class(dt1$Width_t)=="character") dt1$Width_t <-as.numeric(dt1$Width_t)
if (class(dt1$Perp_t)=="factor") dt1$Perp_t <-as.numeric(levels(dt1$Perp_t))[as.integer(dt1$Perp_t) ]               
if (class(dt1$Perp_t)=="character") dt1$Perp_t <-as.numeric(dt1$Perp_t)
if (class(dt1$NS_t)=="factor") dt1$NS_t <-as.numeric(levels(dt1$NS_t))[as.integer(dt1$NS_t) ]               
if (class(dt1$NS_t)=="character") dt1$NS_t <-as.numeric(dt1$NS_t)
if (class(dt1$Goodbuds_t)=="factor") dt1$Goodbuds_t <-as.numeric(levels(dt1$Goodbuds_t))[as.integer(dt1$Goodbuds_t) ]               
if (class(dt1$Goodbuds_t)=="character") dt1$Goodbuds_t <-as.numeric(dt1$Goodbuds_t)
if (class(dt1$TotFlowerbuds_t)=="factor") dt1$TotFlowerbuds_t <-as.numeric(levels(dt1$TotFlowerbuds_t))[as.integer(dt1$TotFlowerbuds_t) ]               
if (class(dt1$TotFlowerbuds_t)=="character") dt1$TotFlowerbuds_t <-as.numeric(dt1$TotFlowerbuds_t)
if (class(dt1$ABFlowerbuds_t)=="factor") dt1$ABFlowerbuds_t <-as.numeric(levels(dt1$ABFlowerbuds_t))[as.integer(dt1$ABFlowerbuds_t) ]               
if (class(dt1$ABFlowerbuds_t)=="character") dt1$ABFlowerbuds_t <-as.numeric(dt1$ABFlowerbuds_t)
if (class(dt1$Antcount_t)=="factor") dt1$Antcount_t <-as.numeric(levels(dt1$Antcount_t))[as.integer(dt1$Antcount_t) ]               
if (class(dt1$Antcount_t)=="character") dt1$Antcount_t <-as.numeric(dt1$Antcount_t)
if (class(dt1$Ant_sp_t)!="factor") dt1$Ant_sp_t<- as.factor(dt1$Ant_sp_t)
if (class(dt1$Newplant)!="factor") dt1$Newplant<- as.factor(dt1$Newplant)
if (class(dt1$Survival_t1)!="factor") dt1$Survival_t1<- as.factor(dt1$Survival_t1)
if (class(dt1$Height_t1)=="factor") dt1$Height_t1 <-as.numeric(levels(dt1$Height_t1))[as.integer(dt1$Height_t1) ]               
if (class(dt1$Height_t1)=="character") dt1$Height_t1 <-as.numeric(dt1$Height_t1)
if (class(dt1$Width_t1)=="factor") dt1$Width_t1 <-as.numeric(levels(dt1$Width_t1))[as.integer(dt1$Width_t1) ]               
if (class(dt1$Width_t1)=="character") dt1$Width_t1 <-as.numeric(dt1$Width_t1)
if (class(dt1$Perp_t1)=="factor") dt1$Perp_t1 <-as.numeric(levels(dt1$Perp_t1))[as.integer(dt1$Perp_t1) ]               
if (class(dt1$Perp_t1)=="character") dt1$Perp_t1 <-as.numeric(dt1$Perp_t1)
if (class(dt1$NS_t1)=="factor") dt1$NS_t1 <-as.numeric(levels(dt1$NS_t1))[as.integer(dt1$NS_t1) ]               
if (class(dt1$NS_t1)=="character") dt1$NS_t1 <-as.numeric(dt1$NS_t1)
if (class(dt1$Goodbuds_t1)=="factor") dt1$Goodbuds_t1 <-as.numeric(levels(dt1$Goodbuds_t1))[as.integer(dt1$Goodbuds_t1) ]               
if (class(dt1$Goodbuds_t1)=="character") dt1$Goodbuds_t1 <-as.numeric(dt1$Goodbuds_t1)
if (class(dt1$TotFlowerbuds_t1)=="factor") dt1$TotFlowerbuds_t1 <-as.numeric(levels(dt1$TotFlowerbuds_t1))[as.integer(dt1$TotFlowerbuds_t1) ]               
if (class(dt1$TotFlowerbuds_t1)=="character") dt1$TotFlowerbuds_t1 <-as.numeric(dt1$TotFlowerbuds_t1)
if (class(dt1$ABFlowerbuds_t1)=="factor") dt1$ABFlowerbuds_t1 <-as.numeric(levels(dt1$ABFlowerbuds_t1))[as.integer(dt1$ABFlowerbuds_t1) ]               
if (class(dt1$ABFlowerbuds_t1)=="character") dt1$ABFlowerbuds_t1 <-as.numeric(dt1$ABFlowerbuds_t1)
if (class(dt1$Antcount_t1)=="factor") dt1$Antcount_t1 <-as.numeric(levels(dt1$Antcount_t1))[as.integer(dt1$Antcount_t1) ]               
if (class(dt1$Antcount_t1)=="character") dt1$Antcount_t1 <-as.numeric(dt1$Antcount_t1)
if (class(dt1$Ant_sp_t1)!="factor") dt1$Ant_sp_t1<- as.factor(dt1$Ant_sp_t1)
if (class(dt1$comments)!="factor") dt1$comments<- as.factor(dt1$comments)

# Convert Missing Values to NA for non-dates

dt1$Plot <- as.factor(ifelse((trimws(as.character(dt1$Plot))==trimws("NA")),NA,as.character(dt1$Plot)))
dt1$TagID <- as.factor(ifelse((trimws(as.character(dt1$TagID))==trimws("NA")),NA,as.character(dt1$TagID)))
dt1$Transplant <- as.factor(ifelse((trimws(as.character(dt1$Transplant))==trimws("NA")),NA,as.character(dt1$Transplant)))
dt1$Height_t <- ifelse((trimws(as.character(dt1$Height_t))==trimws("NA")),NA,dt1$Height_t)               
suppressWarnings(dt1$Height_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Height_t))==as.character(as.numeric("NA"))),NA,dt1$Height_t))
dt1$Width_t <- ifelse((trimws(as.character(dt1$Width_t))==trimws("NA")),NA,dt1$Width_t)               
suppressWarnings(dt1$Width_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Width_t))==as.character(as.numeric("NA"))),NA,dt1$Width_t))
dt1$Perp_t <- ifelse((trimws(as.character(dt1$Perp_t))==trimws("NA")),NA,dt1$Perp_t)               
suppressWarnings(dt1$Perp_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Perp_t))==as.character(as.numeric("NA"))),NA,dt1$Perp_t))
dt1$NS_t <- ifelse((trimws(as.character(dt1$NS_t))==trimws("NA")),NA,dt1$NS_t)               
suppressWarnings(dt1$NS_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$NS_t))==as.character(as.numeric("NA"))),NA,dt1$NS_t))
dt1$Goodbuds_t <- ifelse((trimws(as.character(dt1$Goodbuds_t))==trimws("NA")),NA,dt1$Goodbuds_t)               
suppressWarnings(dt1$Goodbuds_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Goodbuds_t))==as.character(as.numeric("NA"))),NA,dt1$Goodbuds_t))
dt1$TotFlowerbuds_t <- ifelse((trimws(as.character(dt1$TotFlowerbuds_t))==trimws("NA")),NA,dt1$TotFlowerbuds_t)               
suppressWarnings(dt1$TotFlowerbuds_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotFlowerbuds_t))==as.character(as.numeric("NA"))),NA,dt1$TotFlowerbuds_t))
dt1$ABFlowerbuds_t <- ifelse((trimws(as.character(dt1$ABFlowerbuds_t))==trimws("NA")),NA,dt1$ABFlowerbuds_t)               
suppressWarnings(dt1$ABFlowerbuds_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$ABFlowerbuds_t))==as.character(as.numeric("NA"))),NA,dt1$ABFlowerbuds_t))
dt1$Antcount_t <- ifelse((trimws(as.character(dt1$Antcount_t))==trimws("NA")),NA,dt1$Antcount_t)               
suppressWarnings(dt1$Antcount_t <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Antcount_t))==as.character(as.numeric("NA"))),NA,dt1$Antcount_t))
dt1$Ant_sp_t <- as.factor(ifelse((trimws(as.character(dt1$Ant_sp_t))==trimws("NA")),NA,as.character(dt1$Ant_sp_t)))
dt1$Newplant <- as.factor(ifelse((trimws(as.character(dt1$Newplant))==trimws("NA")),NA,as.character(dt1$Newplant)))
dt1$Survival_t1 <- as.factor(ifelse((trimws(as.character(dt1$Survival_t1))==trimws("NA")),NA,as.character(dt1$Survival_t1)))
dt1$Height_t1 <- ifelse((trimws(as.character(dt1$Height_t1))==trimws("NA")),NA,dt1$Height_t1)               
suppressWarnings(dt1$Height_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Height_t1))==as.character(as.numeric("NA"))),NA,dt1$Height_t1))
dt1$Width_t1 <- ifelse((trimws(as.character(dt1$Width_t1))==trimws("NA")),NA,dt1$Width_t1)               
suppressWarnings(dt1$Width_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Width_t1))==as.character(as.numeric("NA"))),NA,dt1$Width_t1))
dt1$Perp_t1 <- ifelse((trimws(as.character(dt1$Perp_t1))==trimws("NA")),NA,dt1$Perp_t1)               
suppressWarnings(dt1$Perp_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Perp_t1))==as.character(as.numeric("NA"))),NA,dt1$Perp_t1))
dt1$NS_t1 <- ifelse((trimws(as.character(dt1$NS_t1))==trimws("NA")),NA,dt1$NS_t1)               
suppressWarnings(dt1$NS_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$NS_t1))==as.character(as.numeric("NA"))),NA,dt1$NS_t1))
dt1$Goodbuds_t1 <- ifelse((trimws(as.character(dt1$Goodbuds_t1))==trimws("NA")),NA,dt1$Goodbuds_t1)               
suppressWarnings(dt1$Goodbuds_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Goodbuds_t1))==as.character(as.numeric("NA"))),NA,dt1$Goodbuds_t1))
dt1$TotFlowerbuds_t1 <- ifelse((trimws(as.character(dt1$TotFlowerbuds_t1))==trimws("NA")),NA,dt1$TotFlowerbuds_t1)               
suppressWarnings(dt1$TotFlowerbuds_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotFlowerbuds_t1))==as.character(as.numeric("NA"))),NA,dt1$TotFlowerbuds_t1))
dt1$ABFlowerbuds_t1 <- ifelse((trimws(as.character(dt1$ABFlowerbuds_t1))==trimws("NA")),NA,dt1$ABFlowerbuds_t1)               
suppressWarnings(dt1$ABFlowerbuds_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$ABFlowerbuds_t1))==as.character(as.numeric("NA"))),NA,dt1$ABFlowerbuds_t1))
dt1$Antcount_t1 <- ifelse((trimws(as.character(dt1$Antcount_t1))==trimws("NA")),NA,dt1$Antcount_t1)               
suppressWarnings(dt1$Antcount_t1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Antcount_t1))==as.character(as.numeric("NA"))),NA,dt1$Antcount_t1))
dt1$Ant_sp_t1 <- as.factor(ifelse((trimws(as.character(dt1$Ant_sp_t1))==trimws("NA")),NA,as.character(dt1$Ant_sp_t1)))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(Plot)
summary(TagID)
summary(Transplant)
summary(Year_t)
summary(Height_t)
summary(Width_t)
summary(Perp_t)
summary(NS_t)
summary(Goodbuds_t)
summary(TotFlowerbuds_t)
summary(ABFlowerbuds_t)
summary(Antcount_t)
summary(Ant_sp_t)
summary(Year_t1)
summary(Newplant)
summary(Survival_t1)
summary(Height_t1)
summary(Width_t1)
summary(Perp_t1)
summary(NS_t1)
summary(Goodbuds_t1)
summary(TotFlowerbuds_t1)
summary(ABFlowerbuds_t1)
summary(Antcount_t1)
summary(Ant_sp_t1)
summary(comments) 
# Get more details on character variables

summary(as.factor(dt1$Plot)) 
summary(as.factor(dt1$TagID)) 
summary(as.factor(dt1$Transplant)) 
summary(as.factor(dt1$Ant_sp_t)) 
summary(as.factor(dt1$Newplant)) 
summary(as.factor(dt1$Survival_t1)) 
summary(as.factor(dt1$Ant_sp_t1)) 
summary(as.factor(dt1$comments))
detach(dt1)               



# store in our local repository ------------------------------------------------
library(tidyverse)

# Store Crotalaria avonensis
# Define head-directory 
v_head    <- c('misc_datasets')
# Define species
v_species <- c('Cylindriopuntia imbricata')
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
           paste0(dir_data, '/cylindriopuntia_imbricata_data.csv'), 
           row.names = F )
