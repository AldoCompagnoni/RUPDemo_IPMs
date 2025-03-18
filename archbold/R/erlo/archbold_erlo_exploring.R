# Exploring data - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.03.03


# Website    : https://portal.edirepository.org/nis/mapbrowse?packageid=edi.219.1
# Publication: https://bioone.org/journals/southeastern-naturalist/volume-15/issue-3/058.015.0318/Ecology-and-Conservation-of-the-Endangered-Legume-Crotalaria-avonensis-in/10.1656/058.015.0318.short

# rm(list = ls())


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle)


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
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
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

if (!dir.exists(dir_R     )) {dir.create(dir_R     )}
if (!dir.exists(dir_data  )) {dir.create(dir_data  )}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Data -------------------------------------------------------------------------
df <- read_csv(file.path(dir_data, 'eriogonum_longifolium_data.csv')) %>% 
  janitor::clean_names() 

df_meta <- data.frame(variable = colnames(df)) %>% 
  mutate(definition = c(
    'study site', 'quadrat number',	'macroplot number',	
    'plant number (within quad)', 'direction within circular quad',	
    'distance from quad center',	'was quad caged from 2012 onward',	
    'vegetation type', 'year quadrat was initiated', 
    'year-month of observation',	'fire severity Dec. 2014 or Jan. 2015',
    'fire severity Aug. 2005',	'fire severity May/June 2009',
    "fire severity B's ridge 2016", 'fire severity Feb. 2017',
    'fire severity Oct. 2017', 'survival code for month', 'number of stems',
    'number of branch tips', 'number of flowers (corolla showing)', 
    'number of developing fruits', 'number of mature fruits', 'herbivory code',
    'plant identification', 'sample year', 'sample month'))

skimr::skim(df)
  

# Data exploration -------------------------------------------------------------
df %>% 
  separate( date, c('year','month'), sep='-') %>% 
  count( site, year ) %>% 
  ggplot() + 
  geom_point(aes(x = year, y = site)) +
  theme_bw() +
  labs(title    = 'Sampling inventory') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))


# Prepare data frames for analysis ---------------------------------------------

# produce "transition" growth data

# IMPORTANT: ignore month for now: vast majority of sampling done in June!

# data from year t0
df_t0 <- df %>% 
          separate( date, c('year','month'), sep='-') %>% 
          mutate( year = as.numeric(year) ) %>% 
          # rename number of basal rosettes to show what time the data is from
          rename( nb_t0 = nb ) %>% 
          # retain only relevant variables
          dplyr::select(year, site, pop, qu, plant, nb_t0 ) 
df_t1 <- df %>% 
          separate( date, c('year','month'), sep='-') %>% 
          mutate( year = as.numeric(year) ) %>% 
          mutate( year = year + 1 ) %>% 
          # rename number of basal rosettes to show what time the data is from
          rename( nb_t1 = nb ) %>% 
          # retain only relevant variables
          dplyr::select(year, site, pop, qu, plant, nb_t1 ) 
  
# growth model
grow_df <- full_join( df_t0, df_t1)

# simple growth visualization (looking good, play around with
#   log and sqrt transformation!)
ggplot(grow_df ) +
  geom_point( aes( x = nb_t0,
                   y = nb_t1 ) )

