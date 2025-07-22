# Exploring data - Archbold - Eriogonum longifolium

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.22

# Study organism: Polygala lewtonii
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.226.1
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.226.1
# Citing publication: 
# Time periode: 2001-2017
# Plants were censused during their peak of reproduction annually in July and August


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4)


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

# Plot subtitle
v_ggp_suffix    <- paste(
  tools::toTitleCase(v_head), '-', v_species)

# Models
v_mod_set_gr <- c()
v_mod_set_su <- c()
v_mod_set_fl <- c()

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


# Functions --------------------------------------------------------------------
# function to plot your survival data 'binned' (instead of 'jittered')
source('helper_functions/plot_binned_prop.R')


# Data -------------------------------------------------------------------------
df_org <- read_csv(file.path(dir_data, 'eriogonum_longifolium_data.csv')) %>% 
  janitor::clean_names() 

df_meta <- data.frame(variable = colnames(df_org)) %>% 
  mutate(definition = c(
    'year and month data were collected', 'study site',	'population number',	
    'quadrat number', 'unique plant ID', 'annual survival',
    'annual life history stage',	'rosette diameter at its widest point',	
    'number of flowering scapes', 'number of flowering involucre', 
    'number of basal leaves (across all rosettes)',	'number of rosettes',
    'mammalian herbivory', 'plant status following prescribed fire',
    'comments'))

skimr::skim(df_org)


# Tails from meta --------------------------------------------------------------
"annually in June/July"

"Unmarked plants within the plot were marked and coded as either 
new seedlings (basal diameter less than or equal to 2cm) or new adults."

"Extended dormancy exist in this plant therefore flags were never removed from the field."

"We stopped counting leaves in 2000 and stopped counting the number of rosettes in 2008."

"Following a prescribed burn, we scored all previous alive plants as unburned or burned."

"Starting in 2002, we specified if a burned plant was scorched (leaf material present but brown) or consumed (leaf material absent)."

"Fire in 1998 burned in May and July therefore an additional census was done in September 1998"


# Sampling desing --------------------------------------------------------------
df_org %>% 
  separate(date, c('year','month'), sep='-') %>% 
  mutate(site = paste0(site, '_', pop, '_', qu)) %>% 
  count(site, year) %>% 
  ggplot() + 
  geom_point(aes(x = year, y = site)) +
  theme_bw() +
  labs(title    = 'Sampling inventory') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))


# Working data -----------------------------------------------------------------
df <- df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(qu = paste0(site, '_', pop, '_', qu))


# Explore qu -------------------------------------------------------------------
df_org %>% 
  group_by(site, pop) %>% 
  count(qu) %>% view()

'It appears that each site contains multiple populations. 
For example, site 2 includes only a single population (pop == 8). 
Within each population, there are multiple plots (qu), 
but these plots do not have unique identifiers across the full dataset.'


# Explore tag ------------------------------------------------------------------
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(plant = paste0(site, '_', pop, '_', qu, '_', plant)) %>% 
  arrange(site, pop, plant, year, month) %>% view()

''

# Data aggregation
df_org %>% 
  # there is dormancy but from what I saw so far it is only coded as s == 1
  #  but definitive we have the information in what month the fire was
  filter(s == 1 | !is.na(dia) | !is.na(burn)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(plant = paste0(site, '_', pop, '_', qu, '_', plant)) %>% 
  arrange(site, pop, plant, year, month) %>% view()

''

df_org %>% 
  filter(s > 1 & s < 9 & is.na(dia)) %>% 
  mutate(id = paste0(site, '_', pop, '_', qu, '_', plant)) %>% 
  pull(id) %>% 
  unique()

df_org %>% 
  mutate(id = paste0(site, '_', pop, '_', qu, '_', plant)) %>% 
  filter(id %in% (
    df_org %>% 
      filter(s > 1 & s < 9 & is.na(dia)) %>% 
      mutate(id = paste0(site, '_', pop, '_', qu, '_', plant)) %>% 
      pull(id) %>% 
      unique())) %>% view()

' s     == 8: dormancy 
  stage == 8: no clue'


# Prepare data frames for analysis ---------------------------------------------
# produce "transition" growth data
# IMPORTANT: ignore month for now: vast majority of sampling done in June!

# data from year t0
df_t0 <- df_org %>% 
          separate( date, c('year','month'), sep='-') %>% 
          mutate( year = as.numeric(year) ) %>% 
          # rename number of basal rosettes to show what time the data is from
          rename( nb_t0 = nb ) %>% 
          # retain only relevant variables
          dplyr::select(year, site, pop, qu, plant, nb_t0 ) 
df_t1 <- df_org %>% 
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

