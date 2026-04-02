# Data - Archbold -  - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.04.02


# Website    : 
# Publication: 


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(
  # negative binomial modeling
  MASS,
  # load tidyverse after MASS to not mask the select function
  tidyverse,
  # bbmle is for AICtab
  bbmle,
  # patchwork plot alingment
  patchwork,
  # binom.cofint for the survival plot
  binom,
  # tidy data
  janitor) # , skimr, ipmr, binom, lme4


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Chrysopsis highlandsensis')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()


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
# fig_gr
v_mod_set_su <- c()
# fig_su
v_mod_set_fl <- c()
# fig_fl
v_mod_set_fr <- c()


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
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    plant_id = identifier) %>% 
  mutate(
    plant_id = as.factor(plant_id),
    site     = as.factor(site))
# %>%
#   arrange(site, quad, quad_id, plant, plant_id, year, month)

# 
# df_meta <- data.frame(variable = colnames(df)) %>% 
#   mutate(definition = c(
#     'study site', 'quadrat number',	'macroplot number',	
#     'plant number (within quad)', 'direction within circular quad',	
#     'distance from quad center',	'was quad caged from 2012 onward',	
#     'vegetation type', 'year quadrat was initiated', 
#     'year-month of observation',	'fire severity Dec. 2014 or Jan. 2015',
#     'fire severity Aug. 2005',	'fire severity May/June 2009',
#     "fire severity B's ridge 2016", 'fire severity Feb. 2017',
#     'fire severity Oct. 2017', 'survival code for month', 'number of stems',
#     'number of branch tips', 'number of flowers (corolla showing)', 
#     'number of developing fruits', 'number of mature fruits', 'herbivory code',
#     'plant identification', 'quadrat identification', 'sample year', 
#     'sample month'))

write.csv(df,    row.names = F,
          file.path(dir_data, paste0('ab_', v_sp_abb, '_df.csv')))
