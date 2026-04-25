# Exploring data - Archbold - Menges 2016 - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.04.02


# Website    : 
# Publication: 

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
v_species <- c('Chrysopsis highlandsensis')
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
# Demographic data
df <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    plant_id = identifier) %>% 
  mutate(
    plant_id = as.factor(plant_id),
    site     = as.factor(site))

# Demographic meta
df_meta <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data_meta.csv')) %>% 
  janitor::clean_names()

# Site fire data
df_fire <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_fire.csv')) %>% 
  janitor::clean_names()

# Site fire meta
df_fire_meta <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_fire_meta.csv')) %>% 
  janitor::clean_names()
  
# Plot data
df_plot <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_plot.csv')) %>% 
  janitor::clean_names()

# Plot meta
df_plot_meta <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_plot_meta.csv')) %>% 
  janitor::clean_names()


# Summary ----------------------------------------------------------------------
df      %>% skim()
df_fire %>% skim()
df_plot %>% skim()


# Graphical exploration --------------------------------------------------------
# Sampling and treatment per site over the years
ggplot(df %>%
         left_join(
           df_fire %>% mutate(site = as.factor(site)),
           by = c("site" = "site", "year0" = "burn_yr"))) +
  geom_point(aes(x = year0, y = site, color = treatment)) +
  theme_minimal()

