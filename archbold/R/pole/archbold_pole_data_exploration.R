# Data exploration - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.22

# Study organism: Polygala lewtonii
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.318.1
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.318.1
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
v_species <- c('Polygala lewtonii')
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
df_org <- read_csv(file.path(dir_data, 'polygala_lewtonii_data.csv')) %>% 
  janitor::clean_names()

df_meta <- data.frame(
  var = names(df_org),
  describtion = c(
    'census year and month', 'management unit', 
    '25 cm circular quadrat unique ID number', 
    'cohort month and year when plant first recruited', 
    'plant id', 'degree angle of plant within plot', 
    'distance (cm) of plant from plot center', 'survival', 
    'life history stage', 'stem height (cm)', 
    'maximum crown diameter (cm)', 'number of stems', 
    'number of flowering stems', 'status of the plant post-fire'))


# Tails from meta --------------------------------------------------------------
" In 220 permanent 25 cm radius circular plots, plants were marked with plastic toothpick"

"Data collection occurred quarterly with survival and recruitment recorded in 
March, June, September and December and more detailed measures of size and fecundity taken in March"

"This species is amphicarphic with both showy-open flowers produced in spring (chasmogamy) 
and closed, self-pollinating flowers both above- and belowground produced in late summer (cleistogamy).
Therefore, our fecundity data are based only on chasmogamy."

"This dataset also captures post-fire mortality and recruitment events following four prescribed burns."

"Plants that had died were given a 'dead' code (0) 
indicating first time dead but toothpicks and marking flags remained in the field. 
At the next census, if the plant was still dead, 
it was given a 'previously dead' code (9) and the toothpick and marking flag were removed from the field."

"Plant dormancy does not occur in this species and very few individuals marked as dead were alive in the subsequent census."

"In March, we took additional size measures"

"four prescribed burns during the 2001-2017 study, although not all plots were affected"

"status of the plant as unburned, singed (some green), 
scorched (no green but brown leaves remained), or consumed (little to no plant material remained)."


# Sampling design --------------------------------------------------------------
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>%
  mutate(quad = paste0(unit, '_', quad)) %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>%
  # Create a scatter plot of quadrat counts over the years
  ggplot() + 
  geom_point(aes(x = year, y = quad)) +
  theme_bw() +
  labs(title    = 'Sampling inventory') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))


# Working data -----------------------------------------------------------------
df <- df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-')


# Explore quad ------------------------------------------------------------------
df_org %>% 
  group_by(unit) %>% 
  count(quad) %>% view()

'Quad is indeed a unique identifier'


# Explore tag ------------------------------------------------------------------
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id)) %>% 
  arrange(unit, id, year, month) %>% view()

''

# Data aggregation
df_org %>% 
  # since there is no dormancy we can just filter for growth
  filter(!is.na(height)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id)) %>% 
  arrange(unit, id, year, month) %>% view()

'We found that there are actually only about 6500 observations in total'


