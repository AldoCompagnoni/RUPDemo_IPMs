# Data exploration - Archbold - Hypericum cumulicola

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.21

# Study organism: Hypericum cumulicola
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.181.2
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.181.2
# Citing publication: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2745.13206
# Time periode: 1994-2015
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
v_species <- c('Hypericum cumulicola')
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
df_org <- read_csv(file.path(dir_data, 'hypericum_cumulicola_data.csv')) %>% 
  janitor::clean_names()

df_meta <- data.frame(
  var = names(df_org),
  describtion = c(
    'year of encounter', 'patch id', 'gap id [plot]', 'plant id', 
    'fire year [previous fire]', 'year of collection', 
    'plant stage [t0]', 'maximun stem height [t0]', 
    'reproductive structures [t0]', 'number of reprductive stems [t0]', 
    'plant stage [t1]', 'maximun stem height [t1]', 
    'reproductive structures [t1]'))


# Tails from meta --------------------------------------------------------------
"censused during their peak of reproduction annually 
in July and August between 1994 and 2015"

"we tagged and evaluated ~100 plants. In years when < 100 plants were found, 
  new tagged plants were added (~20-40 individuals) to the existing sample"

"six possible stages: dead (0), alive (1), 
  yearling (5; new plants just reaching its first census), 
  new plant (3; plants of unknown age), missing (2), and previously dead (9)"

"total number of reproductive structures (flowers and fruits) 
for every plant in 1994-1999, 2001, 2002, 2010-2012, 
and counted reproductive structures for a sample of individuals in other years"

"we did not include any plant dormancy"


# Sampling design --------------------------------------------------------------
df_org %>% 
  group_by(site, time) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>%
  # Create a scatter plot of quadrat counts over the years
  ggplot() + 
  geom_point(aes(x = time, y = site)) +
  theme_bw() +
  labs(title    = 'Sampling inventory') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))


# Working data -----------------------------------------------------------------
df <- df_org %>% 
  # Sort the data frame
  arrange(site, tag, time) %>%
  # We don't need the year of first encounter
  select(!year) %>% 
  # Sampling year 'time' is what we want as year
  rename(year = time) %>% 
  # gap and tag as a unique identifier
  mutate(tag = paste0(site, '_', gap, '_', tag),
         gap = paste0(site, '_', gap))
  

# Explore gap ------------------------------------------------------------------
df_org %>% 
  group_by(site) %>% 
  count(gap) #%>% view()

'Within each site, multiple distinct gap values are present. 
These values generally increase consecutively starting from 1, 
resembling a form of enumeration—
suggesting that gap likely represents individual plots within each site.'


# Explore tag ------------------------------------------------------------------
df_org %>% 
  mutate(tag = paste0(site, '_', gap, '_', tag)) %>% 
  arrange(site, tag, time) #%>% view()

'The tag variable overlaps across different sites and gaps (i.e., plots), 
so we construct a unique identifier by combining site, gap, and tag. 
Additionally, we found that the second measurements for 
height, survival status, and reproductive structures
correspond to time t+1 (i.e., the subsequent time point).'

X_tag <- df_org %>% 
  mutate(tag = paste0(site, '_', gap, '_', tag)) %>% 
  arrange(site, tag, time) %>% 
  group_by(site, tag) %>% 
  count(time) %>% 
  filter(n > 1) #%>% view()

df_org %>% 
  # Sort the data frame
  arrange(site, tag, time) %>% 
  mutate(tag = paste0(site, '_', gap, '_', tag)) %>% 
  filter(tag %in% X_tag$tag) %>% view()

'It appears there was an issue with the original tag variable. 
After creating unique identifiers for each plant, we found that many tags still 
have multiple observations recorded within the same year (i.e. time), 
suggesting possible duplication or repeated measurements.


