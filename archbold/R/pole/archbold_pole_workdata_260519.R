# Data - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.05.19

# Study organism: Polygala lewtonii
# Link: 
# Meta data link: 
# Citing publication: https://doi.org/10.1071/BT11271
# Time periode: 2001-2017


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
load_packages(MASS, tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, GGally)


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
v_mod_set_su <- c()
v_mod_set_gr <- c()
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
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df_og <- read_csv(file.path(dir_data, 'ab_pole_df_original_260519.csv')) %>% 
  janitor::clean_names() %>% 
  rename(id = unique_id)
df_site <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_sitehist_260519.csv')))

df_meta <- data.frame(
  var = names(df_og),
  describtion = c(
    'Unique identifier for each individual plant; Numeric',
    'Name of circular quadrat with radius of 0.25m; Numeric',
    'Name of an aggregate of several quads; Numeric',
    'Year data were collected; Numeric',
    'Month data were collected for March, June, September, and December. April was also used in the first sampling; Numeric',
    'Combination of year and month; MMYY',
    'Quarterly survival code. 0 or 1 refer to survival from the previous quarter. If a plant was dropped from monitoring, survival is NA and these plants should NOT be considered dead. 0 = dead; 1 = alive; 2 = missing; 3 = new adult; 5 = seedling',
    'Stage class from March census. If plant was a seedling during the previous year it is considered a seedling here. 1 = putative seedling; 2 = vegetative; 3 = reproductive adult',
    'Height from March census in cm; Numeric',
    'Maximum crown diameter from March census in cm; Numeric',
    'Total number of stems; Numeric',
    'Number of flowering stems; Numeric'))


# Tails from meta --------------------------------------------------------------
"Plants were censused during their peak of reproduction annually in July and August"

"In 220 permanent 25 cm radius circular plots, plants were marked with plastic toothpick"

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


# Sampling inventory -----------------------------------------------------------
fig_sampling <- df_og %>% 
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
fig_sampling

# ggsave(file.path(dir_result, '0_sampling_inventory.png'), 
#        plot = fig_sampling, width = 10, height = 5, dpi = 300)


# Exploring duplicates
df_og %>%
  group_by(id, year, month) %>%
  # Count number of entries per date per id
  mutate(n_entries = n()) %>%
  ungroup() %>%
  # Keep only those individuals who have >1 record in at least one date
  group_by(id) %>%
  filter(any(n_entries > 1)) %>%
  ungroup() %>%
  arrange(id, year, month) #%>% view()


# Data aggregation
df_og %>% 
  # since there is no dormancy we can just filter for growth
  filter(!is.na(ht)) %>% 
  arrange(id, year, month) #%>% view()
'We found that there are actually only about 10441 observations in total'


# Exploring recruits
df_og %>% 
  filter(qsurv == 5) %>%
  arrange(id, year, month) #%>% view()
'There are recruits that seem to die without ever having a recorded height. 
These individuals should be excluded, as they can distort our analysis of recruit survival. 
Specifically, they contribute to the total count of recruits at the plot level but 
do not provide any data for survival within the recruitment size category.'


# Exploring plant death
df_og %>% 
  filter(qsurv == 0) %>% 
  arrange(id, year, month) #%>% view()
'we need to adapt our survival to t1'


# Generating data --------------------------------------------------------------
df_gen <- df_og %>%
  arrange(id, year, month) %>% #view()
  # Month of death
  group_by(id) %>%
  mutate(death = paste0(year, "_", month)[which(qsurv == 0)[1]]) %>%
  # ungroup() %>% 
  # group_by(id, year) %>%
  # mutate(postburn_plant = if_else(all(is.na(postburn_plant)), NA_real_,
  #                                 max(postburn_plant, na.rm = TRUE))) %>%
  # ungroup() %>% #view()
  # filter(survival < 9) %>% 
  # # Recruits: 5 indicates that the plant was detected for the first time. 
  # #  however, there are plants that were detected in a month without size measurement
  # #  many of them are present in the next year too.. for now I will count them as recruits too
  # #  actually I can give them a 2 - I like it, easy to change later too..
  # # Furthermore, 1_609_101_100-19 is an example of how it is not considered a recruit in the original data set
  # #  if the plant was there already in the fist year of the sampling campaign (i.e. 2001)
  # #  We also don't know if it was or wasn't a recruit and thus it gets an NA
  group_by(id) %>%
  mutate(
    recruit = if_else(
      (lag(qsurv) == 5 & is.na(lag(stg))), 2, 0, missing = 0)) %>% 
  mutate(recruit = if_else(qsurv == 5, 1, recruit)) %>% 
  mutate(recruit = if_else(as.numeric(year) == 2001, NA_real_, recruit)) %>% 
  mutate(recruit = if_else(stg == 1, 1, recruit)) %>% 
  ungroup() %>% #view()
  # Survival: 
  #  Since there is no dormancy I can remove everything that does not have a size 
  filter(!is.na(stg)) %>% #view()
  group_by(id) %>%
  mutate(
    qsurv = if_else(row_number() == n(), 0, qsurv)) %>%
  ungroup() %>% #view()
  # Include the NEW ADDITIONS (qsurv == 3)
  mutate(qsurv = if_else(qsurv == 3, 1, qsurv),
         qsurv = if_else(qsurv == 5, 1, qsurv)) %>% 
  # Growth:
  group_by(id) %>%
  mutate(size_t1   = lead(ht),
         c_dim_t1  = lead(mcd),
         stems_t1  = lead(st),
         volume_t0 = ((mcd / 2) ^ 2) * pi * ht,
         volume_t1 = lead(volume_t0)) %>%
  ungroup() #%>% view()


# Working data -----------------------------------------------------------------
df <- df_gen %>% 
  rename(survives = qsurv,
         size_t0  = ht,
         site     = quad,
         flower   = flst) %>% 
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),    
         logsize_t0_2 = logsize_t0^2,     
         logsize_t0_3 = logsize_t0^3,
         logvol_t0    = log(volume_t0),
         logvol_t1    = log(volume_t1),
         logvol_t0_2  = logvol_t0^2,
         logvol_t0_3  = logvol_t0^3,
         year         = as.numeric(year),
         stage        = as.factor( stg),
         recruits     = if_else(recruit > 0, 1, recruit)) %>%
  dplyr::select(site, id, year, 
                stage, survives, size_t0, flower, recruits, recruit, 
                size_t1, logsize_t1, logsize_t0, logsize_t0_2, logsize_t0_3,
                mcd, st,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3)


# Save data --------------------------------------------------------------------
# write.csv(df_og, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original.csv')))
# write.csv(df_meta, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_meta.csv')))
write.csv(df, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata_260519.csv')))

