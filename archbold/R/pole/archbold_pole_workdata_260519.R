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
# Time periode: 2001-2025


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
load_packages(MASS, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, GGally, tidyverse)


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
df_oglong <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_oglong_260519.csv')))

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


# Build transition disturbance data ----
# Disturbance was recorded at the quadrat level. Missing assessments during
# known fire events are retained as NA rather than coded as undisturbed,
# because absence of data does not confirm absence of fire. This affects
# 35 annual records: 2 in 2019 and 33 in 2021. Years without a fire event
# are coded as 0.
df_dist_event <- df_oglong %>%
  group_by(quad) %>%
  summarise(
    dist_2001 = if_else(
      all(is.na(postburn_status1201)), NA_real_,
      as.numeric(any(postburn_status1201 > 0, na.rm = TRUE))),
    dist_2007 = if_else(
      all(is.na(postburn_pct1207)), NA_real_,
      as.numeric(any(postburn_pct1207 > 0, na.rm = TRUE))),
    dist_2009 = if_else(
      all(is.na(postburn_pct0609)), NA_real_,
      as.numeric(any(postburn_pct0609 > 0, na.rm = TRUE))),
    dist_2015 = if_else(
      all(is.na(postburn_pct0615)), NA_real_,
      as.numeric(any(postburn_pct0615 > 0, na.rm = TRUE))),
    dist_2019 = if_else(
      all(is.na(postburn_pct0319)), NA_real_,
      as.numeric(any(postburn_pct0319 > 0, na.rm = TRUE))),
    dist_2020 = if_else(
      all(is.na(postburn_pct0120)), NA_real_,
      as.numeric(any(postburn_pct0120 > 0, na.rm = TRUE))),
    dist_2022 = if_else(
      all(is.na(postburn_pct0322)), NA_real_,
      as.numeric(any(postburn_pct0322 > 0, na.rm = TRUE))),
    .groups = "drop")


# The postburn variables were recorded after the corresponding fire event.
# Fire occurring in March is assigned to the annual transition beginning in
# the preceding census year.
df_dist_transition <- bind_rows(
  df_dist_event %>%
    transmute(quad, year = 2001, dist_transition = dist_2001),
  df_dist_event %>%
    transmute(quad, year = 2007, dist_transition = dist_2007),
  df_dist_event %>%
    transmute(quad, year = 2009, dist_transition = dist_2009),
  df_dist_event %>%
    transmute(quad, year = 2015, dist_transition = dist_2015),
  df_dist_event %>%
    transmute(quad, year = 2018, dist_transition = dist_2019),
  df_dist_event %>%
    transmute(quad, year = 2019, dist_transition = dist_2020),
  df_dist_event %>%
    transmute(quad, year = 2021, dist_transition = dist_2022))

# Recruits ---------------------------------------------------------------------
df_recruit <- df_og %>%
  filter(qsurv == 5) %>%
  transmute(
    id,
    recruit_year = case_when(
      year == 2001 & month == 4 ~ 2001,
      month > 3 ~ year + 1,
      TRUE ~ year)) %>%
  distinct(id, recruit_year) %>%
  mutate(recruit_qsurv_t0 = 1)


# Death during annual transition ----------------------------------------------
df_death_interval <- df_og %>%
  filter(!is.na(stg)) %>%
  transmute(
    id, year,
    census_time = year * 12 + month,
    target_time = (year + 1) * 12 + 3) %>%
  left_join(
    df_og %>%
      filter(qsurv == 0) %>%
      transmute(id, death_time = year * 12 + month),
    by = "id", relationship = "many-to-many") %>%
  group_by(id, year) %>%
  summarise(
    death_in_interval = any(
      !is.na(death_time) &
        death_time > census_time &
        death_time <= target_time),
    .groups = "drop")


# Generating data --------------------------------------------------------------
df_gen <- df_og %>%
  arrange(id, year, month) %>%
  # Retain annual census records
  filter(!is.na(stg)) %>%
  # Recruitment
  left_join(
    df_recruit,
    by = c("id", "year" = "recruit_year")) %>%
  mutate(
    recruit = case_when(
      year == 2001 ~ NA_real_,
      recruit_qsurv_t0 == 1 ~ 1,
      TRUE ~ 0)) %>%
  # Death during the following annual transition
  left_join(df_death_interval, by = c("id", "year")) %>%
  mutate(
    death_in_interval = replace_na(death_in_interval, FALSE)) %>%
  group_by(id) %>%
  arrange(year, month, .by_group = TRUE) %>%
  mutate(
    next_year = lead(year),
    last_observed_year = max(year, na.rm = TRUE),
    annual_transition = coalesce(next_year == year + 1, FALSE),
    observed_later = last_observed_year > year,
    valid_growth_transition =
      annual_transition & !death_in_interval,
    survives = case_when(
      death_in_interval ~ 0,
      annual_transition ~ 1,
      observed_later ~ 1,
      TRUE ~ NA_real_),
    size_t1 = if_else(
      valid_growth_transition, lead(ht), NA_real_),
    c_dim_t1 = if_else(
      valid_growth_transition, lead(mcd), NA_real_),
    stems_t1 = if_else(
      valid_growth_transition, lead(st), NA_real_),
    volume_t0 = ((mcd / 2)^2) * pi * ht,
    volume_t1 = if_else(
      valid_growth_transition, lead(volume_t0), NA_real_)) %>%
  ungroup() %>%
  # Disturbance
  left_join(df_dist_transition, by = c("quad", "year")) %>%
  mutate(
    dist_transition = case_when(
      year %in% c(2001, 2007, 2009, 2015, 2018, 2019, 2021) ~
        dist_transition,
      TRUE ~ 0))


# Working data -----------------------------------------------------------------
df <- df_gen %>%
  rename(size_t0 = ht, site = patch, fl_nr = flst) %>%
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
         recruits     = if_else(recruit > 0, 1, recruit),
         flower       = if_else(fl_nr > 0, 1, fl_nr)) %>%
  dplyr::select(
    site, quad, id, year, stage, survives, size_t0, flower, fl_nr,
    recruits, recruit, size_t1, logsize_t1, logsize_t0,
    logsize_t0_2, logsize_t0_3, dist_transition,
    mcd, st, volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3)


# Save data --------------------------------------------------------------------
# write.csv(df_og, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original.csv')))
# write.csv(df_meta, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_meta.csv')))
write.csv(df, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata_260519.csv')))

