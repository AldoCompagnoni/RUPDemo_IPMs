# Data exploration - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.22

# Study organism: Polygala lewtonii
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.318.1
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.318.1
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
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, MASS, GGally)


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
df_og <- read_csv(file.path(dir_data, 'polygala_lewtonii_data.csv')) %>% 
  janitor::clean_names()

df_meta <- data.frame(
  var = names(df_og),
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


# Sampling design --------------------------------------------------------------
df_og %>% 
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


# Explore quad ------------------------------------------------------------------
df_og %>% 
  group_by(unit) %>% 
  count(quad) #%>% view()
'Quad is indeed a unique identifier'


# Explore tag ------------------------------------------------------------------
df_og %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id)) %>% 
  arrange(unit, id, year, month) #%>% view()
''

# Exploring duplicates
df_og %>%
  separate(col = date, into = c('year', 'month'), sep = '-') %>%
  mutate(id = paste0(unit, '_', quad, '_', id)) %>%
  group_by(id, year, month) %>%
  # Count number of entries per date per id
  mutate(n_entries = n(),
         unique_angle_dist = n_distinct(paste(angle, dist))) %>%
  ungroup() %>%
  # Keep only those individuals who have >1 record in at least one date
  group_by(id) %>%
  filter(any(n_entries > 1)) %>%
  ungroup() %>%
  arrange(unit, id, year, month) #%>% view()
'The case of individual 1_609_109 indicates that 
some plant IDs appear multiple times in the dataset, 
occasionally with different angle and dist values. 
In some instances, an individual is recorded as dead (survival == 0/ == 9), 
and later, a new record with the same ID appears, 
potentially with different spatial coordinates and its own associated growth measurements. 
This pattern suggests that the same ID may be reused or reassigned to a new plant at a later time. 
Therefore, I suggest to use angle/ distance as part of the ID.

The publication does not mention anything about duplications (in the methods)'

df_og %>%
  separate(date, into = c("year", "month"), sep = "-") %>%
  mutate(id_full = paste0(unit, "_", quad, "_", id)) %>%
  group_by(id_full, year, month) %>%
  filter(n() > 1) %>%                    # Step 1: detect duplicates within same year-month
  ungroup() %>%
  group_by(id_full) %>%
  arrange(year, month, angle, dist, .by_group = TRUE) %>%
  slice_head(n = 2) %>%                 # Step 2: only keep first 2 rows *across all time*
  ungroup() %>%
  group_by(id_full) %>%
  summarise(has_5 = any(survival == 5),
            has_9 = any(survival == 9)) %>%
  filter(!(has_5 & has_9))
'there are 9 cases where the original plant is not dead and they have a new one with the same id'

df_og %>%
  separate(col = date, into = c('year', 'month'), sep = '-') %>%
  mutate(id = paste0(unit, '_', quad, '_', id)) %>%
  group_by(id, year, month) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  distinct(id) %>%
  inner_join(
    df_og %>%
      separate(col = date, into = c('year', 'month'), sep = '-') %>%
      mutate(id = paste0(unit, '_', quad, '_', id)),
    by = "id") %>%
  arrange(unit, id, year, month) #%>% view()

# Data aggregation
df_og %>% 
  # since there is no dormancy we can just filter for growth
  filter(!is.na(height)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) #%>% view()
'We found that there are actually only about 6500 observations in total'

# Exploring recruits
df_og %>% 
  filter(survival < 9,
         !c(survival == 1 & is.na(height))) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) #%>% view()
'There are recruits (assuming survival == 5 indicates a new recruit) 
that die without ever having a recorded height. 
These individuals should be excluded, as they can distort our analysis of recruit survival. 
Specifically, they contribute to the total count of recruits at the plot level but 
do not provide any data for survival within the recruitment size category.'

df_og %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>% 
  filter(id %in% (df_og %>% 
                       mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>%
                       filter(survival == 3) %>% 
                       pull(id) %>% 
                       unique())) #%>% view()
'Trying to figure out what survival == 3 means. 
*** Plant individual id == 3_539_108_130-6, was first encountered in March 2004, 
it first encounter was marked with a survivial == 3.
It does not have anything obvisouly to do with cohort, nor angle or distance.
The individual survives to the next observation (3 months later without size measurement).
Stage is marked as 3 (whatever that means) - but it is not therfore not a recruit?
Also it has a height of 8cm, whereas the other recruits are around 1-4 cm.
It only has 1 reporducing stem, 1 stem and 2 crowndiameter, all relatively low
*** Plant individaul is == 3_553_6_NA-NA, was first encounted in March 2011.
It dies in the next observation. It is at stage == 2. Height is relatively low at 2cm.
Crowndiameter of 1, stems of 2, and reproducitve stem of 0 is also relatively low.
*** Individual 4_502_127_NA-NA, March 2009, also dies with the next cencus. 
It is at stage 3, with 9cm height, crowndiameter of 5, 6 stems, and 3 reproductive stems.
*** many of them survive many years after'
df_og %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>%
  filter(survival == 3) %>% group_by(year, month) %>% count()
'We conculde that it is some measurement of how they were added to the sampling campain,
but sice all of them are not recruits we just convert them to survival == 1 and thats it
-> it is new individuals added (meta-data)'

# Exploring survival == 2
df_og %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>% 
  filter(id %in% (df_og %>% 
                    mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>%
                    filter(survival == 2) %>% 
                    pull(id) %>% 
                    unique())) #%>% view()
'it is 5 individuals in total
-> remove them from the data set!? -> yes!!'

# Exploring plant death
df_og %>% 
  filter(survival == 0 | !is.na(height)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) #%>% view()
''


# Generating data --------------------------------------------------------------
df_gen <- df_og %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>% #view()
  # Month of death
  group_by(id) %>%
  mutate(death = paste0(year, "_", month)[which(survival == 0)[1]]) %>%
  ungroup() %>% #view()
  # Code 2 survival: individual not found, 5 in total. REMOVE
  filter(!id %in% (df_og %>% 
                     mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>%
                     filter(survival == 2) %>% 
                     pull(id) %>% 
                     unique())) %>% #view()
# Fire: For now we will just drag the maximum value of each year with us
  dplyr::select(!c(angle, dist)) %>% #view()
  group_by(id, year) %>%
  mutate(postburn_plant = if_else(all(is.na(postburn_plant)), NA_real_,
                                  max(postburn_plant, na.rm = TRUE))) %>%
  ungroup() %>% #view()
  filter(survival < 9) %>% 
  # Recruits: 5 indicates that the plant was detected for the first time. 
  #  however, there are plants that were detected in a month without size measurement
  #  many of them are present in the next year too.. for now I will count them as recruits too
  #  actually I can give them a 2 - I like it, easy to change later too..
  # Furthermore, 1_609_101_100-19 is an example of how it is not considered a recruit in the original data set
  #  if the plant was there already in the fist year of the sampling campaign (i.e. 2001)
  #  We also don't know if it was or wasn't a recruit and thus it gets an NA
  group_by(id) %>%
  mutate(
    recruit = if_else(
      (lag(survival) == 5 & is.na(lag(stage))), 2, 0, missing = 0)) %>% 
  mutate(recruit = if_else(survival == 5, 1, recruit)) %>% 
  mutate(recruit = if_else(as.numeric(year) == 2001, NA_real_, recruit)) %>% 
  mutate(recruit = if_else(stage == 1, 1, recruit)) %>% 
  ungroup() %>% #view()
  # Survival: 
  #  Since there is no dormancy I can remove everything that does not have a size 
  filter(!is.na(stage)) %>% #view()
  group_by(id) %>%
  mutate(
    survival = if_else(row_number() == n(), 0, survival)) %>%
  ungroup() %>% #view()
  # Include the NEW ADDITIONS (survival == 3)
  mutate(survival = if_else(survival == 3, 1, survival),
         survival = if_else(survival == 5, 1, survival)) %>% 
  # Growth:
  group_by(id) %>%
  mutate(size_t1   = lead(height),
         c_dim_t1  = lead(crown_diameter),
         stems_t1  = lead(stems),
         volume_t0 = ((crown_diameter / 2) ^ 2) * pi * height,
         volume_t1 = lead(volume_t0)) %>%
  ungroup() #%>% view()


# Working data -----------------------------------------------------------------
df <- df_gen %>% 
  rename(survives = survival,
         size_t0  = height,
         site     = unit) %>% 
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),    
         logsize_t0_2 = logsize_t0^2,     
         logsize_t0_3 = logsize_t0^3,
         logvol_t0    = log(volume_t0),
         logvol_t1    = log(volume_t1),
         logvol_t0_2  = logvol_t0^2,
         logvol_t0_3  = logvol_t0^3,
         year         = as.numeric(year),
         stage        = as.factor( stage),
         recruits     = if_else(recruit > 0, 1, recruit),
         flower       = if_else(stage == 3, 1, 0)) %>%
  dplyr::select(site, quad, cohort, id, year, 
         stage, survives, size_t0, flower, flowering_stems, recruits, recruit, 
         size_t1, logsize_t1, logsize_t0, logsize_t0_2, logsize_t0_3,
         crown_diameter, stems, 
         volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3)


# GROWTH -----------------------------------------------------------------------
df_selected <- df %>%
  dplyr::select(height = size_t0, crown_diameter, stems, volume_t0) %>%
  filter(!if_any(everything(), is.na))  # remove rows with any NAs

# Create the pair plot
ggpairs(df_selected)


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  filter(
    !is.na(size_t1),      !is.na(size_t0), 
    !is.na(logsize_t0),   !is.na(logsize_t1), 
    !is.na(logsize_t0_2), !is.na(logsize_t0_3),
    !is.na(volume_t0),    !is.na(volume_t1), 
    !is.na(logvol_t0),    !is.na(logvol_t1), 
    !is.na(logvol_t0_2),  !is.na(logvol_t0_3)
    ) %>% 
  dplyr::select(
    id, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
    stage, recruits)

ggplot(
  data  = df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8)) + 
  facet_wrap('stage')

ggplot(df_gr) + geom_point(aes(logvol_t0, logsize_t1))


# Growth model -----------------------------------------------------------------
mod_gr_0 <- lm(logsize_t1 ~ 1, data = df_gr)
mod_gr_1 <- lm(logsize_t1 ~ logsize_t0, data = df_gr)
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2, data = df_gr)  
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, data = df_gr)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- AICctab(mods_gr, weights = T, sort = F)$dAIC
mods_gr_sorted <- order(mods_gr_dAIC)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr}

mod_gr_bestfit         <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef           <- coef(mod_gr_bestfit)

df_gr$pred <- predict(mod_gr_bestfit, type = 'response')
fig_gr_line <- ggplot(
  df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, mod_gr_ranef), 
                color = line_color_pred_fun(mod_gr_ranef), lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_pred <- ggplot(
  df_gr, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1), color = 'red', lwd = 2) + 
  theme_bw()

fig_gr <- fig_gr_line + fig_gr_pred + plot_layout() 
fig_gr


# Growth model volume -----------------------------------------
mod_gr_vol_0 <- lm(logvol_t1 ~ 1, data = df_gr)
mod_gr_vol_1 <- lm(logvol_t1 ~ logvol_t0, data = df_gr)
mod_gr_vol_2 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr)  
mod_gr_vol_3 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr)

mods_gr_vol      <- list(mod_gr_vol_0, mod_gr_vol_1, mod_gr_vol_2, mod_gr_vol_3)
mods_gr_vol_dAIC <- AICctab(mods_gr_vol, weights = T, sort = F)$dAIC

mods_gr_vol_sorted       <- order(mods_gr_vol_dAIC)
mod_gr_vol_index_bestfit <- mods_gr_vol_sorted[1]
mod_gr_vol_bestfit       <- mods_gr_vol[[mod_gr_vol_index_bestfit]]
mod_gr_vol_ranef         <- coef(mod_gr_vol_bestfit)

df_gr$pred_vol <- predict(mod_gr_vol_bestfit, type = "response")
ggplot(df_gr, aes(x = logvol_t0, y = logvol_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = pred_vol), color = "red", size = 1.2) +
  labs(x = "volumen t0 (log)", y = "volumen t1 (log)")


# Growth volume with stage -----------------------------------------------------
# Stage as a category
mod_gr_vol_01 <- lm(logvol_t1 ~ stage, data = df_gr)
mod_gr_vol_11 <- lm(logvol_t1 ~ logvol_t0 + stage, data = df_gr)
mod_gr_vol_21 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + stage, data = df_gr)  
mod_gr_vol_31 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + stage, data = df_gr)
mod_gr_vol_12 <- lm(logvol_t1 ~ logvol_t0 * stage, data = df_gr)
mod_gr_vol_22 <- lm(logvol_t1 ~ logvol_t0 * stage + logvol_t0_2:stage, data = df_gr)  
mod_gr_vol_32 <- lm(logvol_t1 ~ logvol_t0 * stage + logvol_t0_2:stage + logvol_t0_3:stage, data = df_gr)

mods_gr_vol_2      <- list(
  mod_gr_vol_0,  mod_gr_vol_1,  mod_gr_vol_2,  mod_gr_vol_3,
  mod_gr_vol_01, mod_gr_vol_11, mod_gr_vol_21, mod_gr_vol_31,
                 mod_gr_vol_12, mod_gr_vol_22, mod_gr_vol_32)
mods_gr_vol_2_dAIC <- AICctab(mods_gr_vol_2, weights = T, sort = F)$dAIC

mods_gr_vol_2_sorted       <- order(mods_gr_vol_2_dAIC)
mod_gr_vol_2_index_bestfit <- mods_gr_vol_2_sorted[1]
mod_gr_vol_2_bestfit       <- mods_gr_vol_2[[mod_gr_vol_2_index_bestfit]]
mod_gr_vol_2_ranef         <- coef(mod_gr_vol_2_bestfit)

df_gr$pred_vol_2 <- predict(mod_gr_vol_2_bestfit, type = "response")
ggplot(df_gr, aes(x = logvol_t0, y = logvol_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = pred_vol_2, color = stage), size = 1.2) +
  labs(x = "volumen t0 (log)", y = "volumen t1 (log)", color = "Stage") +
  theme_minimal()


# Growth volume by stage -------------------------------------------------------
df_gr_s1     <- df_gr %>% filter(stage == 1)
df_gr_s2     <- df_gr %>% filter(stage == 2)
df_gr_s3     <- df_gr %>% filter(stage == 3)

mod_gr_vol_s1_0 <- lm(logvol_t1 ~ 1, data = df_gr_s1)
mod_gr_vol_s1_1 <- lm(logvol_t1 ~ logvol_t0, data = df_gr_s1)
mod_gr_vol_s1_2 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr_s1)  
mod_gr_vol_s1_3 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr_s1)

mods_gr_vol_s1              <- list(mod_gr_vol_s1_0, mod_gr_vol_s1_1, mod_gr_vol_s1_2, mod_gr_vol_s1_3)
mods_gr_vol_s1_dAIC         <- AICctab(mods_gr_vol_s1, weights = T, sort = F)$dAIC
mods_gr_vol_s1_sorted       <- order(mods_gr_vol_s1_dAIC)
mod_gr_vol_s1_index_bestfit <- mods_gr_vol_s1_sorted[1]
mod_gr_vol_s1_bestfit       <- mods_gr_vol_s1[[mod_gr_vol_s1_index_bestfit]]

df_gr_s1$pred_vol_s1 <- predict(mod_gr_vol_s1_bestfit, type = "response")
ggplot(df_gr_s1, aes(x = logvol_t0, y = logvol_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = pred_vol_s1), color = "red", size = 1.2) +
  labs(x = "volumen t0 (log)", y = "volumen t1 (log)")

mod_gr_vol_s2_0 <- lm(logvol_t1 ~ 1, data = df_gr_s2)
mod_gr_vol_s2_1 <- lm(logvol_t1 ~ logvol_t0, data = df_gr_s2)
mod_gr_vol_s2_2 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr_s2)  
mod_gr_vol_s2_3 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr_s2)

mods_gr_vol_s2              <- list(mod_gr_vol_s2_0, mod_gr_vol_s2_1, mod_gr_vol_s2_2, mod_gr_vol_s2_3)
mods_gr_vol_s2_dAIC         <- AICctab(mods_gr_vol_s2, weights = T, sort = F)$dAIC
mods_gr_vol_s2_sorted       <- order(mods_gr_vol_s2_dAIC)
mod_gr_vol_s2_index_bestfit <- mods_gr_vol_s2_sorted[1]
mod_gr_vol_s2_bestfit       <- mods_gr_vol_s2[[mod_gr_vol_s2_index_bestfit]]

df_gr_s2$pred_vol_s2 <- predict(mod_gr_vol_s2_bestfit, type = "response")
ggplot(df_gr_s2, aes(x = logvol_t0, y = logvol_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = pred_vol_s2), color = "red", size = 1.2) +
  labs(x = "volumen t0 (log)", y = "volumen t1 (log)")

mod_gr_vol_s3_0 <- lm(logvol_t1 ~ 1, data = df_gr_s3)
mod_gr_vol_s3_1 <- lm(logvol_t1 ~ logvol_t0, data = df_gr_s3)
mod_gr_vol_s3_2 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr_s3)  
mod_gr_vol_s3_3 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr_s3)

mods_gr_vol_s3              <- list(mod_gr_vol_s3_0, mod_gr_vol_s3_1, mod_gr_vol_s3_2, mod_gr_vol_s3_3)
mods_gr_vol_s3_dAIC         <- AICctab(mods_gr_vol_s3, weights = T, sort = F)$dAIC
mods_gr_vol_s3_sorted       <- order(mods_gr_vol_s3_dAIC)
mod_gr_vol_s3_index_bestfit <- mods_gr_vol_s3_sorted[1]
mod_gr_vol_s3_bestfit       <- mods_gr_vol_s3[[mod_gr_vol_s3_index_bestfit]]

df_gr_s3$pred_vol_s3 <- predict(mod_gr_vol_s3_bestfit, type = "response")
ggplot(df_gr_s3, aes(x = logvol_t0, y = logvol_t1)) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1) +
  geom_line(aes(y = pred_vol_s3), color = "red", size = 1.2) +
  labs(x = "volumen t0 (log)", y = "volumen t1 (log)")


# Growth volume with recruits --------------------------------------------------
df_gr_r <- df_gr %>% filter(
  !is.na(recruits)) %>%
  mutate(recruits = as.factor(recruits))

mod_gr_vol_r_0  <- lm(logvol_t1 ~ 1, data = df_gr_r)
mod_gr_vol_r_1  <- lm(logvol_t1 ~ logvol_t0, data = df_gr_r)
mod_gr_vol_r_2  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr_r)  
mod_gr_vol_r_3  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr_r)
mod_gr_vol_r_01 <- lm(logvol_t1 ~ recruits, data = df_gr_r)
mod_gr_vol_r_11 <- lm(logvol_t1 ~ logvol_t0 + recruits, data = df_gr_r)
mod_gr_vol_r_21 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + recruits, data = df_gr_r)  
mod_gr_vol_r_31 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits, data = df_gr_r)
mod_gr_vol_r_12 <- lm(logvol_t1 ~ logvol_t0 * recruits, data = df_gr_r)
mod_gr_vol_r_22 <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits, data = df_gr_r)  
mod_gr_vol_r_32 <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, data = df_gr_r)

mods_gr_vol_r      <- list(
  mod_gr_vol_r_0,  mod_gr_vol_r_1,  mod_gr_vol_r_2,  mod_gr_vol_r_3,
  mod_gr_vol_r_01, mod_gr_vol_r_11, mod_gr_vol_r_21, mod_gr_vol_r_31,
  mod_gr_vol_r_12, mod_gr_vol_r_22, mod_gr_vol_r_32)
mods_gr_vol_r_dAIC <- AICctab(mods_gr_vol_r, weights = T, sort = F)$dAIC

mods_gr_vol_r_sorted       <- order(mods_gr_vol_r_dAIC)
mod_gr_vol_r_index_bestfit <- mods_gr_vol_r_sorted[1]
mod_gr_vol_r_bestfit       <- mods_gr_vol_r[[mod_gr_vol_r_index_bestfit]]
mod_gr_vol_r_ranef         <- coef(mod_gr_vol_r_bestfit)

df_gr_r_newdata <- df_gr_r %>%
  group_by(recruits) %>%
  summarise(range = list(seq(min(logvol_t0), max(logvol_t0), length.out = 100)), .groups = 'drop') %>%
  unnest(range) %>%
  mutate(
    logvol_t0 = range,
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)
  ) %>%
  dplyr::select(-range) %>%
  mutate(predicted = predict(mod_gr_vol_r_bestfit, newdata = .))

ggplot(df_gr_r, aes(x = logvol_t0, y = logvol_t1, color = recruits)) +
  geom_point(alpha = 0.4) +
  geom_line(data = df_gr_r_newdata, aes(y = predicted), size = 1) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  labs(
    title    = 'Growth prediction',
    subtitle = v_ggp_suffix,
    x = 'Volumen t0 (log)',
    y = 'Volumen t1 (log)',
    color = 'Recruit'
  ) +
  theme_minimal()


# Survival ---------------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives), !is.na(logsize_t0), !is.na(logsize_t0_2), !is.na(logsize_t0_3)) %>%
  filter(size_t0 != 0) %>%
  dplyr::select(id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                stage)

fig_su_overall <- ggplot(
  data = plot_binned_prop(df, 10, logsize_t0, survives)) +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logsize_t0, y = survives),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Survival',
       subtitle = v_ggp_suffix,
       x = expression('log(size)'[t0]),
       y = expression('Survival to time t1'))
fig_su_overall 


# Survival per stage -----------------------------------------------------------
fun_make_survival_plot <- function(df_stage, bin_data, stage_label) {
  ggplot(data = bin_data) +
    geom_jitter(data = df_stage, aes(x = logsize_t0, y = survives), 
                position = position_jitter(width = 0.1, height = 0.3),
                alpha = 0.1) +
    geom_point(aes(x = logsize_t0, y = survives),
               alpha = 1, pch = 16, color = 'red') +
    geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                  linewidth = 0.5, width = 0.5) +
    scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          title = element_text(size = 10),
          plot.subtitle = element_text(size = 8)) +
    labs(title = paste("Survival — Stage", stage_label),
         subtitle = v_ggp_suffix,
         x = expression('log(size)'[t0]),
         y = expression('Survival to time t1'))}

# Filter and bin each stage
df_su_s1     <- df_su %>% filter(stage == 1)
df_su_s1_bin <- plot_binned_prop(df_su_s1, 10, logsize_t0, survives)

df_su_s2     <- df_su %>% filter(stage == 2)
df_su_s2_bin <- plot_binned_prop(df_su_s2, 10, logsize_t0, survives)

df_su_s3     <- df_su %>% filter(stage == 3)
df_su_s3_bin <- plot_binned_prop(df_su_s3, 10, logsize_t0, survives)

# Create plots
fig_su_s1 <- fun_make_survival_plot(df_su_s1, df_su_s1_bin, 1)
fig_su_s2 <- fun_make_survival_plot(df_su_s2, df_su_s2_bin, 2)
fig_su_s3 <- fun_make_survival_plot(df_su_s3, df_su_s3_bin, 3)

# Combine with patchwork
fig_su_s <- fig_su_s1 + fig_su_s2 + fig_su_s3 + 
  plot_layout(nrow = 1)


# Survival model ---------------------------------------------------------------
# Logistic regression
mod_su_0 <- glm(survives ~ 1,
                data = df_su, family = 'binomial') 
# Logistic regression
mod_su_1 <- glm(survives ~ logsize_t0,
                data = df_su, family = 'binomial') 
# Quadratic logistic model
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2,
                data = df_su, family = 'binomial')  
# Cubic logistic model
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                data = df_su, family = 'binomial')  


# Compare models using AIC
mods_su      <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- AICctab(mods_su, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_su_sorted <- order(mods_su_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index       <- mod_su_index_bestfit - 1 
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index       <- v_mod_set_su
}

mod_su_bestfit   <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef         <- coef(mod_su_bestfit)

# Generate predictions for survival across a range of sizes
mod_su_x <- seq(
  min(df_su$logsize_t0, na.rm = T),
  max(df_su$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_su_pred <- predictor_fun(mod_su_x, mod_su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_su_x, survives = .)

# Survival plots
fig_su_line <- ggplot() +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives),
            color = line_color_pred_fun(mod_su_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_su_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su

ggsave(file.path(dir_result, 'mean_survival.png'), 
       plot = fig_su, width = 10, height = 5, dpi = 300)


# Survival by volume -----------------------------------------------------------
df_su_vol <- df %>% 
  filter(!is.na(survives), !is.na(logvol_t0), !is.na(logvol_t0_2), !is.na(logvol_t0_3),
         !is.na(recruits)) %>%
  filter(size_t0 != 0) %>%
  mutate(recruits = as.factor(recruits)) %>% 
  dplyr::select(id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
                stage, recruits)

fig_su_vol <- ggplot(
  data = plot_binned_prop(df_su_vol, 10, logvol_t0, survives)) +
  geom_jitter(data = df_su_vol, aes(x = logvol_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logvol_t0, y = survives),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logvol_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Survival',
       subtitle = v_ggp_suffix,
       x = expression('log(volume)'[t0]),
       y = expression('Survival to time t1'))
fig_su_vol


# Survival by volume model -----------------------------------------------------
mod_su_vol_0 <- glm(survives ~ 1, data = df_su_vol, family = 'binomial') 
mod_su_vol_1 <- glm(survives ~ logvol_t0, data = df_su_vol, family = 'binomial') 
mod_su_vol_2 <- glm(survives ~ logvol_t0 + logvol_t0_2, data = df_su_vol, family = 'binomial')  
mod_su_vol_3 <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_su_vol, family = 'binomial')  

mods_su_vol        <- list(mod_su_vol_0, mod_su_vol_1, mod_su_vol_2, mod_su_vol_3)
mods_su_vol_dAIC   <- AICctab(mods_su_vol, weights = T, sort = F)$dAIC
mods_su_vol_sorted <- order(mods_su_vol_dAIC)

mod_su_vol_index_bestfit <- mods_su_vol_sorted[1]
mod_su_vol_bestfit       <- mods_su_vol[[mod_su_vol_index_bestfit]]
mod_su_vol_ranef         <- coef(mod_su_vol_bestfit)

mod_su_vol_x <- seq(
  min(df_su_vol$logvol_t0, na.rm = T),
  max(df_su_vol$logvol_t0, na.rm = T), length.out = 100)
df_su_vol_pred <- predictor_fun(mod_su_vol_x, mod_su_vol_ranef) %>% 
  boot::inv.logit() %>% 
  data.frame(logvol_t0 = mod_su_vol_x, survives = .)

fig_su_vol_line <- ggplot() +
  geom_jitter(data = df_su_vol, aes(x = logvol_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_vol_pred, aes(x = logvol_t0, y = survives),
            color = line_color_pred_fun(mod_su_vol_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_su_vol_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_su_vol, 10, logvol_t0, survives), 
    aes(x = logvol_t0, y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df_su_vol, 10, logvol_t0, survives), 
    aes(x = logvol_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_su_vol_pred, aes(x = logvol_t0, y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_su_vol <- fig_su_vol_line + fig_su_vol_bin + plot_layout()
fig_su_vol

fig_su_line + fig_su_bin + fig_su_vol_line + fig_su_vol_bin + plot_layout()


# Survival volume with stage -----------------------------------------------------
# Stage as a category
mod_su_vol_01 <- glm(survives ~ stage, data = df_su_vol, family = 'binomial')
mod_su_vol_11 <- glm(survives ~ logvol_t0 + stage, data = df_su_vol, family = 'binomial')
mod_su_vol_21 <- glm(survives ~ logvol_t0 + logvol_t0_2 + stage, data = df_su_vol, family = 'binomial')  
mod_su_vol_31 <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + stage, data = df_su_vol, family = 'binomial')
mod_su_vol_12 <- glm(survives ~ logvol_t0 * stage, data = df_su_vol, family = 'binomial')
mod_su_vol_22 <- glm(survives ~ logvol_t0 * stage + logvol_t0_2:stage, data = df_su_vol, family = 'binomial')  
mod_su_vol_32 <- glm(survives ~ logvol_t0 * stage + logvol_t0_2:stage + logvol_t0_3:stage, data = df_su_vol, family = 'binomial')

mods_su_vol_2      <- list(
  mod_su_vol_0,  mod_su_vol_1,  mod_su_vol_2,  mod_su_vol_3,
  mod_su_vol_01, mod_su_vol_11, mod_su_vol_21, mod_su_vol_31,
  mod_su_vol_12, mod_su_vol_22, mod_su_vol_32)
mods_su_vol_2_dAIC <- AICctab(mods_su_vol_2, weights = T, sort = F)$dAIC

mods_su_vol_2_sorted       <- order(mods_su_vol_2_dAIC)
mod_su_vol_2_index_bestfit <- mods_su_vol_2_sorted[1]
mod_su_vol_2_bestfit       <- mods_su_vol_2[[mod_su_vol_2_index_bestfit]]
mod_su_vol_2_ranef         <- coef(mod_su_vol_2_bestfit)

mod_su_vol_2_x <- seq(
  min(df_su_vol$logvol_t0, na.rm = TRUE),
  max(df_su_vol$logvol_t0, na.rm = TRUE), length.out = 100)
df_su_vol_2_pred <- bind_rows(
  data.frame(logvol_t0 = mod_su_vol_2_x, stage = "1"),
  data.frame(logvol_t0 = mod_su_vol_2_x, stage = "2"),
  data.frame(logvol_t0 = mod_su_vol_2_x, stage = "3")
) %>%
  mutate(
    logvol_t0_2 = logvol_t0^2)

df_su_vol_2_pred$survives <- predict(mod_su_vol_2_bestfit, newdata = df_su_vol_2_pred, type = "response")

fig_su_vol_2_line <- ggplot() +
  geom_jitter(data = df_su_vol, aes(x = logvol_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_vol_2_pred, aes(x = logvol_t0, y = survives, color = stage), size = 1.2) +
  labs(title = "Survival prediction",
       subtitle = v_ggp_suffix,
       color = "Stage") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 8))

fig_su_vol_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_su_vol, 10, logvol_t0, survives), 
    aes(x = logvol_t0, y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df_su_vol, 10, logvol_t0, survives), 
    aes(x = logvol_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_su_vol_pred, aes(x = logvol_t0, y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_su_vol <- fig_su_vol_line + fig_su_vol_bin + plot_layout()
fig_su_vol

fig_su_line + fig_su_bin + fig_su_vol_line + fig_su_vol_bin + plot_layout()


# Survival volume with recruits ------------------------------------------------
mod_su_vol_r_0  <- glm(survives ~ 1, data = df_su_vol, family = 'binomial') 
mod_su_vol_r_1  <- glm(survives ~ logvol_t0, data = df_su_vol, family = 'binomial') 
mod_su_vol_r_2  <- glm(survives ~ logvol_t0 + logvol_t0_2, data = df_su_vol, family = 'binomial')  
mod_su_vol_r_3  <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_su_vol, family = 'binomial')  
mod_su_vol_r_01 <- glm(survives ~ recruits, data = df_su_vol, family = 'binomial')
mod_su_vol_r_11 <- glm(survives ~ logvol_t0 + recruits, data = df_su_vol, family = 'binomial')
mod_su_vol_r_21 <- glm(survives ~ logvol_t0 + logvol_t0_2 + recruits, data = df_su_vol, family = 'binomial')  
mod_su_vol_r_31 <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits, data = df_su_vol, family = 'binomial')
mod_su_vol_r_12 <- glm(survives ~ logvol_t0 * recruits, data = df_su_vol, family = 'binomial')
mod_su_vol_r_22 <- glm(survives ~ logvol_t0 * recruits + logvol_t0_2:recruits, data = df_su_vol, family = 'binomial')  
mod_su_vol_r_32 <- glm(survives ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, data = df_su_vol, family = 'binomial')

mods_su_vol_r      <- list(
  mod_su_vol_r_0,  mod_su_vol_r_1,  mod_su_vol_r_2,  mod_su_vol_r_3,
  mod_su_vol_r_01, mod_su_vol_r_11, mod_su_vol_r_21, mod_su_vol_r_31,
                   mod_su_vol_r_12, mod_su_vol_r_22, mod_su_vol_r_32)
mods_su_vol_r_dAIC <- AICctab(mods_su_vol_r, weights = T, sort = F)$dAIC

mods_su_vol_r_sorted       <- order(mods_su_vol_r_dAIC)
mod_su_vol_r_index_bestfit <- mods_su_vol_r_sorted[1]
mod_su_vol_r_bestfit       <- mods_su_vol_r[[mod_su_vol_r_index_bestfit]]
mod_su_vol_r_ranef         <- coef(mod_su_vol_r_bestfit)

df_su_vol_r_newdata <- df_su_vol %>%
  group_by(recruits) %>%
  summarise(x = list(seq(min(logvol_t0), max(logvol_t0), length.out = 100)), .groups = "drop") %>%
  unnest(x) %>%
  mutate(
    logvol_t0 = x,
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)) %>% 
  mutate(predicted = predict(mod_su_vol_r_bestfit, newdata = ., type = "response"))

fig_su_vol_r_line <- ggplot(df_su_vol, aes(x = logvol_t0, y = survives, color = factor(recruits))) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3) +
  geom_line(data = df_su_vol_r_newdata, aes(y = predicted), size = 1) +
  labs(
    title    = "Survival probability by volume and recruits",
    x        = "Volumen t0 (log)",
    y        = "Probability of survival",
    color    = "Recruits"
  ) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  theme_bw() +
  theme(legend.position = "none")

df_su_vol_r_bindata <- df_su_vol %>%
  group_by(recruits) %>%
  group_modify(~ plot_binned_prop(.x, 10, logvol_t0, survives)) %>%
  ungroup() %>%
  mutate(recruits = factor(recruits))

df_su_vol_r_pred <- df_su_vol %>%
  group_by(recruits) %>%
  reframe(logvol_t0 = seq(min(logvol_t0), max(logvol_t0), length.out = 100)) %>%
  mutate(
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)
  ) %>%
  mutate(
    survives = predict(mod_su_vol_r_bestfit, newdata = ., type = 'response')
  )

fig_su_vol_r_bin <- ggplot() +
  geom_point(data = df_su_vol_r_bindata, aes(x = logvol_t0, y = survives, color = recruits)) +
  geom_errorbar(
    data = df_su_vol_r_bindata,
    aes(x = logvol_t0, ymin = lwr, ymax = upr, color = recruits),
    width = 0.1
  ) +
  geom_line(
    data = df_su_vol_r_pred,
    aes(x = logvol_t0, y = survives, color = recruits),
    size = 1.2
  ) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  labs(
    x      = 'Volumen t0 (log)',
    y      = '',
    color  = 'Recruit'
  ) +
  theme_bw() +
  ylim(0, 1)

fig_su_vol_r <- fig_su_vol_r_line + fig_su_vol_r_bin + plot_layout()
fig_su_vol_r


# Flower data ------------------------------------------------------------------
# We exclude all recruits form the probablity to flower
df_fl <- df %>% 
  filter(!is.na(flower), !is.na(logvol_t0), !is.na(logvol_t0_2), !is.na(logvol_t0_3)) %>%
  filter(recruits == 0) %>% 
  dplyr::select(id, year, size_t0, flower, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
                stage)

fig_fl_overall <- ggplot(
  data = plot_binned_prop(df_fl, 10, logsize_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower), alpha = 0.1, 
              position = position_jitter(width = 0.1, height = 0.3)) +
  geom_point(aes(x = logsize_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text     = element_text(size = 8),
        title         = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Flowering probability in t0'))
fig_fl_overall


# Flower by volume -------------------------------------------------------------
ggplot(
  data = plot_binned_prop(df_fl, 10, logvol_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logvol_t0, y = flower), alpha = 0.1, 
              position = position_jitter(width = 0.1, height = 0.3)) +
  geom_point(aes(x = logvol_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logvol_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text     = element_text(size = 8),
        title         = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('log(volume)'[t0]),
       y        = expression('Flowering probability in t0'))


# Flower by volume model -------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ 1,
                data = df_fl, family = 'binomial') 
# Logistic regression
mod_fl_1 <- glm(flower ~ logvol_t0,
                data = df_fl, family = 'binomial') 
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logvol_t0 + logvol_t0_2,
                data = df_fl, family = 'binomial')  
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
                data = df_fl, family = 'binomial')  


# Compare models using AIC
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- AICctab(mods_fl, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fl_sorted <- order(mods_fl_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index       <- mod_fl_index_bestfit - 1 
} else {
  mod_fl_index_bestfit <- v_mod_set_fl +1
  v_mod_fl_index       <- v_mod_set_fl
}

mod_fl_bestfit <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef   <- coef(mod_fl_bestfit)

# Generate predictions for survival across a range of sizes
mod_fl_x <- seq(
  min(df_fl$logvol_t0, na.rm = T),
  max(df_fl$logvol_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_fl_pred <- predictor_fun(mod_fl_x, mod_fl_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logvol_t0 = mod_fl_x, flower = .)

# Survival plots
fig_fl_line <- ggplot() +
  geom_jitter(data = df_fl, aes(x = logvol_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred, aes(x = logvol_t0, y = flower),
            color = line_color_pred_fun(mod_fl_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Flowering prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_fl_bin <- ggplot() +
  geom_point(data = plot_binned_prop(df_fl, 10, logvol_t0, flower), 
             aes(x = logvol_t0, y = flower)) +
  geom_errorbar(
    data = plot_binned_prop(df_fl, 10, logvol_t0, flower), 
    aes(x = logvol_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_fl_pred, aes(x = logvol_t0, y = flower),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_fl <- fig_fl_line + fig_fl_bin + plot_layout()
fig_fl


# Fecundity volume -------------------------------------------------------------
# Conditional on flowering
df_fec <- df %>%
  filter(flower == 1)

# Since there are no 0s in the dataset we go for a truncated nb model
df_fec$flowering_stems %>% summary()
'I couldnt find a functioning truncated nb function'

mod_fe_nb <- glm.nb(flowering_stems ~ logvol_t0, data = df_fec)

df_fec_preddata <- tibble(
  logvol_t0 = seq(min(df_fec$logvol_t0, na.rm = TRUE),
                  max(df_fec$logvol_t0, na.rm = TRUE),
                  length.out = 100))
df_fec_preddata$predicted_stems <- predict(mod_fe_nb, newdata = df_fec_preddata, type = 'response')

# Plot
ggplot(df_fec, aes(x = logvol_t0, y = flowering_stems)) +
  geom_jitter(width = 0.1, height = 0.2, alpha = 0.4) +
  geom_line(data = df_fec_preddata, aes(y = predicted_stems), color = 'darkgreen', size = 1.2) +
  labs(
    title = 'Fecundity Model (Negative Binomial)',
    x = 'Volume t0 (log)',
    y = 'Number of Flowering Stems'
  ) +
  theme_minimal()


# Flowering stock to Recruit transition ----------------------------------------
stocks_by_plot <- df %>%
  filter(!is.na(flowering_stems)) %>%
  group_by(site, quad, year) %>%
  summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = "drop") %>%
  mutate(year_recruits = year + 1)

recruits_by_plot <- df %>%
  filter(recruits == 1) %>%
  group_by(site, quad, year) %>%
  summarise(recruits_present = 1, .groups = "drop")  # binary presence

stocks_to_recruits <- stocks_by_plot %>%
  left_join(recruits_by_plot, by = c("site", "quad", "year_recruits" = "year")) %>%
  mutate(recruits_present = ifelse(is.na(recruits_present), 0, recruits_present))

mod_stocks_to_recruits <- glm(recruits_present ~ total_stocks, data = stocks_to_recruits, family = binomial)
summary(mod_stocks_to_recruits)

ggplot(stocks_to_recruits, aes(x = total_stocks, y = recruits_present)) +
  geom_jitter(height = 0.05, alpha = 0.4) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "blue", se = TRUE) +
  labs(
    title = "Probability of Recruitment by Total Flowering Stems",
    x = "Total Flowering Stems (per site/plot/year)",
    y = "Recruitment Presence (Next Year)"
  ) +
  theme_minimal()


mod_log <- glm(recruits_present ~ log1p(total_stocks), data = stocks_to_recruits, family = binomial)
summary(mod_log)




# Aggregate flowering stems by year and shift forward
stock_to_recruit_by_year <- {
  
  # Total number of flowering stems per year
  stocks_by_year <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE)) %>%
    mutate(year = year + 1)  # Shift forward: flowering stems in year → recruits in year + 1
  
  # Count recruits by year
  recruits_by_year <- df %>%
    filter(recruit == 1) %>%
    group_by(year) %>%
    summarise(n_recruits = n())
  
  # Combine and compute per-stock recruitment
  recruits_by_year %>%
    left_join(stocks_by_year, by = "year") %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)
}

stock_to_recruit_by_year %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE)
  )

hist(stock_to_recruit_by_year$recruits_per_stock,
     main = "Recruitment Efficiency (Recruits per Flowering Stem)",
     xlab = "Recruits per Stock", col = "lightblue", border = "white")


repr_pc_mean_stock   <- mean(stock_to_recruit_by_year$recruits_per_stock, na.rm = TRUE)
repr_pc_median_stock <- median(stock_to_recruit_by_year$recruits_per_stock, na.rm = TRUE)


stock_to_recruit_by_year_site <- {
  # Total flowering stems by site and year
  stocks_by_year_site <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(site, year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = "drop") %>%
    mutate(year = year + 1)
  
  # Recruits by site and year
  recruits_by_year_site <- df %>%
    filter(recruit == 1) %>%
    group_by(site, year) %>%
    summarise(n_recruits = n(), .groups = "drop")
  
  # Join and compute ratio
  recruits_by_year_site %>%
    left_join(stocks_by_year_site, by = c("site", "year")) %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)
}

stock_to_recruit_by_year_plot <- {
  # Flowering stems per plot and year
  stocks_by_year_plot <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(site, quad, year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = "drop") %>%
    mutate(year = year + 1)
  
  # Recruits per plot and year
  recruits_by_year_plot <- df %>%
    filter(recruit == 1) %>%
    group_by(site, quad, year) %>%
    summarise(n_recruits = n(), .groups = "drop")
  
  # Join and compute
  recruits_by_year_plot %>%
    left_join(stocks_by_year_plot, by = c("site", "quad", "year")) %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)
}

stock_to_recruit_by_year_site <- stock_to_recruit_by_year_site %>%
  filter(!is.na(total_stocks) & total_stocks > 0)

stock_to_recruit_by_year_plot <- stock_to_recruit_by_year_plot %>%
  filter(!is.na(total_stocks) & total_stocks > 0)

# Summary for site-level
stock_to_recruit_by_year_site %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE)
  )

# Summary for plot-level
stock_to_recruit_by_year_plot %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE)
  )
