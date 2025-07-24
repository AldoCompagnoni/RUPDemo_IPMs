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


# Explore quad ------------------------------------------------------------------
df_org %>% 
  group_by(unit) %>% 
  count(quad) %>% view()
'Quad is indeed a unique identifier'


# Explore tag ------------------------------------------------------------------
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id)) %>% 
  arrange(unit, id, year, month) #%>% view()
''

# Exploring duplicates
df_org %>%
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
Therefore, I suggest to use angle/ distance as part of the ID.'

# Data aggregation
df_org %>% 
  # since there is no dormancy we can just filter for growth
  filter(!is.na(height)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) #%>% view()
'We found that there are actually only about 6500 observations in total'

# Exploring recruits
df_org %>% 
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

df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>% 
  filter(id %in% (df_org %>% 
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
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>%
  filter(survival == 3) %>% group_by(year, month) %>% count()
'We conculde that it is some measurement of how they were added to the sampling campain,
but sice all of them are not recruits we just convert them to survival == 1 and thats it
-> it is new individuals added (meta-data)'

# Exploring survival == 2
df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>% 
  filter(id %in% (df_org %>% 
                    mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>%
                    filter(survival == 2) %>% 
                    pull(id) %>% 
                    unique())) #%>% view()
'it is 5 individuals in total
-> remove them from the data set!?'

# Exploring plant death
df_org %>% 
  filter(survival == 0 | !is.na(height)) %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) #%>% view()
''


# Generating data --------------------------------------------------------------
df_gen <- df_org %>% 
  separate(col = date, into = c('year', 'month'), sep = '-') %>% 
  mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>% 
  arrange(unit, id, year, month) %>%
  # Code 2 survival: individual not found, 5 in total. REMOVE
  filter(!id %in% (df_org %>% 
                     mutate(id = paste0(unit, '_', quad, '_', id, '_', angle, '-', dist)) %>%
                     filter(survival == 2) %>% 
                     pull(id) %>% 
                     unique())) %>% #view()
# Fire: For now we will just drag the maximum value of each year with us
  select(!c(angle, dist)) %>% #view()
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
  mutate(recruit = if_else(as.numeric(year) == 2001, NA_real_, recruit)) %>% #view()
  filter(survival < 5) %>% 
  ungroup() %>%
  # Survival: 
  #  Since there is no dormancy I can remove everything that does not have a size 
  filter(!is.na(stage)) %>% #view()
  group_by(id) %>%
  mutate(
    survival = if_else(row_number() == n(), 0, survival)) %>%
  ungroup() %>% #view()
  # Include the NEW ADDITIONS (survival == 3)
  mutate(survival = if_else(survival == 3, 1, survival)) %>% 
  # Growth:
  group_by(id) %>%
  mutate(size_t1    = lead(height),
         c_dim_t1   = lead(crown_diameter),
         stems_t1   = lead(stems),
         volumen_t0 = ((crown_diameter / 2) ^ 2) * pi * height,
         volumen_t1 = lead(volumen_t0)) %>%
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
         logvol_t0    = log(volumen_t0),
         logvol_t1    = log(volumen_t1)) %>% 
  select(site, quad, cohort, id, year, 
         stage, survives, size_t0, flowering_stems, recruit, 
         size_t1, logsize_t1, logsize_t0, logsize_t0_2, logsize_t0_3,
         crown_diameter, stems, 
         volumen_t0, volumen_t1, logvol_t0, logvol_t1)


# GROWTH -----------------------------------------------------------------------
df_selected <- df %>%
  select(height = size_t0, crown_diameter, stems, volumen) %>%
  filter(!if_any(everything(), is.na))  # remove rows with any NAs

# Create the pair plot
ggpairs(df_selected)


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  dplyr::select(id, year, size_t0, size_t1,
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volumen_t0, volumen_t1, logvol_t0, logvol_t1,
                stage)

ggplot(
  data  = df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8)) + 
  facet_wrap('stage')

ggplot(df_gr) + geom_point(aes(logsize_t0, size_t1))
ggplot(df_gr) + geom_point(aes(logvol_t0, logsize_t1))
ggplot(df_gr) + geom_point(aes(logvol_t0, volumen_t1))


# Growth model -----------------------------------------------------------------
# Intercept model 
mod_gr_0 <- lm(logsize_t1 ~ 1,
               data = df_gr)
# Linear model
mod_gr_1 <- lm(logsize_t1 ~ logsize_t0, 
               data = df_gr)
# Quadratic model
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2, 
               data = df_gr)  
# Cubic model
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
               data = df_gr)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- AICtab(mods_gr, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_gr_sorted <- order(mods_gr_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr
}

mod_gr_bestfit         <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef           <- coef(mod_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_gr$pred <- predict(mod_gr_bestfit, type = 'response')

# Growth plot
fig_gr_line <- ggplot(
  df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
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


# Growth model negative binomial -----------------------------------------------
mod_gr_nb_0 <- glm.nb(size_t1 ~ 1, data = df_gr)
mod_gr_nb_1 <- glm.nb(size_t1 ~ logsize_t0, data = df_gr)
mod_gr_nb_2 <- glm.nb(size_t1 ~ logsize_t0 + logsize_t0_2, data = df_gr)
mod_gr_nb_3 <- glm.nb(size_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, data = df_gr)

mods_gr_nb        <- list(  mod_gr_nb_0, mod_gr_nb_1, mod_gr_nb_2, mod_gr_nb_3)
mods_gr_nb_dAIC   <- AICtab(mods_gr_nb, weights = T, sort = F)$dAIC
mods_gr_nb_sorted <- order(mods_gr_nb_dAIC)

mod_gr_nb_index_bestfit <- mods_gr_nb_sorted[1]
mod_gr_nb_bestfit       <- mods_gr_nb[[mod_gr_nb_index_bestfit]]
mod_gr_nb_ranef         <- coef(mod_gr_nb_bestfit)

df_gr$nb_predicted <- predict(mod_gr_nb_2, type = "response")
ggplot(df_gr, aes(x = logsize_t0, y = size_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = nb_predicted), color = "red", size = 1.2) +
  labs(x = "log(size_t0)",
       y = "size_t1")

# Survival ---------------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  dplyr::select(id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

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



