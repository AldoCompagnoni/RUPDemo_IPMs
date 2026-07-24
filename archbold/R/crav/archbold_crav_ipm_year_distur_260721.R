# by year IPM with distrubance - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.03


# Website    : 
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
  # GLMM: glmer function
  lme4, 
  janitor,
  glmmTMB,
  scales) # , skimr, ipmr, binom


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head             <- c('archbold')
# Define species
v_species          <- c('Crotalaria avonensis')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()
# Years that we want to remove from the analysis
v_years_re         <- c()


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
v_mod_set_fl_n <- c()
v_mod_set_fr <- c()
v_mod_set_fr_n <- c()
v_mod_set_fr2re <- c()



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
# function to plot your survival data 'binned' per year (instead of 'jittered')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df_og <- read_delim(
  file.path(dir_data, 'crotalaria_avonensis_data_v260515.csv'),
  delim = ';', escape_double = FALSE, trim_ws = TRUE)

df_meta <- data.frame(variable = colnames(df_og)) %>% 
  mutate(definition = c(
    'Unique identifier for each plant',
    'site # with the Carter Creek population',
    "macro plot #. Typically 5 quads per macroplot, and 3 mp's per site",	
    'Circular quadrat number, 0.5m diameter',
    'Year data were collected',
    'Survival code	0 = dead or dormant, 1 = alive, 2 = missing, 3 = new adult, 5 = new seedling, 6 = dead/run over, 8= dormant',
    'Whether plant was alive or not	0=dead, 1=alive',
    'Stage class	0=dead, 1=seedling, 2=vegetative and not flowering, 3=flowering',
    'Maximum number of branches observed at one of 3 time points (Feb, Apr, Jun)',
    'Maximum number of flowers observed at one of 3 time points (Feb, Apr, Jun)',
    'Maximum number of fruits (developing or ripe) observed at one of 3 time points (Feb, Apr, Jun)'))

df_og_quad <- read_delim(
  file.path(dir_data, 'crotalaria_avonensis_data_quad_v260515.csv'),
  delim = ';', escape_double = FALSE, trim_ws = TRUE)

df_disturbance <- df_og_quad %>%
  # keep only ID columns + burn columns
  select(site, mp, quad,
         burn221, burn1222, burn0123, burn05, burn2014_2015, burn2016,
         burn0217, burnoct2017, burnjun2018) %>%
  pivot_longer(
    cols = starts_with('burn'), names_to = 'burn_event', values_to = 'value') %>%
  mutate(
    disturbance = case_when(
      is.na(value) ~ 0,
      value == FALSE ~ 0,
      value == TRUE ~ 1,
      suppressWarnings(as.numeric(value)) > 0 ~ 1,
      TRUE ~ 0)) %>%
  filter(disturbance == 1) %>%
  # assign actual years
  mutate(
    year = case_when(
      burn_event == 'burn221' ~ 2021,
      burn_event == 'burn1222' ~ 2022,
      burn_event == 'burn0123' ~ 2023,
      burn_event == 'burn05' ~ 2005,
      burn_event == 'burn2014_2015' ~ 2015,
      burn_event == 'burn2016' ~ 2016,
      burn_event == 'burn0217' ~ 2017,
      burn_event == 'burnoct2017' ~ 2017,
      burn_event == 'burnjun2018' ~ 2018),
    disturbance = 1) %>%
  mutate(site = as.factor(site),
         mp   = as.factor(mp),
         quad = as.factor(quad)) %>%
  # final columns
  select(site, mp, quad, year, disturbance) %>%
  distinct() %>%
  arrange(site, mp, quad, year)


# Mean data frame --------------------------------------------------------------
df <- df_og %>% 
  janitor::clean_names() %>%  
  rename(plant_id = id,
         size_t0  = maxbr,
         flower   = maxfl,
         fruit    = maxfr) %>% 
  mutate(
    plant_id = as.factor(plant_id),
    quad_id  = as.factor(paste(site, mp, quad, sep = '_')),
    site     = as.factor(site),
    quad     = as.factor(quad),
    mp       = as.factor(mp)) %>%
  arrange(site, quad, quad_id, plant_id, year) %>%
  
  # Survival
  group_by(plant_id) %>% 
  mutate(survives = lead(alive)) %>% 
  ungroup %>% 
  
  # Growth based on survival; propagate size_t0 if survives, otherwise set to NA
  group_by(plant_id) %>%
  mutate(
    size_t1 = case_when(
      survives == 1 ~ lead(size_t0), 
      TRUE          ~ NA_real_)) %>%
  ungroup() %>% 
  
  # Recruits
  mutate(
    latest_alive_date      = max(year[(s > 0 & s != 2 & s < 6)], na.rm = TRUE),
    earliest_recorded_date = min(year[s > 0 & s != 2 & s != 6], na.rm = TRUE)) %>%
  mutate(
    recruit = case_when(
      stage == 1 ~ 1,
      year == earliest_recorded_date & earliest_recorded_date > 1999 + 4 ~ 1,
      TRUE ~ NA_real_)) %>% 
  
  # Dormancy
  group_by(plant_id) %>% 
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0)              ~ 1,
      size_t0  >  0 & !is.na(survives)            ~ 0,
      survives == 0 & year <= max(df_og$year) - 4 ~ 0, 
      TRUE                                        ~ NA_real_ ),
    # Generate a new column 'dormancy_count' that counts consecutive 1s
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 & 
        lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1 & lag(dormancy, 5) == 1 ~ 6,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1                         ~ 5,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1                                                 ~ 4,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1           ~ 3,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 0           ~ 2,
      dormancy == 1                                                           ~ 1,
      TRUE ~ dormancy)) %>%
  ungroup() %>% 
  # adjust the survival for potential dormancy
  mutate(survives = if_else(survives == 0 & year > max(df_og$year) - 4, NA, survives)) %>%
  
  # Log-transformed sizes
  mutate(
    logsize_t0   = log(size_t0),     
    logsize_t1   = log(size_t1),    
    logsize_t0_2 = logsize_t0^2,     
    logsize_t0_3 = logsize_t0^3) %>% 
  
  # Disturbance form the quad data
  left_join(df_disturbance, by = c('site', 'mp', 'quad', 'year')) %>%
  mutate(disturbance = ifelse(is.na(disturbance), 0, disturbance)) %>% 
  
  select(site, mp, quad, quad_id, plant_id, year, 
         s, survives, size_t0, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         flower, fruit, recruit, dormancy, disturbance)


# Recruit validity: previous-year plot sampling -------------------------------
df_plot_year_sampled <- df %>%
  mutate(
    year = as.integer(as.character(year)),
    quad_id = factor(quad_id)) %>%
  group_by(site, mp, quad, quad_id, year) %>%
  summarise(
    n_records_plot_year = n(),
    n_plants_plot_year = n_distinct(plant_id),
    .groups = "drop")

df_recruit_plot_check <- df %>%
  filter(recruit == 1) %>%
  mutate(
    year = as.integer(as.character(year)),
    year_prev = year - 1,
    quad_id = factor(quad_id)) %>%
  left_join(
    df_plot_year_sampled %>%
      mutate(sampled_prev_year = TRUE) %>%
      select(
        quad_id,
        year_prev = year,
        sampled_prev_year,
        n_records_prev_year = n_records_plot_year,
        n_plants_prev_year = n_plants_plot_year),
    by = c("quad_id", "year_prev")) %>%
  mutate(
    sampled_prev_year = replace_na(sampled_prev_year, FALSE),
    recruit_valid_plot_previous_year = sampled_prev_year)

df_recruit_valid_ids <- df_recruit_plot_check %>%
  filter(recruit_valid_plot_previous_year) %>%
  select(site, mp, quad, quad_id, plant_id, year) %>%
  mutate(recruit_plot_valid = 1)

df <- df %>%
  mutate(year = as.integer(as.character(year))) %>%
  left_join(
    df_recruit_valid_ids,
    by = c("site", "mp", "quad", "quad_id", "plant_id", "year")) %>%
  mutate(
    recruit_plot_valid = case_when(
      recruit == 1 & recruit_plot_valid == 1 ~ 1,
      TRUE ~ NA_real_))




# Survival data ----------------------------------------------------------------
df_su <- df %>%
  filter(
    size_t0 != 0,
    !is.na(survives),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, survives, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    disturbance)


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(
    size_t0 != 0,
    survives == 1,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, disturbance)


# Flower data ------------------------------------------------------------------
df_fl <- df %>%
  filter(
    !is.na(flower),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance),
    flower = if_else(flower > 0, 1, 0)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, flower, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    disturbance)


# Flower-number data -----------------------------------------------------------
# Flower number is conditional on flowering because flower > 0.
# Expected flowers per plant in the IPM:
# Pr(flowering) * E(number of flowers | flowering).

df_fl_n <- df %>%
  filter(
    flower > 0,
    is.finite(flower),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, flower,
    logsize_t0, logsize_t0_2, logsize_t0_3,
    disturbance)


# Fruiting-probability data ----------------------------------------------------
# Fruiting probability is conditional on flowering because flower > 0.

df_fr <- df %>%
  filter(
    flower > 0,
    is.finite(fruit),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance),
    fruiting = if_else(fruit > 0, 1, 0)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, flower, fruit, fruiting,
    logsize_t0, logsize_t0_2, logsize_t0_3, disturbance)


# Fruit-number data ------------------------------------------------------------
# Fruit number is conditional on fruiting because fruit > 0.

df_fr_n <- df %>%
  filter(
    fruit > 0,
    is.finite(fruit),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    plant_id, site, year, size_t0, flower, fruit,
    logsize_t0, logsize_t0_2, logsize_t0_3, disturbance)


# Fruit-to-recruit data --------------------------------------------------------
# Total fruits in year t predict recruits in year t + 1.
# Only quadrats sampled in consecutive years are retained.

df_fr2re_quad <- df %>%
  mutate(year = as.integer(as.character(year))) %>%
  group_by(site, quad_id, year) %>%
  summarise(
    fruits_t0 = sum(fruit, na.rm = TRUE),
    recruits_t0 = sum(recruit_plot_valid == 1, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(site, quad_id, year) %>%
  group_by(site, quad_id) %>%
  mutate(
    recruits_t1 = lead(recruits_t0),
    year_t1 = lead(year),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1)

df_fr2re <- df_fr2re_quad %>%
  group_by(site, year, year_t1) %>%
  summarise(
    fruits_t0 = sum(fruits_t0, na.rm = TRUE),
    recruits_t1 = sum(recruits_t1, na.rm = TRUE),
    n_quads = n(),
    .groups = "drop") %>%
  mutate(
    site = factor(site),
    year_t1 = factor(year_t1),
    logfruits_t0 = log1p(fruits_t0),
    logrecruits_t1 = log1p(recruits_t1))


# Removing years with too few data --------------------------------------------
v_years_og <- sort(unique(df$year))

df <- df %>% filter(!is.na(year), !(year %in% v_years_re))
df_su <- df_su %>% filter(!is.na(year), !(year %in% v_years_re)) %>%
  droplevels()
df_gr <- df_gr %>% filter(!is.na(year), !(year %in% v_years_re)) %>%
  droplevels()
df_fl <- df_fl %>%
  filter(!is.na(year), !(year %in% v_years_re)) %>% droplevels()
df_fl_n <- df_fl_n %>% filter(!is.na(year), !(year %in% v_years_re)) %>%
  droplevels()
df_fr <- df_fr %>% filter(!is.na(year), !(year %in% v_years_re)) %>%
  droplevels()
df_fr_n <- df_fr_n %>% filter(!is.na(year), !(year %in% v_years_re)) %>%
  droplevels()
df_fr2re <- df_fr2re %>%
  filter(
    !is.na(year),
    !is.na(year_t1),
    !(year %in% v_years_re),
    !(as.integer(as.character(year_t1)) %in% v_years_re)) %>%
  droplevels()

v_years <- sort(unique(df$year))


# Survival model ---------------------------------------------------------------
ctrl_glmer <- glmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

ctrl_lmer <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5))

# GLMM; binomial
mod_su_0 <- glmer(
  survives ~ disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_1 <- glmer(
  survives ~ logsize_t0 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mods_su <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- bbmle::AICctab(mods_su, weights = TRUE, sort = FALSE)$dAIC
mods_su_sorted <- order(mods_su_dAIC)

if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index <- mod_su_index_bestfit - 1
} else {
  mod_su_index_bestfit <- v_mod_set_su + 1
  v_mod_su_index <- v_mod_set_su
}

mod_su_best <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef <- coef(mod_su_best)$year

summary(mod_su_best)


# Survival plots by year -------------------------------------------------------
make_su_year_plot <- function(year_i) {
  df_i <- df_su %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_su$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      survives = predict(
        mod_su_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, survives, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, survives, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, survives, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = 'Survival probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

su_yrs <- lapply(levels(df_su$year), make_su_year_plot)

fig_su_years <- wrap_plots(su_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Survival - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_su_years


# Growth model -----------------------------------------------------------------
mod_gr_0 <- lmer(
  logsize_t1 ~ disturbance + (logsize_t0 | year),
  data = df_gr, control = ctrl_lmer)

mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + disturbance + (logsize_t0 | year),
  data = df_gr, control = ctrl_lmer)

mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance + (logsize_t0 | year),
  data = df_gr, control = ctrl_lmer)

mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (logsize_t0 | year),
  data = df_gr, control = ctrl_lmer)

mods_gr <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- bbmle::AICctab(mods_gr, weights = TRUE, sort = FALSE)$dAIC
mods_gr_sorted <- order(mods_gr_dAIC)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index <- mod_gr_index_bestfit - 1
} else {
  mod_gr_index_bestfit <- v_mod_set_gr + 1
  v_mod_gr_index <- v_mod_set_gr
}

mod_gr_best <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef <- coef(mod_gr_best)$year

summary(mod_gr_best)


# Growth plots by year ---------------------------------------------------------
make_gr_year_plot <- function(year_i) {
  df_i <- df_gr %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_gr$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      logsize_t1 = predict(
        mod_gr_best, newdata = pred_i, re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      logsize_t1 = mean(logsize_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, logsize_t1, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, logsize_t1, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, logsize_t1, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = expression('log(size)'[t1])) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

gr_yrs <- lapply(levels(df_gr$year), make_gr_year_plot)

fig_gr_years <- wrap_plots(gr_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Growth - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_gr_years


# Growth variance --------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_best)
mod_gr_y <- resid(mod_gr_best)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Flower model -----------------------------------------------------------------
mod_fl_0 <- glmer(
  flower ~ disturbance + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_1 <- glmer(
  flower ~ logsize_t0 + disturbance + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mods_fl <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- bbmle::AICctab(mods_fl, weights = TRUE, sort = FALSE)$dAIC
mods_fl_sorted <- order(mods_fl_dAIC)

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index <- mod_fl_index_bestfit - 1
} else {
  mod_fl_index_bestfit <- v_mod_set_fl + 1
  v_mod_fl_index <- v_mod_set_fl
}

mod_fl_best <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef <- coef(mod_fl_best)$year

summary(mod_fl_best)


# Flowering plots by year ------------------------------------------------------
make_fl_year_plot <- function(year_i) {
  df_i <- df_fl %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fl$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      flower = predict(
        mod_fl_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      flower = mean(flower, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, flower, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, flower, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, flower, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = 'Flowering probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fl_yrs <- lapply(levels(df_fl$year), make_fl_year_plot)

fig_fl_years <- wrap_plots(fl_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Flowering probability - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fl_years


# Flower-number model ----------------------------------------------------------
mod_fl_n_0 <- glmer.nb(
  flower ~ disturbance + (1 | year),
  data = df_fl_n, control = ctrl_glmer)

mod_fl_n_1 <- glmer.nb(
  flower ~ logsize_t0 + disturbance + (1 | year),
  data = df_fl_n, control = ctrl_glmer)

mod_fl_n_2 <- glmer.nb(
  flower ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fl_n, control = ctrl_glmer)

mod_fl_n_3 <- glmer.nb(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance +
    (1 | year),
  data = df_fl_n, control = ctrl_glmer)

mods_fl_n <- list(
  mod_fl_n_0, mod_fl_n_1, mod_fl_n_2, mod_fl_n_3)

mods_fl_n_dAIC <- bbmle::AICctab(
  mods_fl_n, weights = TRUE, sort = FALSE)$dAIC

mods_fl_n_sorted <- order(mods_fl_n_dAIC)

if (length(v_mod_set_fl_n) == 0) {
  mod_fl_n_index_bestfit <- mods_fl_n_sorted[1]
  v_mod_fl_n_index <- mod_fl_n_index_bestfit - 1
} else {
  mod_fl_n_index_bestfit <- v_mod_set_fl_n + 1
  v_mod_fl_n_index <- v_mod_set_fl_n
}

mod_fl_n_best <- mods_fl_n[[mod_fl_n_index_bestfit]]
mod_fl_n_ranef <- coef(mod_fl_n_best)$year

summary(mod_fl_n_best)


# Flower-number plots by year --------------------------------------------------
make_fl_n_year_plot <- function(year_i) {
  df_i <- df_fl_n %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fl_n$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      flower = predict(
        mod_fl_n_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      flower = mean(flower, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, flower, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, flower, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, flower, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = 'Number of flowers') +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fl_n_yrs <- lapply(levels(df_fl_n$year), make_fl_n_year_plot)

fig_fl_n_years <- wrap_plots(fl_n_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Flower number given flowering - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fl_n_years


# Fruiting-probability model ---------------------------------------------------
mod_fr_0 <- glmer(
  fruiting ~ disturbance + (1 | year),
  data = df_fr, family = binomial, control = ctrl_glmer)

mod_fr_1 <- glmer(
  fruiting ~ logsize_t0 + disturbance + (1 | year),
  data = df_fr, family = binomial, control = ctrl_glmer)

mod_fr_2 <- glmer(
  fruiting ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fr, family = binomial, control = ctrl_glmer)

mod_fr_3 <- glmer(
  fruiting ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance +
    (1 | year),
  data = df_fr, family = binomial, control = ctrl_glmer)

mods_fr <- list(mod_fr_0, mod_fr_1, mod_fr_2, mod_fr_3)
mods_fr_dAIC <- bbmle::AICctab(
  mods_fr, weights = TRUE, sort = FALSE)$dAIC
mods_fr_sorted <- order(mods_fr_dAIC)

if (length(v_mod_set_fr) == 0) {
  mod_fr_index_bestfit <- mods_fr_sorted[1]
  v_mod_fr_index <- mod_fr_index_bestfit - 1
} else {
  mod_fr_index_bestfit <- v_mod_set_fr + 1
  v_mod_fr_index <- v_mod_set_fr
}

mod_fr_best <- mods_fr[[mod_fr_index_bestfit]]
mod_fr_ranef <- coef(mod_fr_best)$year

summary(mod_fr_best)


# Fruiting-probability plots by year ------------------------------------------
make_fr_year_plot <- function(year_i) {
  df_i <- df_fr %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fr$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      fruiting = predict(
        mod_fr_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      fruiting = mean(fruiting, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, fruiting, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, fruiting, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, fruiting, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = 'Fruiting probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fr_yrs <- lapply(levels(df_fr$year), make_fr_year_plot)

fig_fr_years <- wrap_plots(fr_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Fruiting probability given flowering - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fr_years


# Fruit-number model -----------------------------------------------------------
mod_fr_n_0 <- glmer.nb(
  fruit ~ disturbance + (1 | year),
  data = df_fr_n, control = ctrl_glmer)

mod_fr_n_1 <- glmer.nb(
  fruit ~ logsize_t0 + disturbance + (1 | year),
  data = df_fr_n, control = ctrl_glmer)

mod_fr_n_2 <- glmer.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fr_n, control = ctrl_glmer)

mod_fr_n_3 <- glmer.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance +
    (1 | year),
  data = df_fr_n, control = ctrl_glmer)

mods_fr_n <- list(mod_fr_n_0, mod_fr_n_1, mod_fr_n_2, mod_fr_n_3)
mods_fr_n_dAIC <- bbmle::AICctab(
  mods_fr_n, weights = TRUE, sort = FALSE)$dAIC
mods_fr_n_sorted <- order(mods_fr_n_dAIC)

if (length(v_mod_set_fr_n) == 0) {
  mod_fr_n_index_bestfit <- mods_fr_n_sorted[1]
  v_mod_fr_n_index <- mod_fr_n_index_bestfit - 1
} else {
  mod_fr_n_index_bestfit <- v_mod_set_fr_n + 1
  v_mod_fr_n_index <- v_mod_set_fr_n
}

mod_fr_n_best <- mods_fr_n[[mod_fr_n_index_bestfit]]
mod_fr_n_ranef <- coef(mod_fr_n_best)$year

summary(mod_fr_n_best)


# Fruit-number plots by year ---------------------------------------------------
make_fr_n_year_plot <- function(year_i) {
  df_i <- df_fr_n %>%
    filter(as.character(year) == year_i)
  
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fr_n$year)),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred_i <- pred_i %>%
    mutate(
      fruit = predict(
        mod_fr_n_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(
      bin = cut(logsize_t0, breaks = 8),
      disturbance_plot = factor(
        disturbance, levels = c(0, 1),
        labels = c('No fire', 'Fire'))) %>%
    group_by(bin, disturbance_plot) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      fruit = mean(fruit, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logsize_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i %>%
        mutate(
          disturbance_plot = factor(
            disturbance, levels = c(0, 1),
            labels = c('No fire', 'Fire'))),
      aes(logsize_t0, fruit, color = disturbance_plot),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logsize_t0, fruit, color = disturbance_plot),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, fruit, color = disturbance_plot),
      linewidth = 0.7) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(size)'[t0]),
      y = 'Number of fruits') +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fr_n_yrs <- lapply(levels(df_fr_n$year), make_fr_n_year_plot)

fig_fr_n_years <- wrap_plots(fr_n_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Fruit number given fruiting - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fr_n_years


# Fruit-to-recruit models ------------------------------------------------------
# Zero fruits predict zero recruits.
# Recruitment varies among years, but not among sites.

mod_fr2re_0 <- glmmTMB(
  logrecruits_t1 ~ 0 + logfruits_t0,
  data = df_fr2re, family = gaussian, REML = FALSE)

mod_fr2re_1 <- glmmTMB(
  logrecruits_t1 ~ 0 + logfruits_t0 +
    (0 + logfruits_t0 | year_t1),
  data = df_fr2re, family = gaussian, REML = FALSE)

mods_fr2re <- list(mod_fr2re_0, mod_fr2re_1)

mods_fr2re_dAIC <- bbmle::AICctab(
  mods_fr2re, weights = TRUE, sort = FALSE)$dAIC
mods_fr2re_sorted <- order(mods_fr2re_dAIC)

if (length(v_mod_set_fr2re) == 0) {
  mod_fr2re_index_bestfit <- mods_fr2re_sorted[1]
  v_mod_fr2re_index <- mod_fr2re_index_bestfit - 1
} else {
  mod_fr2re_index_bestfit <- v_mod_set_fr2re + 1
  v_mod_fr2re_index <- v_mod_set_fr2re
}

mod_fr2re_best <- mods_fr2re[[mod_fr2re_index_bestfit]]
mod_fr2re_fixef <- fixef(mod_fr2re_best)$cond
mod_fr2re_ranef <- ranef(mod_fr2re_best)$cond

summary(mod_fr2re_best)


# Fruit-to-recruit plots by year ----------------------------------------------
make_fr2re_year_plot <- function(year_i) {
  df_i <- df_fr2re %>%
    filter(as.character(year_t1) == year_i)
  
  x <- seq(
    min(df_i$logfruits_t0, na.rm = TRUE),
    max(df_i$logfruits_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- tibble(
    logfruits_t0 = x,
    year_t1 = factor(year_i, levels = levels(df_fr2re$year_t1)))
  
  pred_i <- pred_i %>%
    mutate(
      logrecruits_t1 = predict(
        mod_fr2re_best, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE),
      recruits_t1 = pmax(0, expm1(logrecruits_t1)))
  
  pts_i <- df_i %>%
    mutate(bin = cut(logfruits_t0, breaks = 8)) %>%
    group_by(bin) %>%
    summarise(
      logfruits_t0 = mean(logfruits_t0, na.rm = TRUE),
      recruits_t1 = mean(recruits_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(is.finite(logfruits_t0), n > 0)
  
  ggplot() +
    geom_jitter(
      data = df_i,
      aes(logfruits_t0, recruits_t1),
      alpha = 0.2, size = 0.7, width = 0.2, height = 0.2) +
    geom_point(
      data = pts_i,
      aes(logfruits_t0, recruits_t1),
      shape = 15, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(logfruits_t0, recruits_t1),
      linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(number of fruits + 1)'[t]),
      y = expression('Number of recruits'[t+1])) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

fr2re_yrs <- lapply(
  levels(df_fr2re$year_t1), make_fr2re_year_plot)

fig_fr2re_years <- wrap_plots(fr2re_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Fruit-to-recruit relationship - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fr2re_years


# Recruit size distribution ----------------------------------------------------
df_recr <- df %>%
  filter(
    recruit_plot_valid == 1,
    size_t0 > 0,
    is.finite(logsize_t0))

rc_sz <- df_recr %>%
  summarise(
    mean = mean(logsize_t0, na.rm = TRUE),
    sd = sd(logsize_t0, na.rm = TRUE))


# Parameter extraction ---------------------------------------------------------
# pars_mean contains fixed-effect coefficients.
# pars_year contains complete conditional coefficients for each year.

get_model_fixef <- function(model) {
  if (inherits(model, 'glmmTMB')) {
    return(glmmTMB::fixef(model)$cond)
  }
  
  if (inherits(model, 'merMod')) {
    return(lme4::fixef(model))
  }
  
  stats::coef(model)
}


get_group_coef <- function(model, group_var) {
  if (inherits(model, 'glmmTMB')) {
    coef_list <- stats::coef(model)$cond
  } else if (inherits(model, 'merMod')) {
    coef_list <- stats::coef(model)
  } else {
    return(NULL)
  }
  
  if (is.null(coef_list) || !(group_var %in% names(coef_list))) {
    return(NULL)
  }
  
  as.data.frame(coef_list[[group_var]])
}


size_term_map <- c(
  b0 = '(Intercept)',
  b1 = 'logsize_t0',
  b2 = 'logsize_t0_2',
  b3 = 'logsize_t0_3',
  bd = 'disturbance')


extract_terms <- function(coefs, prefix, term_map = size_term_map) {
  values <- vapply(term_map, function(term_i) {
    if (term_i %in% names(coefs)) {
      unname(coefs[[term_i]])
    } else {
      0
    }
  }, numeric(1))
  
  names(values) <- paste0(prefix, names(term_map))
  values
}


extract_fixed_rate <- function(model, prefix) {
  extract_terms(
    coefs = get_model_fixef(model),
    prefix = prefix)
}


extract_year_rate <- function(
    model, prefix, years, group_var = 'year') {
  
  fixed_coef <- get_model_fixef(model)
  group_coef <- get_group_coef(model, group_var)
  
  purrr::map_dfr(years, function(year_i) {
    level_i <- as.character(year_i)
    coef_i <- fixed_coef
    
    if (!is.null(group_coef) &&
        level_i %in% rownames(group_coef)) {
      coef_i <- unlist(
        group_coef[level_i, , drop = FALSE],
        use.names = FALSE)
      
      names(coef_i) <- colnames(group_coef)
    }
    
    tibble::as_tibble_row(c(
      year = year_i,
      extract_terms(coef_i, prefix)))
  })
}


extract_fr2re_year <- function(model, years) {
  fixed_coef <- get_model_fixef(model)
  group_coef <- get_group_coef(model, 'year_t1')
  
  purrr::map_dfr(years, function(year_i) {
    year_t1_i <- year_i + 1
    level_i <- as.character(year_t1_i)
    coef_i <- fixed_coef
    
    if (!is.null(group_coef) &&
        level_i %in% rownames(group_coef)) {
      coef_i <- unlist(
        group_coef[level_i, , drop = FALSE],
        use.names = FALSE)
      
      names(coef_i) <- colnames(group_coef)
    }
    
    fr2re_b1 <- if ('logfruits_t0' %in% names(coef_i)) {
      unname(coef_i[['logfruits_t0']])
    } else {
      0
    }
    
    tibble(
      year = year_i,
      year_t1 = year_t1_i,
      fr2re_b1 = fr2re_b1)
  })
}


# Mean parameters --------------------------------------------------------------
grow_var_coef <- stats::coef(mod_gr_var)

mesh_limits <- range(
  c(
    df_gr$logsize_t0,
    df_gr$logsize_t1,
    df_recr$logsize_t0),
  na.rm = TRUE,
  finite = TRUE)

pars_mean <- as.list(c(
  extract_fixed_rate(mod_su_best, 'surv_'),
  extract_fixed_rate(mod_gr_best, 'grow_'),
  extract_fixed_rate(mod_fl_best, 'fl_'),
  extract_fixed_rate(mod_fl_n_best, 'fln_'),
  extract_fixed_rate(mod_fr_best, 'fr_'),
  extract_fixed_rate(mod_fr_n_best, 'frn_'),
  fr2re_b1 = unname(
    get_model_fixef(mod_fr2re_best)[['logfruits_t0']]),
  recr_mean = rc_sz$mean,
  recr_sd = rc_sz$sd,
  grow_var_a = unname(grow_var_coef[['a']]),
  grow_var_b = unname(grow_var_coef[['b']]),
  fruits_ref = mean(df_fr2re$fruits_t0, na.rm = TRUE),
  fr2re_sigma = stats::sigma(mod_fr2re_best),
  L = mesh_limits[1],
  U = mesh_limits[2],
  mat_siz = 200))


pars_mean_table <- tibble(
  parameter = names(pars_mean),
  value = as.numeric(unlist(pars_mean)))

pars_mean_table


# Years with all required vital-rate models -----------------------------------
get_year_values <- function(x) {
  sort(unique(as.integer(as.character(x))))
}

ipm_years <- Reduce(
  intersect,
  list(
    get_year_values(df_su$year),
    get_year_values(df_gr$year),
    get_year_values(df_fl$year),
    get_year_values(df_fl_n$year),
    get_year_values(df_fr$year),
    get_year_values(df_fr_n$year),
    get_year_values(df_fr2re$year_t1) - 1))

if (length(ipm_years) == 0) {
  stop('No years have complete data for all IPM processes')
}

ipm_years


# Year-specific parameter extraction ------------------------------------------
pars_su_year <- extract_year_rate(
  model = mod_su_best, prefix = 'surv_', years = ipm_years)

pars_gr_year <- extract_year_rate(
  model = mod_gr_best, prefix = 'grow_', years = ipm_years)

pars_fl_year <- extract_year_rate(
  model = mod_fl_best, prefix = 'fl_', years = ipm_years)

pars_fl_n_year <- extract_year_rate(
  model = mod_fl_n_best, prefix = 'fln_', years = ipm_years)

pars_fr_year <- extract_year_rate(
  model = mod_fr_best, prefix = 'fr_', years = ipm_years)

pars_fr_n_year <- extract_year_rate(
  model = mod_fr_n_best, prefix = 'frn_', years = ipm_years)

pars_fr2re_year <- extract_fr2re_year(
  model = mod_fr2re_best, years = ipm_years)

pars_year <- list(
  pars_su_year,
  pars_gr_year,
  pars_fl_year,
  pars_fl_n_year,
  pars_fr_year,
  pars_fr_n_year,
  pars_fr2re_year) %>%
  purrr::reduce(left_join, by = 'year') %>%
  arrange(year)

pars_year


# Create mean or year-specific parameter list ---------------------------------
make_ipm_pars <- function(year_i = NULL) {
  pars <- pars_mean
  
  if (is.null(year_i)) {
    return(pars)
  }
  
  year_i <- as.integer(year_i)
  
  pars_i <- pars_year %>%
    filter(.data$year == .env$year_i)
  
  if (nrow(pars_i) != 1) {
    stop(
      'Expected one parameter row for year ', year_i,
      ', but found ', nrow(pars_i))
  }
  
  pars_i <- pars_i %>%
    select(-year, -any_of('year_t1')) %>%
    as.list()
  
  for (par_name in names(pars_i)) {
    pars[[par_name]] <- pars_i[[par_name]]
  }
  
  pars
}


# IPM functions ----------------------------------------------------------------
get_par <- function(pars, par_name, default = 0) {
  if (!is.null(pars[[par_name]])) {
    pars[[par_name]]
  } else {
    default
  }
}


vital_lp <- function(x, pars, prefix, disturbance = 0) {
  get_par(pars, paste0(prefix, 'b0')) +
    get_par(pars, paste0(prefix, 'b1')) * x +
    get_par(pars, paste0(prefix, 'b2')) * x^2 +
    get_par(pars, paste0(prefix, 'b3')) * x^3 +
    get_par(pars, paste0(prefix, 'bd')) * disturbance
}


# Survival ---------------------------------------------------------------------
sx <- function(x, pars, disturbance = 0) {
  plogis(vital_lp(
    x = x,
    pars = pars,
    prefix = 'surv_',
    disturbance = disturbance))
}


# Growth -----------------------------------------------------------------------
grow_mu <- function(x, pars, disturbance = 0) {
  vital_lp(
    x = x,
    pars = pars,
    prefix = 'grow_',
    disturbance = disturbance)
}


grow_sd <- function(mean_value, pars) {
  variance <- get_par(pars, 'grow_var_a') *
    exp(get_par(pars, 'grow_var_b') * mean_value)
  
  sqrt(pmax(variance, .Machine$double.eps))
}


gxy <- function(x, y, pars, disturbance = 0) {
  mean_value <- grow_mu(
    x = x,
    pars = pars,
    disturbance = disturbance)
  
  sd_value <- grow_sd(
    mean_value = mean_value,
    pars = pars)
  
  dnorm(y, mean = mean_value, sd = sd_value)
}


pxy <- function(x, y, pars, disturbance = 0) {
  sx(x, pars, disturbance) *
    gxy(x, y, pars, disturbance)
}


# Flowering probability --------------------------------------------------------
fl_x <- function(x, pars, disturbance = 0) {
  plogis(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fl_',
    disturbance = disturbance))
}


# Flower number conditional on flowering --------------------------------------
fln_x <- function(x, pars, disturbance = 0) {
  exp(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fln_',
    disturbance = disturbance))
}


flowers_x <- function(x, pars, disturbance = 0) {
  fl_x(x, pars, disturbance) *
    fln_x(x, pars, disturbance)
}


# Fruiting probability conditional on flowering -------------------------------
fr_x <- function(x, pars, disturbance = 0) {
  plogis(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fr_',
    disturbance = disturbance))
}


# Fruit number conditional on fruiting ----------------------------------------
frn_x <- function(x, pars, disturbance = 0) {
  exp(vital_lp(
    x = x,
    pars = pars,
    prefix = 'frn_',
    disturbance = disturbance))
}


# Expected fruits per individual ----------------------------------------------
# Flower number is not multiplied here because the fruit models describe
# plant-level fruiting and total fruit number.

fruits_x <- function(x, pars, disturbance = 0) {
  fl_x(x, pars, disturbance) *
    fr_x(x, pars, disturbance) *
    frn_x(x, pars, disturbance)
}


# Recruitment at the reference fruit production -------------------------------
re_total_ref <- function(pars) {
  eta <- get_par(pars, 'fr2re_b1') *
    log1p(get_par(pars, 'fruits_ref'))
  
  recruits <- expm1(
    eta + 0.5 * get_par(pars, 'fr2re_sigma')^2)
  
  pmax(recruits, 0)
}


re_per_fruit <- function(pars) {
  re_total_ref(pars) /
    pmax(
      get_par(pars, 'fruits_ref'),
      .Machine$double.eps)
}


# Expected recruits contributed by an individual ------------------------------
rx <- function(x, pars, disturbance = 0) {
  fruits_x(
    x = x,
    pars = pars,
    disturbance = disturbance) *
    re_per_fruit(pars)
}


# Recruit size distribution ----------------------------------------------------
recr_y <- function(y, pars) {
  dnorm(
    y,
    mean = get_par(pars, 'recr_mean'),
    sd = get_par(pars, 'recr_sd'))
}


# Fertility matrix -------------------------------------------------------------
fy <- function(y, x, pars, h, disturbance = 0) {
  Pvec <- recr_y(y, pars)
  
  recruit_mass <- sum(Pvec * h)
  
  if (!is.finite(recruit_mass) || recruit_mass <= 0) {
    stop('Recruit size distribution has no mass inside the mesh')
  }
  
  Pvec <- Pvec / recruit_mass
  
  Rvec <- rx(
    x = x,
    pars = pars,
    disturbance = disturbance)
  
  outer(Pvec, Rvec) * h
}


# Kernel -----------------------------------------------------------------------
kernel <- function(pars, disturbance = 0) {
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  # Fertility matrix
  Fmat <- fy(
    y = y,
    x = y,
    pars = pars,
    h = h,
    disturbance = disturbance)
  
  # Survival vector
  Smat <- sx(
    x = y,
    pars = pars,
    disturbance = disturbance)
  
  # Growth matrix
  Gmat <- matrix(0, n, n)
  
  Gmat[] <- t(
    outer(
      y,
      y,
      Vectorize(function(x, y) {
        gxy(
          x = x,
          y = y,
          pars = pars,
          disturbance = disturbance)
      }))) * h
  
  # Growth-survival matrix
  Tmat <- matrix(0, n, n)
  
  # Correct eviction of small individuals
  for (i in 1:(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }
  
  # Correct eviction of large individuals
  for (i in (n / 2 + 1):n) {
    Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }
  
  k_yx <- Fmat + Tmat
  
  list(
    k_yx = k_yx,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    Smat = Smat,
    meshpts = y,
    h = h)
}


# Population growth rate -------------------------------------------------------
lambda_ipm <- function(pars, disturbance = 0) {
  K <- kernel(
    pars = pars,
    disturbance = disturbance)$k_yx
  
  Re(eigen(K, only.values = TRUE)$values[1])
}


# Mean IPM ---------------------------------------------------------------------
lambda_mean_no_fire <- lambda_ipm(
  pars = pars_mean,
  disturbance = 0)

lambda_mean_fire <- lambda_ipm(
  pars = pars_mean,
  disturbance = 1)

lambda_mean_no_fire
lambda_mean_fire


# Years with complete year-specific parameter sets ----------------------------
ipm_years <- pars_year %>%
  filter(!is.na(year)) %>%
  distinct(year) %>%
  arrange(year) %>%
  pull(year)

ipm_years


# Year-specific IPM ------------------------------------------------------------
lambda_ipm_year <- function(year, disturbance = 0) {
  pars_i <- make_ipm_pars(year_i = year)
  
  lambda_ipm(
    pars = pars_i,
    disturbance = disturbance)
}


lambda_year <- tibble(
  year = ipm_years,
  lambda_no_fire = purrr::map_dbl(
    ipm_years,
    ~ lambda_ipm_year(.x, disturbance = 0)),
  lambda_fire = purrr::map_dbl(
    ipm_years,
    ~ lambda_ipm_year(.x, disturbance = 1)))

lambda_year


# Observed PGR vs asymptotic and projected lambda -----------------------------

# Site-year disturbance lookup ------------------------------------------------
disturbance_lookup_site <- df %>%
  filter(!is.na(site), !is.na(year)) %>%
  group_by(site, year) %>%
  summarise(
    disturbance = as.integer(any(disturbance == 1, na.rm = TRUE)),
    .groups = 'drop')


get_disturbance_site_year <- function(site_i, year_i) {
  disturbance_i <- disturbance_lookup_site %>%
    filter(site == site_i, year == year_i) %>%
    pull(disturbance)
  
  if (length(disturbance_i) == 0) {
    disturbance_i <- 0
  }
  
  disturbance_i[1]
}


# Observed site-year abundance -------------------------------------------------
df_counts_site <- df %>%
  filter(
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  group_by(site, year) %>%
  summarise(n = n(), .groups = 'drop')


# Pair abundance in year t with abundance in year t + 1 -----------------------
df_obs_pgr_site <- df_counts_site %>%
  rename(n_t0 = n) %>%
  left_join(
    df_counts_site %>%
      transmute(
        site,
        year = year - 1,
        n_t1 = n),
    by = c('site', 'year')) %>%
  filter(
    !is.na(n_t1),
    year %in% ipm_years) %>%
  left_join(
    disturbance_lookup_site,
    by = c('site', 'year')) %>%
  mutate(
    disturbance = replace_na(disturbance, 0),
    obs_pgr_site = n_t1 / n_t0)


# Initial observed size distribution ------------------------------------------
make_initial_n_site <- function(year_i, site_i, pars) {
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  df_i <- df %>%
    filter(
      year == year_i,
      site == site_i,
      size_t0 > 0,
      is.finite(logsize_t0))
  
  counts_i <- hist(
    pmin(pmax(df_i$logsize_t0, L), U),
    breaks = seq(L, U, length.out = n + 1),
    include.lowest = TRUE,
    plot = FALSE)$counts
  
  counts_i / h
}


# Project one site from year t to year t + 1 ----------------------------------
project_one_site_year <- function(year_i, site_i) {
  pars_i <- make_ipm_pars(year_i)
  
  disturbance_i <- get_disturbance_site_year(
    site_i = site_i,
    year_i = year_i)
  
  n_obs <- make_initial_n_site(
    year_i = year_i,
    site_i = site_i,
    pars = pars_i)
  
  kernel_i <- kernel(
    pars = pars_i,
    disturbance = disturbance_i)
  
  K <- kernel_i$k_yx
  h <- kernel_i$h
  
  n_proj <- K %*% n_obs
  
  n_obs_total <- sum(n_obs) * h
  n_proj_total <- sum(n_proj) * h
  
  tibble(
    year = year_i,
    site = site_i,
    disturbance = disturbance_i,
    n_obs_model = n_obs_total,
    n_proj_model = n_proj_total,
    asym_lambda = Re(
      eigen(K, only.values = TRUE)$values[1]),
    proj_lambda = n_proj_total / n_obs_total)
}


# Project all comparable site-years -------------------------------------------
df_proj_site <- bind_rows(
  lapply(seq_len(nrow(df_obs_pgr_site)), function(i) {
    project_one_site_year(
      year_i = df_obs_pgr_site$year[i],
      site_i = df_obs_pgr_site$site[i])
  }))


# Combine observed and projected site values ----------------------------------
df_compare_site <- df_obs_pgr_site %>%
  left_join(
    df_proj_site,
    by = c('year', 'site', 'disturbance'))


# Whole-population annual comparison ------------------------------------------
df_compare <- df_compare_site %>%
  group_by(year) %>%
  summarise(
    asym_lambda = weighted.mean(
      asym_lambda, w = n_obs_model, na.rm = TRUE),
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    obs_pgr = n_t1 / n_t0,
    n_obs_model = sum(n_obs_model, na.rm = TRUE),
    n_proj_model = sum(n_proj_model, na.rm = TRUE),
    proj_lambda = n_proj_model / n_obs_model,
    disturbance = if_else(
      any(disturbance == 1, na.rm = TRUE), 1, 0),
    .groups = 'drop') %>%
  mutate(
    disturbance = factor(
      disturbance, levels = c(0, 1),
      labels = c('No fire', 'Fire')))


# Long format for plotting -----------------------------------------------------
df_plot <- df_compare %>%
  select(
    year,
    obs_pgr,
    asym_lambda,
    proj_lambda,
    disturbance) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = 'lambda_type',
    values_to = 'lambda') %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = 'Asymptotic lambda',
      proj_lambda = 'Projected lambda'))


# Observed population growth versus modeled lambda ----------------------------
g_mod_vs_obs <- ggplot(
  df_plot,
  aes(
    x = lambda,
    y = obs_pgr,
    color = disturbance)) +
  geom_point(size = 2.5) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2) +
  facet_wrap(
    ~ lambda_type,
    scales = 'free_x') +
  scale_color_manual(
    values = c(
      'No fire' = 'black',
      'Fire' = 'red')) +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = v_ggp_suffix,
    x = expression('Modeled ' * lambda),
    y = 'Observed population growth rate',
    color = NULL) +
  theme_classic()

g_mod_vs_obs
