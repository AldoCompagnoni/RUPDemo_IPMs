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
# fig_gr
v_mod_set_su <- c()
# fig_su
v_mod_set_fl <- c()
# fig_fl
v_mod_set_fr <- c()
v_mod_set_fr2re_site <- c()



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
  


# Survival data ----------------------------------------------------------------
df_su <- df %>%
  filter(
    size_t0 != 0,
    !is.na(survives),
    is.finite(logsize_t0),
    !is.na(disturbance)
  ) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)
  ) %>%
  dplyr::select(
    plant_id, site, year, size_t0, survives, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    disturbance
  )


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(
    size_t0 != 0,
    survives == 1,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    !is.na(disturbance)
  ) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)
  ) %>%
  dplyr::select(
    plant_id, site, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    disturbance
  )


# Flower data ------------------------------------------------------------------
df_fl <- df %>%
  filter(
    !is.na(flower),
    is.finite(logsize_t0),
    !is.na(disturbance)
  ) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance),
    flower = if_else(flower > 0, 1, 0)
  ) %>%
  dplyr::select(
    plant_id, site, year, size_t0, flower, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    disturbance
  )


# Fruit data -------------------------------------------------------------------
# Biological interpretation:
# Flowering model estimates Pr(flowering | size, disturbance).
# Fruit model estimates E(fruits | flowering, size, disturbance).
# Therefore the IPM uses:
# expected fruits per plant = Pr(flowering) * E(fruits | flowering).

df_fr <- df %>%
  filter(
    flower > 0,
    is.finite(fruit),
    is.finite(logsize_t0),
    !is.na(disturbance)
  ) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)
  ) %>%
  dplyr::select(
    plant_id, site, year, size_t0, fruit,
    logsize_t0, logsize_t0_2, logsize_t0_3,
    disturbance
  )


# Recruitment data -------------------------------------------------------------
df_re <- df %>%
  group_by(year, site, quad_id) %>%
  summarise(
    disturbance = if_else(max(disturbance, na.rm = TRUE) == 1, 1, 0),
    tot_p_area = sum(size_t0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  {
    df_quad <- .
    
    df_group <- df_quad %>%
      group_by(year) %>%
      summarise(
        g_cov = mean(tot_p_area, na.rm = TRUE),
        .groups = "drop"
      )
    
    df_cover <- left_join(df_quad, df_group, by = "year") %>%
      mutate(year = as.integer(year + 1)) %>%
      drop_na()
    
    df_re <- df %>%
      group_by(year, site, quad_id) %>%
      summarise(
        nr_quad = sum(recruit == 1, na.rm = TRUE),
        .groups = "drop"
      )
    
    left_join(df_cover, df_re, by = c("year", "site", "quad_id"))
  } %>%
  mutate(
    year = factor(year),
    site = factor(site),
    disturbance = as.numeric(disturbance)
  )


# Removing year with too few data ----------------------------------------------
# Years original
v_years_og <- sort(unique(df$year))

df    <- df    %>% filter(!is.na(year) & !(year %in% v_years_re))
df_su <- df_su %>% filter(!is.na(year) & !(year %in% v_years_re))
df_gr <- df_gr %>% filter(!is.na(year) & !(year %in% v_years_re))
df_fl <- df_fl %>% filter(!is.na(year) & !(year %in% v_years_re))
df_fr <- df_fr %>% filter(!is.na(year) & !(year %in% v_years_re))
df_re <- df_re %>% filter(!is.na(year) & !(year %in% v_years_re))

# Years analysed
v_years      <- sort(unique(df$year))



# Survival model ---------------------------------------------------------------

ctrl_glmer <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

ctrl_lmer <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# GLMM; binomial
mod_su_0 <- glmer(
  survives ~ disturbance + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer
)

mod_su_1 <- glmer(
  survives ~ logsize_t0 + disturbance + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer
)

mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer
)

mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer
)

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

mod_su_best
summary(mod_su_best)
mods_su_dAIC


# Growth model -----------------------------------------------------------------

mod_gr_0 <- lmer(
  logsize_t1 ~ disturbance + (logsize_t0 | year),
  data = df_gr,
  control = ctrl_lmer
)

mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + disturbance + (logsize_t0 | year),
  data = df_gr,
  control = ctrl_lmer
)

mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance + (logsize_t0 | year),
  data = df_gr,
  control = ctrl_lmer
)

mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (logsize_t0 | year),
  data = df_gr,
  control = ctrl_lmer
)

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

mod_gr_best
summary(mod_gr_best)
mods_gr_dAIC


# Growth variance --------------------------------------------------------------

mod_gr_x <- fitted(mod_gr_best)
mod_gr_y <- resid(mod_gr_best)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE)
)


# Flower model -----------------------------------------------------------------

mod_fl_0 <- glmer(
  flower ~ disturbance + (1 | year),
  data = df_fl,
  family = binomial,
  control = ctrl_glmer
)

mod_fl_1 <- glmer(
  flower ~ logsize_t0 + disturbance + (1 | year),
  data = df_fl,
  family = binomial,
  control = ctrl_glmer
)

mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fl,
  family = binomial,
  control = ctrl_glmer
)

mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance + (1 | year),
  data = df_fl,
  family = binomial,
  control = ctrl_glmer
)

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

mod_fl_best
summary(mod_fl_best)
mods_fl_dAIC


# Fruit model ------------------------------------------------------------------
# Fruit is conditional on flowering because df_fr was filtered to flower > 0.
# Therefore the IPM will later use:
# Pr(flowering) * E(fruits | flowering).

mod_fr_0 <- glmer.nb(
  fruit ~ disturbance + (1 | year),
  data = df_fr,
  control = ctrl_glmer
)

mod_fr_1 <- glmer.nb(
  fruit ~ logsize_t0 + disturbance + (1 | year),
  data = df_fr,
  control = ctrl_glmer
)

mod_fr_2 <- glmer.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_fr,
  control = ctrl_glmer
)

# Do not add the cubic fruit model yet.
# Disturbed fruiting observations are sparse, so keep this conservative.

mods_fr <- list(mod_fr_0, mod_fr_1, mod_fr_2)
mods_fr_dAIC <- bbmle::AICctab(mods_fr, weights = TRUE, sort = FALSE)$dAIC
mods_fr_i_sort <- order(mods_fr_dAIC)

if (length(v_mod_set_fr) == 0) {
  mod_fr_i_best <- mods_fr_i_sort[1]
  v_mod_fr_i <- mod_fr_i_best - 1
} else {
  mod_fr_i_best <- v_mod_set_fr + 1
  v_mod_fr_i <- v_mod_set_fr
}

mod_fr_best <- mods_fr[[mod_fr_i_best]]
mod_fr_ranef <- coef(mod_fr_best)$year
mod_fr_ranef$year <- rownames(mod_fr_ranef)

mod_fr_best
summary(mod_fr_best)
mods_fr_dAIC


# Fruit to recruit: site-year transition ---------------------------------------
# Biological interpretation:
# Total fruits at site i in year t produce recruits at site i in year t + 1.
# This is the reproductive transition used by the IPM.
#
# Disturbance is not included directly here because disturbance already affects:
# survival, growth, flowering, and fruit production.
# Its effect on recruitment should therefore flow through fruit production.

df_fr2re_site_year <- df %>%
  mutate(
    site = factor(site),
    year = as.integer(as.character(year))
  ) %>%
  group_by(site, year) %>%
  summarise(
    fruits_site_t0 = sum(fruit, na.rm = TRUE),
    
    # recruit is coded NA / 1, so count true recruit records explicitly
    recruits_site_t0 = sum(recruit == 1, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(
    recruits_site_t1 = lead(recruits_site_t0),
    year_t1 = lead(year),
    year_gap = year_t1 - year
  ) %>%
  ungroup() %>%
  filter(year_gap == 1) %>%
  mutate(
    year_t0 = factor(year),
    year_t1 = factor(year_t1),
    site = factor(site),
    logfruits_site_t0 = log1p(fruits_site_t0),
    logrecruits_site_t1 = log1p(recruits_site_t1)
  )

df_fr2re_site_year %>%
  summarise(
    n_site_year_transitions = n(),
    n_sites = n_distinct(site),
    n_year_t1 = n_distinct(year_t1),
    n_positive_recruit_transitions = sum(recruits_site_t1 > 0),
    n_zero_recruit_transitions = sum(recruits_site_t1 == 0),
    total_fruits_t0 = sum(fruits_site_t0, na.rm = TRUE),
    total_recruits_t1 = sum(recruits_site_t1, na.rm = TRUE),
    max_fruits_site_t0 = max(fruits_site_t0, na.rm = TRUE),
    max_recruits_site_t1 = max(recruits_site_t1, na.rm = TRUE),
    proportion_zero_recruit = mean(recruits_site_t1 == 0, na.rm = TRUE)
  )


# Fruit-to-recruit candidate models -------------------------------------------

ctrl_re <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# 0. Fixed-effect model
fr2re_mod_00 <- lm(
  logrecruits_site_t1 ~ logfruits_site_t0,
  data = df_fr2re_site_year
)

# 1. Site random intercept
fr2re_mod_10 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 + (1 | site),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 2. Year random intercept
fr2re_mod_11 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 + (1 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 3. Site + year random intercepts
fr2re_mod_12 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 + (1 | site) + (1 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 4. Site random slope
fr2re_mod_20 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 + (logfruits_site_t0 | site),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 5. Year random slope
fr2re_mod_21 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 + (logfruits_site_t0 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 6. Site + year random slopes
fr2re_mod_22 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 +
    (logfruits_site_t0 | site) +
    (logfruits_site_t0 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 7. Site random slope only
fr2re_mod_30 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 +
    (0 + logfruits_site_t0 | site),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 8. Year random slope only
fr2re_mod_31 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 +
    (0 + logfruits_site_t0 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)

# 9. Site + year random slopes only
fr2re_mod_32 <- lmer(
  logrecruits_site_t1 ~ logfruits_site_t0 +
    (0 + logfruits_site_t0 | site) +
    (0 + logfruits_site_t0 | year_t1),
  data = df_fr2re_site_year,
  REML = FALSE,
  control = ctrl_re
)


fr2re_mods <- list(
  fr2re_mod_00,
  fr2re_mod_10,
  fr2re_mod_11,
  fr2re_mod_12,
  fr2re_mod_20,
  fr2re_mod_21,
  fr2re_mod_22,
  fr2re_mod_30,
  fr2re_mod_31,
  fr2re_mod_32
)

fr2re_mods_dAIC <- bbmle::AICctab(
  fr2re_mods,
  weights = TRUE,
  sort = FALSE
)$dAIC

fr2re_mods_i_sort <- order(fr2re_mods_dAIC)

if (length(v_mod_set_fr2re_site) == 0) {
  fr2re_mod_i_best <- fr2re_mods_i_sort[1]
  v_mod_fr2re_site_index <- fr2re_mod_i_best
} else {
  fr2re_mod_i_best <- v_mod_set_fr2re_site
  v_mod_fr2re_site_index <- v_mod_set_fr2re_site
}

fr2re_mod_best <- fr2re_mods[[fr2re_mod_i_best]]

fr2re_mod_best
summary(fr2re_mod_best)
fr2re_mods_dAIC


# Fruit-to-recruit convergence diagnostics -------------------------------------

fr2re_model_diagnostics <- tibble(
  model_id = seq_along(fr2re_mods),
  dAIC = as.numeric(fr2re_mods_dAIC),
  singular = purrr::map_lgl(
    fr2re_mods,
    ~ if (inherits(.x, "merMod")) isSingular(.x) else FALSE
  ),
  convergence_message = purrr::map_chr(
    fr2re_mods,
    function(m) {
      if (!inherits(m, "merMod")) return(NA_character_)
      msg <- summary(m)$optinfo$conv$lme4$messages
      if (is.null(msg)) return(NA_character_)
      paste(msg, collapse = "; ")
    }
  )
)

fr2re_model_diagnostics


# Extract fixed and random effects --------------------------------------------

get_fixef_re <- function(model) {
  if (inherits(model, "merMod")) {
    fixef(model)
  } else {
    coef(model)
  }
}

fr2re_fixef <- get_fixef_re(fr2re_mod_best)

get_ranef_re <- function(model, group_var, levels_vec) {
  
  if (!inherits(model, "merMod")) {
    out <- tibble(!!group_var := levels_vec)
  } else if (!group_var %in% names(ranef(model))) {
    out <- tibble(!!group_var := levels_vec)
  } else {
    out <- ranef(model)[[group_var]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column(group_var)
  }
  
  if (!"(Intercept)" %in% names(out)) {
    out[["(Intercept)"]] <- 0
  }
  
  if (!"logfruits_site_t0" %in% names(out)) {
    out[["logfruits_site_t0"]] <- 0
  }
  
  out %>%
    mutate(
      across(c(`(Intercept)`, logfruits_site_t0), ~ replace_na(.x, 0))
    )
}


fr2re_site_ranef <- get_ranef_re(
  model = fr2re_mod_best,
  group_var = "site",
  levels_vec = levels(df_fr2re_site_year$site)
) %>%
  rename(
    site_u0 = `(Intercept)`,
    site_u1 = logfruits_site_t0
  )

fr2re_year_ranef <- get_ranef_re(
  model = fr2re_mod_best,
  group_var = "year_t1",
  levels_vec = levels(df_fr2re_site_year$year_t1)
) %>%
  rename(
    year_u0 = `(Intercept)`,
    year_u1 = logfruits_site_t0
  )


df_fr2re_site_year_coef <- tidyr::expand_grid(
  site = levels(df_fr2re_site_year$site),
  year_t1 = levels(df_fr2re_site_year$year_t1)
) %>%
  mutate(
    site = factor(site, levels = levels(df_fr2re_site_year$site)),
    year_t1 = factor(year_t1, levels = levels(df_fr2re_site_year$year_t1))
  ) %>%
  left_join(fr2re_site_ranef, by = "site") %>%
  left_join(fr2re_year_ranef, by = "year_t1") %>%
  mutate(
    site_u0 = replace_na(site_u0, 0),
    site_u1 = replace_na(site_u1, 0),
    year_u0 = replace_na(year_u0, 0),
    year_u1 = replace_na(year_u1, 0),
    
    fr2re_b0_site_year =
      unname(fr2re_fixef["(Intercept)"]) + site_u0 + year_u0,
    
    fr2re_b1_site_year =
      unname(fr2re_fixef["logfruits_site_t0"]) + site_u1 + year_u1
  ) %>%
  select(
    site,
    year_t1,
    fr2re_b0_site_year,
    fr2re_b1_site_year,
    site_u0,
    site_u1,
    year_u0,
    year_u1
  )

df_fr2re_site_year_coef %>%
  print(n = 100)


# Backtransform predicted log recruits to recruit counts -----------------------

predict_site_year_recruits <- function(logfruits_site_t0, fr2re_b0, fr2re_b1) {
  
  pred_log_recruits <- fr2re_b0 + fr2re_b1 * logfruits_site_t0
  
  pred_recruits <- expm1(pred_log_recruits)
  
  pmax(pred_recruits, 0)
}


# Diagnostic prediction table and plot -----------------------------------------

df_fr2re_site_year_plot <- df_fr2re_site_year %>%
  left_join(
    df_fr2re_site_year_coef,
    by = c("site", "year_t1")
  ) %>%
  mutate(
    pred_recruits_site_t1 = predict_site_year_recruits(
      logfruits_site_t0 = logfruits_site_t0,
      fr2re_b0 = fr2re_b0_site_year,
      fr2re_b1 = fr2re_b1_site_year
    )
  )

df_fr2re_site_year_plot %>%
  select(
    site,
    year,
    year_t1,
    fruits_site_t0,
    recruits_site_t1,
    pred_recruits_site_t1,
    fr2re_b0_site_year,
    fr2re_b1_site_year
  ) %>%
  print(n = 100)


fig_fr2re_site_year <- ggplot(
  df_fr2re_site_year_plot,
  aes(x = fruits_site_t0, y = recruits_site_t1)
) +
  geom_jitter(
    width = 0,
    height = 0.08,
    alpha = 0.45,
    size = 1
  ) +
  geom_point(
    aes(y = pred_recruits_site_t1),
    colour = "red",
    alpha = 0.6,
    size = 1.2
  ) +
  facet_wrap(~ site, scales = "free_y") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1),
    name = "Total fruits at site in year t"
  ) +
  theme_bw() +
  labs(
    title = "Site-year fruit-to-recruit transition",
    subtitle = "Observed recruits and model predictions by site",
    y = "Total recruits at site in year t + 1"
  )

fig_fr2re_site_year


# Exporting parameter estimates ------------------------------------------------
# Goal:
# 1. pars_mean = fixed-effect / mean parameters
# 2. pars_year = year-specific vital-rate parameters
# 3. pars_site_year_recruit = site-year fruit-to-recruit parameters


# Helper functions -------------------------------------------------------------

empty_coef_df <- function() {
  data.frame(
    coefficient = character(),
    value = numeric()
  )
}


is_mixed_model <- function(model) {
  inherits(model, "merMod")
}


get_fixef_safe <- function(model) {
  if (is_mixed_model(model)) {
    lme4::fixef(model)
  } else {
    coef(model)
  }
}


get_group_coef_safe <- function(model, group_var) {
  
  if (!is_mixed_model(model)) {
    return(NULL)
  }
  
  coef_list <- coef(model)
  
  if (!(group_var %in% names(coef_list))) {
    return(NULL)
  }
  
  coef_list[[group_var]]
}


find_model_term <- function(model_terms, aliases) {
  
  aliases <- unlist(aliases)
  
  for (alias in aliases) {
    
    if (grepl("^regex:", alias)) {
      pattern <- sub("^regex:", "", alias)
      hit <- grep(pattern, model_terms, value = TRUE)
    } else {
      hit <- model_terms[model_terms == alias]
    }
    
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  NA_character_
}


bind_coef_rows <- function(x) {
  
  x <- x[!vapply(x, is.null, logical(1))]
  
  if (length(x) == 0) {
    return(empty_coef_df())
  }
  
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  
  out
}


check_duplicate_coefficients <- function(x, object_name = "parameter table") {
  
  dup <- x$coefficient[duplicated(x$coefficient)]
  
  if (length(dup) > 0) {
    stop(
      object_name,
      " contains duplicated coefficient names: ",
      paste(unique(dup), collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


coef_df_to_list <- function(x) {
  
  if (nrow(x) == 0) {
    return(list())
  }
  
  check_duplicate_coefficients(x)
  
  x %>%
    mutate(coefficient = as.character(coefficient)) %>%
    tidyr::pivot_wider(
      names_from = coefficient,
      values_from = value
    ) %>%
    as.list()
}


extract_fixed_pars <- function(model, term_map, fill_missing = TRUE) {
  
  model_coefs <- get_fixef_safe(model)
  model_terms <- names(model_coefs)
  
  out <- lapply(names(term_map), function(coef_name) {
    
    model_term <- find_model_term(
      model_terms = model_terms,
      aliases = term_map[[coef_name]]
    )
    
    if (is.na(model_term)) {
      
      if (fill_missing) {
        value <- 0
      } else {
        return(NULL)
      }
      
    } else {
      value <- unname(model_coefs[model_term])
    }
    
    data.frame(
      coefficient = coef_name,
      value = value
    )
  })
  
  bind_coef_rows(out)
}


extract_year_pars <- function(model, group_var, term_map, fill_missing = TRUE) {
  
  coef_matrix <- get_group_coef_safe(
    model = model,
    group_var = group_var
  )
  
  if (is.null(coef_matrix)) {
    return(empty_coef_df())
  }
  
  model_terms <- colnames(coef_matrix)
  
  out <- lapply(names(term_map), function(coef_prefix) {
    
    model_term <- find_model_term(
      model_terms = model_terms,
      aliases = term_map[[coef_prefix]]
    )
    
    if (is.na(model_term)) {
      
      if (fill_missing) {
        value <- rep(0, nrow(coef_matrix))
      } else {
        return(NULL)
      }
      
    } else {
      value <- coef_matrix[, model_term]
    }
    
    data.frame(
      coefficient = paste0(coef_prefix, "_", rownames(coef_matrix)),
      value = value
    )
  })
  
  bind_coef_rows(out)
}


# Term maps --------------------------------------------------------------------
# These maps allow the code to work whether the selected model contains
# intercept only, linear, quadratic, or cubic size terms.

size_term_map <- list(
  b0 = c("(Intercept)"),
  b1 = c("logsize_t0"),
  b2 = c("logsize_t0_2", "I(logsize_t0^2)"),
  b3 = c("logsize_t0_3", "I(logsize_t0^3)"),
  bf = c("disturbance", "regex:^disturbance")
)


make_size_term_map <- function(prefix) {
  
  out <- size_term_map
  names(out) <- paste0(prefix, names(out))
  
  out
}


su_term_map <- make_size_term_map("surv_")
gr_term_map <- make_size_term_map("grow_")
fl_term_map <- make_size_term_map("fl_")
fr_term_map <- make_size_term_map("fr_")


# Fixed-effect / mean parameters ----------------------------------------------

su_fe <- extract_fixed_pars(
  model = mod_su_best,
  term_map = su_term_map,
  fill_missing = TRUE
)

gr_fe <- extract_fixed_pars(
  model = mod_gr_best,
  term_map = gr_term_map,
  fill_missing = TRUE
)

fl_fe <- extract_fixed_pars(
  model = mod_fl_best,
  term_map = fl_term_map,
  fill_missing = TRUE
)

fr_fe <- extract_fixed_pars(
  model = mod_fr_best,
  term_map = fr_term_map,
  fill_missing = TRUE
)


# Recruit size distribution ----------------------------------------------------

df_re_size <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0)
  )

recr_sz_mean <- mean(df_re_size$logsize_t0, na.rm = TRUE)
recr_sz_sd <- sd(df_re_size$logsize_t0, na.rm = TRUE)


# Growth variance and IPM constants -------------------------------------------

gr_var_coef <- coef(mod_gr_var)

constants <- tibble::tribble(
  ~coefficient, ~value,
  "recr_sz", recr_sz_mean,
  "recr_sd", recr_sz_sd,
  "a", as.numeric(gr_var_coef[1]),
  "b", as.numeric(gr_var_coef[2]),
  "L", min(df_gr$logsize_t0, na.rm = TRUE),
  "U", max(df_gr$logsize_t0, na.rm = TRUE),
  "mat_siz", 200,
  "mod_su_index", v_mod_su_index,
  "mod_gr_index", v_mod_gr_index,
  "mod_fl_index", v_mod_fl_index,
  "mod_fr_index", v_mod_fr_i,
  "mod_fr2re_site_year_index", v_mod_fr2re_site_index
) %>%
  mutate(
    coefficient = as.character(coefficient),
    value = as.numeric(value)
  )


# Fruit-to-recruit fixed effects ----------------------------------------------

fr2re_fe <- tibble::tribble(
  ~coefficient, ~value,
  "fr2re_b0", unname(fr2re_fixef["(Intercept)"]),
  "fr2re_b1", unname(fr2re_fixef["logfruits_site_t0"])
) %>%
  mutate(
    coefficient = as.character(coefficient),
    value = as.numeric(value)
  )


# Mean parameter object --------------------------------------------------------

pars_cons <- bind_coef_rows(
  list(
    su_fe,
    gr_fe,
    fl_fe,
    fr_fe,
    fr2re_fe,
    constants
  )
)

check_duplicate_coefficients(pars_cons, object_name = "pars_cons")

pars_mean <- coef_df_to_list(pars_cons)

pars_mean


# Year-specific vital-rate parameters -----------------------------------------
# These are conditional coefficients from coef(model)$year.
# For example, if survival has a year random intercept, then surv_b0_2007
# is the year-specific intercept for 2007.

su_out_yr <- extract_year_pars(
  model = mod_su_best,
  group_var = "year",
  term_map = su_term_map,
  fill_missing = TRUE
)

gr_out_yr <- extract_year_pars(
  model = mod_gr_best,
  group_var = "year",
  term_map = gr_term_map,
  fill_missing = TRUE
)

fl_out_yr <- extract_year_pars(
  model = mod_fl_best,
  group_var = "year",
  term_map = fl_term_map,
  fill_missing = TRUE
)

fr_out_yr <- extract_year_pars(
  model = mod_fr_best,
  group_var = "year",
  term_map = fr_term_map,
  fill_missing = TRUE
)

pars_var <- bind_coef_rows(
  list(
    su_out_yr,
    gr_out_yr,
    fl_out_yr,
    fr_out_yr
  )
)

check_duplicate_coefficients(pars_var, object_name = "pars_var")

pars_var_wide <- coef_df_to_list(pars_var)

pars_var_wide


# Build one vital-rate parameter object per year -------------------------------

make_year_pars <- function(year_i, pars_mean, pars_var_wide) {
  
  p <- pars_mean
  
  yr <- as.character(year_i)
  
  possible_names <- names(pars_var_wide)
  
  for (nm in possible_names) {
    
    suffix <- paste0("_", yr)
    
    if (endsWith(nm, suffix)) {
      
      base_nm <- sub(paste0("_", yr, "$"), "", nm)
      p[[base_nm]] <- pars_var_wide[[nm]]
    }
  }
  
  p$year <- yr
  
  p
}


v_years_model <- sort(unique(as.character(df$year)))

pars_year <- tibble(
  year = v_years_model,
  pars_year = purrr::map(
    year,
    ~ make_year_pars(
      year_i = .x,
      pars_mean = pars_mean,
      pars_var_wide = pars_var_wide
    )
  )
)

pars_year


# Site-year fruit-to-recruit parameter table -----------------------------------
# This uses the selected fruit-to-recruit model.
# In your current run, model 10 was selected:
#   logrecruits ~ logfruits + (0 + logfruits | site) +
#                              (0 + logfruits | year_t1)
# but this table also works if a simpler model is selected.

pars_site_year_recruit <- df_fr2re_site_year_coef %>%
  mutate(
    site = factor(site),
    year_t1 = factor(year_t1),
    year = as.character(as.integer(as.character(year_t1)) - 1),
    
    fr2re_b0 = fr2re_b0_site_year,
    fr2re_b1 = fr2re_b1_site_year
  ) %>%
  select(
    site,
    year,
    year_t1,
    fr2re_b0,
    fr2re_b1,
    site_u0,
    site_u1,
    year_u0,
    year_u1
  )

pars_site_year_recruit %>%
  print(n = 100)


# Building the year- and site-specific IPMs ------------------------------------
# These functions use:
# 1. year-specific survival, growth, flowering, and fruiting parameters
# 2. site-year fruit-to-recruit parameters
# 3. disturbance as a state applied to the vital-rate intercepts


# Core functions ---------------------------------------------------------------

inv_logit <- function(x) {
  plogis(x)
}


apply_disturbance_state <- function(p, disturbance_value) {
  
  p2 <- p
  
  if (!is.null(p$surv_bf)) {
    p2$surv_b0 <- p$surv_b0 + p$surv_bf * disturbance_value
  }
  
  if (!is.null(p$grow_bf)) {
    p2$grow_b0 <- p$grow_b0 + p$grow_bf * disturbance_value
  }
  
  if (!is.null(p$fl_bf)) {
    p2$fl_b0 <- p$fl_b0 + p$fl_bf * disturbance_value
  }
  
  if (!is.null(p$fr_bf)) {
    p2$fr_b0 <- p$fr_b0 + p$fr_bf * disturbance_value
  }
  
  p2
}


# Survival ---------------------------------------------------------------------

sx <- function(x, pars, num_pars = pars$mod_su_index) {
  
  survival_value <- pars$surv_b0
  
  for (i in seq_len(num_pars)) {
    param_name <- paste0("surv_b", i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^i
    }
  }
  
  inv_logit(survival_value)
}


# Growth -----------------------------------------------------------------------

grow_mu <- function(x, pars, num_pars = pars$mod_gr_index) {
  
  mean_value <- 0
  
  for (i in 0:num_pars) {
    param_name <- paste0("grow_b", i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  
  mean_value
}


grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * x))
}


gxy <- function(x, y, pars) {
  
  dnorm(
    y,
    mean = grow_mu(x, pars),
    sd = grow_sd(x, pars)
  )
}


pxy <- function(x, y, pars) {
  sx(x, pars) * gxy(x, y, pars)
}


# Flowering --------------------------------------------------------------------

fl_x <- function(x, pars, num_pars = pars$mod_fl_index) {
  
  flower_value <- pars$fl_b0
  
  for (i in seq_len(num_pars)) {
    param_name <- paste0("fl_b", i)
    if (!is.null(pars[[param_name]])) {
      flower_value <- flower_value + pars[[param_name]] * x^i
    }
  }
  
  inv_logit(flower_value)
}


# Fruiting conditional on flowering --------------------------------------------

fr_x <- function(x, pars, num_pars = pars$mod_fr_index) {
  
  fruit_value <- pars$fr_b0
  
  for (i in seq_len(num_pars)) {
    param_name <- paste0("fr_b", i)
    if (!is.null(pars[[param_name]])) {
      fruit_value <- fruit_value + pars[[param_name]] * x^i
    }
  }
  
  exp(fruit_value)
}


# Expected fruits per plant ----------------------------------------------------
# Because the fruit model was fit only to flowering plants:
# expected fruits per plant = Pr(flowering) * E(fruits | flowering)

fruit_x <- function(x, pars) {
  fl_x(x, pars) * fr_x(x, pars)
}


# Recruitment size distribution ------------------------------------------------
# Normalize over the IPM bounds so recruitment density integrates to one
# within [L, U].

re_y_dist <- function(y, pars) {
  
  dens <- dnorm(
    y,
    mean = pars$recr_sz,
    sd = pars$recr_sd
  )
  
  norm_const <- pnorm(
    pars$U,
    mean = pars$recr_sz,
    sd = pars$recr_sd
  ) -
    pnorm(
      pars$L,
      mean = pars$recr_sz,
      sd = pars$recr_sd
    )
  
  dens / norm_const
}


# Fruit-to-recruit backtransformation ------------------------------------------

predict_site_year_recruits_from_fruits <- function(site_fruits, pars) {
  
  if (is.na(site_fruits) || site_fruits <= 0) {
    return(0)
  }
  
  pred_log_recruits <- pars$fr2re_b0 +
    pars$fr2re_b1 * log1p(site_fruits)
  
  pred_recruits <- expm1(pred_log_recruits)
  
  pmax(pred_recruits, 0)
}


site_year_recruits_per_fruit <- function(pars) {
  
  site_fruits <- pars$site_fruits_ref
  
  if (is.null(site_fruits) || is.na(site_fruits) || site_fruits <= 0) {
    return(0)
  }
  
  site_recruits <- predict_site_year_recruits_from_fruits(
    site_fruits = site_fruits,
    pars = pars
  )
  
  site_recruits / site_fruits
}


# Fertility kernel --------------------------------------------------------------
# This replaces the old:
# fl_x * fr_x * fecu_b0 * recruit_size_distribution
#
# with:
# expected fruits per plant *
# site-year recruits per fruit *
# recruit size distribution

fyx <- function(y, x, pars) {
  
  fruit_x(x, pars) *
    site_year_recruits_per_fruit(pars) *
    re_y_dist(y, pars)
}


# Kernel -----------------------------------------------------------------------

kernel <- function(pars) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  # Survival vector
  Smat <- sx(y, pars)
  
  # Growth matrix
  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Growth/survival transition matrix
  Tmat <- matrix(0, n, n)
  
  # Correct eviction of small plants
  for (i in 1:(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for (i in (n / 2 + 1):n) {
    Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }
  
  # Fertility matrix
  Fmat <- outer(
    y, y,
    Vectorize(function(x, y) fyx(x, y, pars))
  ) * h
  
  k_yx <- Fmat + Tmat
  
  list(
    k_yx = k_yx,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    meshpts = y
  )
}


lambda_ipm <- function(pars) {
  Re(eigen(kernel(pars)$k_yx)$values[1])
}


# Site-year projection ---------------------------------------------------------
# One kernel per site-year transition.
#
# Vital rates:
#   survival, growth, flowering, fruiting = year-specific
#
# Recruitment:
#   site-year fruit-to-recruit model
#
# Disturbance:
#   kernel mixture:
#   K = (1 - p_disturbance) * K_no_disturbance +
#       p_disturbance       * K_disturbance


# Observed disturbance exposure by site-year -----------------------------------

p_disturbance_site_year <- df %>%
  filter(
    !is.na(logsize_t0),
    is.finite(logsize_t0),
    !is.na(disturbance)
  ) %>%
  mutate(
    site = factor(site),
    year = as.character(year)
  ) %>%
  group_by(site, year) %>%
  summarise(
    p_disturbance = mean(disturbance, na.rm = TRUE),
    .groups = "drop"
  )

p_disturbance_site_year %>%
  print(n = 100)


# Initial observed size distribution for one site-year -------------------------

make_initial_n_site_year <- function(pars, site_i, year_i, df_init = df) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  breaks <- seq(L, U, length.out = n + 1)
  
  x <- df_init %>%
    filter(
      site == site_i,
      as.character(year) == as.character(year_i),
      !is.na(logsize_t0),
      is.finite(logsize_t0)
    ) %>%
    pull(logsize_t0)
  
  if (length(x) == 0) {
    return(rep(0, n))
  }
  
  adult_counts <- hist(
    pmin(pmax(x, L), U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE
  )$counts
  
  adult_density <- adult_counts / h
  
  adult_density
}


# Expected total fruits from one site's observed size distribution -------------

expected_site_year_fruits <- function(pars, site_i, year_i, df_init = df) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  n_site_year <- make_initial_n_site_year(
    pars = pars,
    site_i = site_i,
    year_i = year_i,
    df_init = df_init
  )
  
  sum(fruit_x(y, pars) * n_site_year * h, na.rm = TRUE)
}


# Get one year-specific parameter object ---------------------------------------

get_year_pars <- function(year_i) {
  
  out <- pars_year %>%
    filter(year == as.character(year_i)) %>%
    pull(pars_year)
  
  if (length(out) != 1) {
    stop("Could not find exactly one parameter object for year: ", year_i)
  }
  
  out[[1]]
}


# Project one site-year kernel -------------------------------------------------

project_site_year_ipm <- function(
    site_i,
    year_i,
    year_t1_i,
    p_disturbance_i,
    fr2re_b0_i,
    fr2re_b1_i,
    df_init = df
) {
  
  # Year-specific vital-rate parameters
  p_base <- get_year_pars(year_i)
  
  p_base$site <- site_i
  p_base$year <- as.character(year_i)
  p_base$year_t1 <- as.character(year_t1_i)
  p_base$p_disturbance <- p_disturbance_i
  
  # No-disturbance and disturbance state kernels
  p_no <- apply_disturbance_state(
    p = p_base,
    disturbance_value = 0
  )
  
  p_dist <- apply_disturbance_state(
    p = p_base,
    disturbance_value = 1
  )
  
  # Expected fruits under the site-year disturbance regime
  site_fruits_no <- expected_site_year_fruits(
    pars = p_no,
    site_i = site_i,
    year_i = year_i,
    df_init = df_init
  )
  
  site_fruits_dist <- expected_site_year_fruits(
    pars = p_dist,
    site_i = site_i,
    year_i = year_i,
    df_init = df_init
  )
  
  site_fruits_ref <-
    (1 - p_disturbance_i) * site_fruits_no +
    p_disturbance_i       * site_fruits_dist
  
  # Add site-year fruit-to-recruit parameters
  p_no$fr2re_b0 <- fr2re_b0_i
  p_no$fr2re_b1 <- fr2re_b1_i
  p_no$site_fruits_ref <- site_fruits_ref
  
  p_dist$fr2re_b0 <- fr2re_b0_i
  p_dist$fr2re_b1 <- fr2re_b1_i
  p_dist$site_fruits_ref <- site_fruits_ref
  
  # Observed initial size distribution
  n_obs <- make_initial_n_site_year(
    pars = p_base,
    site_i = site_i,
    year_i = year_i,
    df_init = df_init
  )
  
  h <- (p_base$U - p_base$L) / p_base$mat_siz
  
  n_initial <- sum(n_obs) * h
  
  if (n_initial <= 0) {
    return(
      tibble(
        lambda_asymptotic = NA_real_,
        lambda_projected = NA_real_,
        n_initial = n_initial,
        n_projected = NA_real_,
        p_disturbance = p_disturbance_i,
        site_fruits_no = site_fruits_no,
        site_fruits_dist = site_fruits_dist,
        site_fruits_ref = site_fruits_ref,
        site_recruits_ref = NA_real_,
        recruits_per_fruit_ref = NA_real_
      )
    )
  }
  
  K_no <- kernel(p_no)$k_yx
  K_dist <- kernel(p_dist)$k_yx
  
  K_regime <-
    (1 - p_disturbance_i) * K_no +
    p_disturbance_i       * K_dist
  
  n_proj <- K_regime %*% n_obs
  
  site_recruits_ref <- predict_site_year_recruits_from_fruits(
    site_fruits = site_fruits_ref,
    pars = p_no
  )
  
  recruits_per_fruit_ref <- ifelse(
    site_fruits_ref > 0,
    site_recruits_ref / site_fruits_ref,
    0
  )
  
  tibble(
    lambda_asymptotic = Re(eigen(K_regime)$values[1]),
    lambda_projected = as.numeric((sum(n_proj) * h) / n_initial),
    n_initial = n_initial,
    n_projected = sum(n_proj) * h,
    p_disturbance = p_disturbance_i,
    site_fruits_no = site_fruits_no,
    site_fruits_dist = site_fruits_dist,
    site_fruits_ref = site_fruits_ref,
    site_recruits_ref = site_recruits_ref,
    recruits_per_fruit_ref = recruits_per_fruit_ref
  )
}


# Projection grid --------------------------------------------------------------

df_site_year_projection_grid <- pars_site_year_recruit %>%
  mutate(
    site = factor(site),
    year = as.character(year),
    year_t1 = as.character(year_t1)
  ) %>%
  left_join(
    p_disturbance_site_year,
    by = c("site", "year")
  ) %>%
  mutate(
    p_disturbance = replace_na(p_disturbance, 0)
  )


# Run all site-year IPMs --------------------------------------------------------

df_lambda_site_year <- df_site_year_projection_grid %>%
  mutate(
    lambda_data = purrr::pmap(
      list(
        site,
        year,
        year_t1,
        p_disturbance,
        fr2re_b0,
        fr2re_b1
      ),
      function(site, year, year_t1, p_disturbance, fr2re_b0, fr2re_b1) {
        project_site_year_ipm(
          site_i = site,
          year_i = year,
          year_t1_i = year_t1,
          p_disturbance_i = p_disturbance,
          fr2re_b0_i = fr2re_b0,
          fr2re_b1_i = fr2re_b1,
          df_init = df
        )
      }
    )
  ) %>%
  select(
    site,
    year,
    year_t1,
    fr2re_b0,
    fr2re_b1,
    lambda_data
  ) %>%
  tidyr::unnest(lambda_data)

df_lambda_site_year %>%
  print(n = 100)


# Observed site-year population growth ----------------------------------------
# Compare the model-projected site-year lambda to observed change in the number
# of sized plants from year t to year t + 1.

df_obs_counts_site_year <- df %>%
  filter(
    !is.na(logsize_t0),
    is.finite(logsize_t0)
  ) %>%
  mutate(
    site = factor(site),
    year = as.character(year)
  ) %>%
  group_by(site, year) %>%
  summarise(
    n_obs = n_distinct(plant_id),
    .groups = "drop"
  ) %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(
    year_t1 = lead(year),
    n_obs_t1 = lead(n_obs),
    year_gap = as.integer(year_t1) - as.integer(year),
    lambda_obs = n_obs_t1 / n_obs
  ) %>%
  ungroup() %>%
  filter(
    year_gap == 1,
    n_obs > 0,
    n_obs_t1 > 0
  )

df_obs_counts_site_year %>%
  print(n = 100)

# Compare modeled and observed site-year lambda --------------------------------

df_lambda_site_year_compare <- df_lambda_site_year %>%
  left_join(
    df_obs_counts_site_year %>%
      select(site, year, year_t1, n_obs, n_obs_t1, lambda_obs),
    by = c("site", "year", "year_t1")
  ) %>%
  mutate(
    error_asymptotic_vs_obs = lambda_asymptotic - lambda_obs,
    error_projected_vs_obs = lambda_projected - lambda_obs
  )

df_lambda_site_year_compare %>%
  select(
    site,
    year,
    year_t1,
    lambda_asymptotic,
    lambda_projected,
    lambda_obs,
    n_initial,
    n_projected,
    n_obs,
    n_obs_t1,
    p_disturbance,
    error_projected_vs_obs
  ) %>%
  print(n = 100)

# Summary of modeled versus observed growth ------------------------------------

df_lambda_site_year_compare_summary <- df_lambda_site_year_compare %>%
  summarise(
    n_rows = n(),
    n_with_observed_lambda = sum(!is.na(lambda_obs)),
    n_missing_modeled_projected = sum(is.na(lambda_projected)),
    n_missing_observed = sum(is.na(lambda_obs)),
    
    mean_lambda_obs = mean(lambda_obs, na.rm = TRUE),
    geometric_lambda_obs = exp(mean(log(lambda_obs), na.rm = TRUE)),
    
    mean_lambda_asymptotic = mean(lambda_asymptotic, na.rm = TRUE),
    mean_lambda_projected = mean(lambda_projected, na.rm = TRUE),
    
    mean_error_asymptotic_vs_obs = mean(error_asymptotic_vs_obs, na.rm = TRUE),
    mean_error_projected_vs_obs = mean(error_projected_vs_obs, na.rm = TRUE),
    
    rmse_asymptotic_vs_obs = sqrt(mean(error_asymptotic_vs_obs^2, na.rm = TRUE)),
    rmse_projected_vs_obs = sqrt(mean(error_projected_vs_obs^2, na.rm = TRUE))
  )

df_lambda_site_year_compare_summary


# Plot observed versus modeled site-year lambda --------------------------------

df_lambda_site_year_compare_plot <- df_lambda_site_year_compare %>%
  select(
    site,
    year,
    year_t1,
    lambda_obs,
    lambda_asymptotic,
    lambda_projected
  ) %>%
  filter(!is.na(lambda_obs)) %>%
  pivot_longer(
    cols = c(lambda_asymptotic, lambda_projected),
    names_to = "modeled_lambda_type",
    values_to = "modeled_lambda"
  ) %>%
  mutate(
    modeled_lambda_type = recode(
      modeled_lambda_type,
      lambda_asymptotic = "Asymptotic lambda",
      lambda_projected = "Projected lambda"
    )
  )

fig_lambda_site_year_compare <- ggplot(
  df_lambda_site_year_compare_plot,
  aes(x = modeled_lambda, y = lambda_obs)
) +
  geom_point(alpha = 0.65, size = 1.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~ modeled_lambda_type) +
  theme_bw() +
  labs(
    title = "Observed versus modeled site-year population growth",
    subtitle = v_ggp_suffix,
    x = expression("Modeled " * lambda),
    y = expression("Observed " * lambda)
  )

fig_lambda_site_year_compare



# Site-level summary -----------------------------------------------------------

df_lambda_site_summary <- df_lambda_site_year_complete %>%
  group_by(site) %>%
  summarise(
    n_year_transitions = n(),
    
    lambda_obs_geometric = exp(mean(log(lambda_obs), na.rm = TRUE)),
    lambda_obs_arithmetic = mean(lambda_obs, na.rm = TRUE),
    
    lambda_asymptotic_geometric = exp(mean(log(lambda_asymptotic), na.rm = TRUE)),
    lambda_asymptotic_arithmetic = mean(lambda_asymptotic, na.rm = TRUE),
    
    lambda_projected_geometric = exp(mean(log(lambda_projected), na.rm = TRUE)),
    lambda_projected_arithmetic = mean(lambda_projected, na.rm = TRUE),
    
    error_projected_geo_vs_obs_geo =
      lambda_projected_geometric - lambda_obs_geometric,
    
    rmse_projected_vs_obs =
      sqrt(mean((lambda_projected - lambda_obs)^2, na.rm = TRUE)),
    
    mean_n_initial = mean(n_initial, na.rm = TRUE),
    mean_p_disturbance = mean(p_disturbance, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(as.numeric(as.character(site)))

df_lambda_site_summary %>%
  print(n = 100, width = Inf)