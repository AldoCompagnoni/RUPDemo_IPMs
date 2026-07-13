# IPM by site with disturbance - Archbold -  - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.01 (y-m-d)


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
  # bbmle is for AICctab
  bbmle,
  # patchwork plot alingment
  patchwork,
  # binom.cofint for the survival plot
  binom,
  # mixed models
  glmmTMB,
  skimr,
  scales,
  lubridate, 
  lme4) # , skimr, ipmr, binom, janitor, lme4


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Crotalaria avonensis')
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
# function to plot your survival data 'binned' (instead of 'jittered')
source('helper_functions/plot_binned_prop.R')
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
  filter(size_t0 != 0) %>%
  mutate(dist_label = ifelse(disturbance == 1, 'Disturbance', 'No disturbance')) %>%
  select(plant_id, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         disturbance, dist_label)

# Create binned summary using group_split()
df_su_binned <- df_su %>%
  group_split(dist_label) %>%
  purrr::map_df(~ plot_binned_prop(.x, 10, logsize_t0, survives) %>%
                  mutate(dist_label = unique(.x$dist_label)))

fig_su_overall <- ggplot(df_su_binned, aes(x = logsize_t0, y = survives, color = dist_label)) +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives, color = dist_label),
              position = position_jitter(width = 0.1, height = 0.3), alpha = 0.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c('No disturbance' = 'black', 'Disturbance' = 'red')) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = 'top') +
  labs(
    title = 'Survival',
    subtitle = v_ggp_suffix,
    x = expression('log(size)'[t0]),
    y = 'Survival Probability')

fig_su_overall


# Survival model ---------------------------------------------------------------
# Logistic regression
mod_su_0 <- glm(survives ~ disturbance,
                data = df_su, family = 'binomial')
# Logistic regression
mod_su_1 <- glm(survives ~ logsize_t0 + disturbance,
                data = df_su, family = 'binomial')
# Quadratic logistic model
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2 + disturbance,
                data = df_su, family = 'binomial')
# Cubic logistic model
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
                data = df_su, family = 'binomial')


# Compare models using AICc
mods_su      <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAICc <- AICctab(mods_su, weights = T, sort = F)$dAICc

# Get the sorted indices of dAICc values
mods_su_sorted <- order(mods_su_dAICc)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index       <- mod_su_index_bestfit - 1
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index       <- v_mod_set_su
}

mod_su_bestfit <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef   <- coef(mod_su_bestfit)

# Create prediction data frames for disturbance = 0 and disturbance = 1
df_su_pred <- data.frame(
  logsize_t0 = rep(seq(min(df_su$logsize_t0), max(df_su$logsize_t0), length.out = 100), 2),
  disturbance = rep(c(0, 1), each = 100)) %>%
  mutate(logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3) %>% 
  mutate(dist_label = ifelse(disturbance == 1, 'Disturbance', 'No disturbance'),
         survives = predict(mod_su_bestfit, newdata = ., type = 'response'))

# Binned observed data for both disturbance levels
df_su_binned <- bind_rows(
  plot_binned_prop(filter(df_su, disturbance == 0), 10, logsize_t0, survives) %>%
    mutate(dist_label = 'No disturbance'),
  plot_binned_prop(filter(df_su, disturbance == 1), 10, logsize_t0, survives) %>%
    mutate(dist_label = 'Disturbance'))

# Plot 1: Raw data + prediction lines
fig_su_line_combined <- ggplot() +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives, color = dist_label),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives, color = dist_label),
            linewidth = 0.9) +
  scale_color_manual(values = c('No disturbance' = 'black', 'Disturbance' = 'red')) +
  theme_bw() +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Survival to time t1') +
  theme(legend.position = 'none')

# Plot 2: Binned proportions + prediction lines
fig_su_bin_combined <- ggplot() +
  geom_point(data = df_su_binned, aes(x = logsize_t0, y = survives, color = dist_label)) +
  geom_errorbar(data = df_su_binned, aes(x = logsize_t0, ymin = lwr, ymax = upr, color = dist_label),
                width = 0.2) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives, color = dist_label),
            linewidth = 0.9) +
  scale_color_manual(values = c('No disturbance' = 'black', 'Disturbance' = 'red')) +
  theme_bw() +
  ylim(0, 1) +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Survival to time t1') +
  theme(legend.title = element_blank(), legend.position = 'top')

# Combine the two plots
fig_su_all <- fig_su_line_combined + fig_su_bin_combined +
  plot_annotation(
    title = 'Survival',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = 'bold'),
      plot.subtitle = element_text(size = 10, face = 'italic')))

fig_su_all
# ggsave(file.path(dir_result, 'survival_by_size.png'),
#        plot = fig_su_all, width = 10, height = 5, dpi = 300)


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(size_t0 != 0,
         is.finite(logsize_t1)) %>%
  mutate(dist_label = ifelse(disturbance == 1, 'disturbance', 'No disturbance')) %>%
  select(plant_id, year, size_t0, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         disturbance, dist_label)

fig_gr_overall <- ggplot(df_gr, aes(x = logsize_t0, y = logsize_t1, color = dist_label)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 8),
    title           = element_text(size = 10),
    plot.subtitle   = element_text(size = 8),
    legend.title    = element_blank(),
    legend.position = 'top') +
  labs(
    title    = 'Growth',
    subtitle = v_ggp_suffix,
    x        = expression('log(size)'[t0]),
    y        = expression('log(size)'[t1]))

fig_gr_overall


# Growth model -----------------------------------------------------------------
# Intercept model
mod_gr_0 <- lm(logsize_t1 ~ disturbance,
               data = df_gr, na.action = na.omit)
# Linear model
mod_gr_1 <- lm(logsize_t1 ~ logsize_t0 + disturbance,
               data = df_gr)
# Quadratic model
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance,
               data = df_gr)
# Cubic model
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
               data = df_gr)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAICc <- AICctab(mods_gr, weights = T, sort = F)$dAICc

# Get the sorted indices of dAICc values
mods_gr_sorted <- order(mods_gr_dAICc)

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

# Create prediction data for disturbance = 0 and disturbance = 1
df_gr_pred <- data.frame(
  logsize_t0 = rep(seq(min(df_gr$logsize_t0, na.rm = TRUE),
                       max(df_gr$logsize_t0, na.rm = TRUE), length.out = 100), 2),
  disturbance = rep(c(0, 1), each = 100))

df_gr_pred <- df_gr_pred %>%
  mutate(logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3,
         dist_label   = ifelse(disturbance == 1, 'disturbance', 'No disturbance'))

df_gr_pred$logsize_t1 <- predict(mod_gr_bestfit, newdata = df_gr_pred)

# Plot 1: Observed points + prediction lines
fig_gr_line_combined <- ggplot(df_gr, aes(x = logsize_t0, y = logsize_t1, color = dist_label)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_line(aes(x = logsize_t0, y = logsize_t1, color = dist_label),
            data = df_gr_pred, linewidth = 1) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  labs(x = expression('log(size)'[t0]),
       y = expression('log(size)'[t1])) +
  theme(
    legend.position = 'top',
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank())

# Plot 2: Predicted vs observed
fig_gr_pred_combined <- ggplot(df_gr, aes(x = predict(mod_gr_bestfit, newdata = df_gr),
                                          y = logsize_t1, color = dist_label)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', linewidth = 1) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  labs(x = 'Predicted', y = 'Observed') +
  theme(
    legend.position = 'top',
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank())

# Combine plots
fig_gr_all_combined <- fig_gr_line_combined + fig_gr_pred_combined +
  plot_annotation(
    title = 'Growth Prediction',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = 'bold'),
      plot.subtitle = element_text(size = 10, face = 'italic')))

fig_gr_all_combined
# ggsave(file.path(dir_result, 'growth_by_size.png'),
#        plot = fig_gr_all_combined, width = 10, height = 5, dpi = 300)


# Growth variance --------------------------------------------------------------
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Flowering data ----------------------------------------------------------------
df_fl <- df %>%
  filter(!is.na(flower), is.finite(logsize_t0), 
         !is.na(disturbance) & !is.na(size_t0) & !is.na(survives)) %>%
  select(plant_id, year, size_t0, flower, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, disturbance) %>%
  mutate(
    flower = if_else(flower > 0, 1, flower),
    dist_label = if_else(disturbance == 1, 'disturbance', 'No disturbance'))

# Create binned flowering probability by disturbance group
df_fl_binned <- df_fl %>%
  filter(!is.na(disturbance) & !is.na(size_t0)) %>%
  group_split(dist_label) %>%
  map_df(~ plot_binned_prop(.x, 10, logsize_t0, flower) %>%
           mutate(dist_label = unique(.x$dist_label)))

# Plot overlapped flowering probability
fig_fl_overall <- ggplot(df_fl_binned, aes(x = logsize_t0, y = flower, color = dist_label)) +
  geom_jitter(
    data = df_fl,
    aes(x = logsize_t0, y = flower, color = dist_label),
    position = position_jitter(width = 0.1, height = 0.3),
    alpha = 0.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = 'top') +
  labs(
    title = 'Flowering probability by disturbance status',
    subtitle = v_ggp_suffix,
    x = expression('log(size)'[t0]),
    y = 'Flowering Probability')

fig_fl_overall


# Flower model -----------------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ disturbance,
                data = df_fl, family = 'binomial')
# Logistic regression
mod_fl_1 <- glm(flower ~ logsize_t0 + disturbance,
                data = df_fl, family = 'binomial')
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logsize_t0 + logsize_t0_2 + disturbance,
                data = df_fl, family = 'binomial')
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
                data = df_fl, family = 'binomial')

# Compare models using AICc
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAICc <- AICctab(mods_fl, weights = T, sort = F)$dAICc

# Get the sorted indices of dAICc values
mods_fl_sorted <- order(mods_fl_dAICc)

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

# Create prediction data frames for disturbance = 0 and disturbance = 1
df_fl_pred_0 <- data.frame(
  logsize_t0 = seq(
    min(df_fl$logsize_t0, na.rm = T), max(df_fl$logsize_t0, na.rm = T), length.out = 100),
  disturbance = 0) %>% mutate(logsize_t0_2 = logsize_t0^2)

df_fl_pred_1 <- data.frame(
  logsize_t0 = seq(
    min(df_fl$logsize_t0, na.rm = T), max(df_fl$logsize_t0, na.rm = T), length.out = 100),
  disturbance = 1) %>% mutate(logsize_t0_2 = logsize_t0^2)

# Add predictions using model
df_fl_pred <- data.frame(
  logsize_t0 = rep(seq(min(df_fl$logsize_t0, na.rm = TRUE), max(df_fl$logsize_t0, na.rm = TRUE), length.out = 100), 2),
  disturbance = rep(c(0, 1), each = 100)) %>%
  mutate(
    dist_label = ifelse(disturbance == 1, 'disturbance', 'No disturbance'),
    logsize_t0_2 = logsize_t0^2)

df_fl_pred$flower <- predict(mod_fl_bestfit, newdata = df_fl_pred, type = 'response')

# Binned observed data for both disturbance levels
df_fl_binned <- bind_rows(
  plot_binned_prop(filter(df_fl, disturbance == 0), 10, logsize_t0, flower) %>%
    mutate(dist_label = 'No disturbance'),
  plot_binned_prop(filter(df_fl, disturbance == 1), 10, logsize_t0, flower) %>%
    mutate(dist_label = 'disturbance'))

# Plot 1: Raw jitter + prediction lines
fig_fl_line_combined <- ggplot() +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower, color = dist_label),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower, color = dist_label),
            linewidth = 0.9) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Flowering Probability') +
  theme(legend.position = 'none')

# Plot 2: Binned points + error bars + prediction lines
fig_fl_bin_combined <- ggplot() +
  geom_point(data = df_fl_binned, aes(x = logsize_t0, y = flower, color = dist_label)) +
  geom_errorbar(data = df_fl_binned, aes(x = logsize_t0, ymin = lwr, ymax = upr, color = dist_label),
                width = 0.2) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower, color = dist_label),
            linewidth = 0.9) +
  scale_color_manual(values = c('No disturbance' = 'black', 'disturbance' = 'red')) +
  theme_bw() +
  ylim(0, 1) +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Flowering Probability') +
  theme(legend.title = element_blank(), legend.position = 'top')

# Combine the two plots with patchwork
fig_fl_all <- fig_fl_line_combined + fig_fl_bin_combined +
  plot_annotation(
    title = 'Flowering',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = 'bold'),
      plot.subtitle = element_text(size = 10, face = 'italic')))

fig_fl_all
# ggsave(file.path(dir_result, 'flowering_by_size.png'), 
#        plot = fig_fl_all, width = 10, height = 5, dpi = 300)


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
    dist_label = ifelse(disturbance == 1, "Disturbance", "No disturbance")
  )

# Variance > Mean -> overdispersed -> negative binomial
mean(df_fr$fruit)
var(df_fr$fruit)


# Fruit model ------------------------------------------------------------------
# Fruit count conditional on flowering.
# Disturbance is included in all candidate models because disturbance is part
# of the biological model we want.

fr_mod_0 <- glm.nb(
  fruit ~ disturbance,
  data = df_fr
)

fr_mod_1 <- glm.nb(
  fruit ~ logsize_t0 + disturbance,
  data = df_fr
)

fr_mod_2 <- glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_fr
)

fr_mod_3 <- glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
  data = df_fr
)

fr_mods <- list(
  fr_mod_0,
  fr_mod_1,
  fr_mod_2,
  fr_mod_3
)

fr_mods_dAICc <- AICctab(
  fr_mods,
  weights = TRUE,
  sort = FALSE
)$dAICc

fr_mods_i_sort <- order(fr_mods_dAICc)

if (length(v_mod_set_fr) == 0) {
  fr_mod_i_best <- fr_mods_i_sort[1]
  v_mod_fr_index <- fr_mod_i_best - 1
} else {
  fr_mod_i_best <- v_mod_set_fr + 1
  v_mod_fr_index <- v_mod_set_fr
}

fr_mod_best  <- fr_mods[[fr_mod_i_best]]
fr_mod_ranef <- coef(fr_mod_best)

fr_mod_best
summary(fr_mod_best)
fr_mods_dAICc


# Prediction data for fruit model ---------------------------------------------

df_fr_pred <- expand.grid(
  logsize_t0 = seq(
    min(df_fr$logsize_t0, na.rm = TRUE),
    max(df_fr$logsize_t0, na.rm = TRUE),
    length.out = 200
  ),
  disturbance = c(0, 1)
) %>%
  as_tibble() %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    dist_label = ifelse(disturbance == 1, "Disturbance", "No disturbance"),
    pred = predict(fr_mod_best, newdata = ., type = "response")
  )


# Plot conditional fruit production -------------------------------------------

fig_fr <- ggplot(df_fr, aes(x = logsize_t0, y = fruit, colour = dist_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_line(
    data = df_fr_pred,
    aes(x = logsize_t0, y = pred, colour = dist_label),
    linewidth = 1
  ) +
  scale_colour_manual(
    values = c("No disturbance" = "black", "Disturbance" = "red")
  ) +
  theme_bw() +
  labs(
    title = "Fruit production conditional on flowering",
    subtitle = v_ggp_suffix,
    x = "log(Size at t0)",
    y = "Fruits per flowering plant",
    colour = NULL
  )

fig_fr

# ggsave(file.path(dir_result, 'fruits_conditional_on_flowering.png'),
#        plot = fig_fr, width = 8, height = 5, dpi = 300)

# Fruit to recruit: site-level transition --------------------------------------
# This model links total fruits at site in year t
# to total recruits at the same site in year t + 1.
# Fire/disturbance is intentionally excluded.
df_fr2re_site <- df %>%
  mutate(
    site = as.factor(site),
    year = as.numeric(year)
  ) %>%
  group_by(site, year) %>%
  summarise(
    fruits_site_t0 = sum(fruit, na.rm = TRUE),
    
    # recruit is coded NA / 1, so count actual recruit records
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
    logfruits_site_t0 = log1p(fruits_site_t0),
    logrecruits_site_t1 = log1p(recruits_site_t1)
  )


# Quick diagnostic -------------------------------------------------------------

df_fr2re_site %>%
  summarise(
    n_site_year_transitions = n(),
    n_positive_recruit_transitions = sum(recruits_site_t1 > 0),
    n_zero_recruit_transitions = sum(recruits_site_t1 == 0),
    total_fruits_t0 = sum(fruits_site_t0),
    total_recruits_t1 = sum(recruits_site_t1),
    max_fruits_site_t0 = max(fruits_site_t0),
    max_recruits_site_t1 = max(recruits_site_t1),
    proportion_zero_recruit = mean(recruits_site_t1 == 0)
  )


# Fruit-to-recruit models ------------------------------------------------------
# These models are forced through 0 on the log1p scale:
# log1p(recruits) = b1 * log1p(fruits).
#
# This means:
# fruits = 0 -> log1p(fruits) = 0
# predicted log1p(recruits) = 0
# predicted recruits = expm1(0) = 0

# 1. Fixed-effect slope-only model
fr2re_site_mod_1 <- glmmTMB(
  logrecruits_site_t1 ~ 0 + logfruits_site_t0,
  family = gaussian,
  REML = FALSE,
  data = df_fr2re_site)

# 2. Site-specific random slope only
fr2re_site_mod_2 <- glmmTMB(
  logrecruits_site_t1 ~ 0 + logfruits_site_t0 +
    (0 + logfruits_site_t0 | site),
  family = gaussian,
  REML = FALSE,
  data = df_fr2re_site)

fr2re_site_mods <- list(
  fr2re_site_mod_1,
  fr2re_site_mod_2)

fr2re_site_mods_dAICc <- AICctab(
  fr2re_site_mods,
  weights = TRUE,
  sort = FALSE)$dAICc

fr2re_site_mods_i_sort <- order(fr2re_site_mods_dAICc)

if (length(v_mod_set_fr2re_site) == 0) {
  fr2re_site_mod_i_best <- fr2re_site_mods_i_sort[1]
  v_mod_fr2re_site_index <- fr2re_site_mod_i_best
} else {
  fr2re_site_mod_i_best <- v_mod_set_fr2re_site
  v_mod_fr2re_site_index <- v_mod_set_fr2re_site
}

fr2re_site_mod_best <- fr2re_site_mods[[fr2re_site_mod_i_best]]


# Model convergence check ------------------------------------------------------

lapply(fr2re_site_mods, function(m) {
  list(
    convergence = m$fit$convergence,
    message = m$fit$message
  )
})


# Extract fruit-to-recruit fixed and site-specific effects ----------------------

fr2re_site_fixef <- fixef(fr2re_site_mod_best)$cond

fr2re_site_ranef <- tryCatch(
  {
    ranef(fr2re_site_mod_best)$cond$site %>%
      as.data.frame() %>%
      tibble::rownames_to_column("site")
  },
  error = function(e) {
    tibble(site = levels(df_fr2re_site$site))
  })

# Add missing random slope column as zero.
# There should be no intercept column because the model is forced through 0.
if (!"logfruits_site_t0" %in% names(fr2re_site_ranef)) {
  fr2re_site_ranef[["logfruits_site_t0"]] <- 0
}

df_fr2re_site_coef <- tibble(site = levels(df_fr2re_site$site)) %>%
  left_join(fr2re_site_ranef, by = "site") %>%
  mutate(
    logfruits_site_t0 = replace_na(logfruits_site_t0, 0),
    
    # Keep this only for compatibility with the rest of the IPM.
    # It is forced to zero.
    fr2re_b0_site = 0,
    
    fr2re_b1_site =
      unname(fr2re_site_fixef["logfruits_site_t0"]) + logfruits_site_t0) %>%
  select(site, fr2re_b0_site, fr2re_b1_site)

df_fr2re_site_coef


# Backtransform predicted log recruits to recruit counts -----------------------

predict_site_recruits <- function(logfruits_site_t0, fr2re_b0, fr2re_b1) {
  
  pred_log_recruits <- fr2re_b1 * logfruits_site_t0
  
  # Model was fit to log1p(recruits), so backtransform with expm1()
  pred_recruits <- expm1(pred_log_recruits)
  
  # Avoid impossible negative predictions
  pmax(pred_recruits, 0)
}


# Diagnostic prediction plot ---------------------------------------------------

df_fr2re_site_plot <- df_fr2re_site %>%
  left_join(df_fr2re_site_coef, by = "site") %>%
  mutate(
    pred_logrecruits_site_t1 = fr2re_b1_site * logfruits_site_t0
  ) %>%
  arrange(site, logfruits_site_t0)

fig_fr2re_site <- ggplot(
  df_fr2re_site_plot,
  aes(x = logfruits_site_t0, y = logrecruits_site_t1)
) +
  geom_point(
    alpha = 0.45,
    size = 1
  ) +
  geom_line(
    aes(y = pred_logrecruits_site_t1),
    linewidth = 0.9
  ) +
  facet_wrap(~ site, scales = "free_y") +
  labs(
    x = "log1p(total fruits at site in year t)",
    y = "log1p(total recruits at site in year t + 1)",
    title = "Site-level fruit-to-recruit transition",
    subtitle = "log1p(recruits) ~ 0 + log1p(fruits), forced through zero"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

fig_fr2re_site


# Site-specific coefficient table ----------------------------------------------

pars_by_site <- df_fr2re_site_coef %>%
  mutate(site = as.factor(site))

pars_by_site


# Extracting parameter estimates -----------------------------------------------
# Survival
coef_su_fe  <- data.frame(coefficient = names(coef(mod_su_bestfit)),
                          value       =       coef(mod_su_bestfit))

coef_su <- Reduce(function(...) rbind(...), list(coef_su_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Growth
coef_gr_fe  <- data.frame(coefficient = names(coef(mod_gr_bestfit)),
                          value       =       coef(mod_gr_bestfit))
coef_gr_var <- data.frame(coefficient = names(coef(mod_gr_var)),
                          value       =       coef(mod_gr_var))

coef_gr <- Reduce(function(...) rbind(...), list(coef_gr_fe, coef_gr_var)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Flower
coef_fl_fe  <- data.frame(coefficient = names(coef(mod_fl_bestfit)),
                          value       =       coef(mod_fl_bestfit))

coef_fl <- Reduce(function(...) rbind(...), list(coef_fl_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Fruit
coef_fr_fe  <- data.frame(coefficient = names(coef(fr_mod_best)),
                          value       =       coef(fr_mod_best))

coef_fr <- Reduce(function(...) rbind(...), list(coef_fr_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))


# Recruitment
df_re_size <- df %>% subset(recruit == 1)

# Miscellany
coef_misc   <- data.frame(
  coefficient = c('rec_siz', 'rec_sd', 'max_siz', 'min_siz'),
  value       = c(mean(log(df_re_size$size_t0), na.rm = T),
                  sd(  log(df_re_size$size_t0), na.rm = T),
                  df_gr$logsize_t0 %>% max,
                  df_gr$logsize_t0 %>% min))

extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = v_script_prefix,
  species = v_species,
  surv_b0 = extr_value(coef_su, 'b0'),
  surv_b1 = extr_value(coef_su, 'logsize_t0'),
  surv_b2 = extr_value(coef_su, 'logsize_t0_2'),
  surv_b3 = extr_value(coef_su, 'logsize_t0_3'),
  surv_bf = extr_value(coef_su, 'disturbance'),
  grow_b0 = extr_value(coef_gr, 'b0'),
  grow_b1 = extr_value(coef_gr, 'logsize_t0'),
  grow_b2 = extr_value(coef_gr, 'logsize_t0_2'),
  grow_b3 = extr_value(coef_gr, 'logsize_t0_3'),
  grow_bf = extr_value(coef_gr, 'disturbance'),
  a       = extr_value(coef_gr, 'a'),
  b       = extr_value(coef_gr, 'b'),
  fl_b0   = extr_value(coef_fl, 'b0'),
  fl_b1   = extr_value(coef_fl, 'logsize_t0'),
  fl_b2   = extr_value(coef_fl, 'logsize_t0_2'),
  fl_b3   = extr_value(coef_fl, 'logsize_t0_3'),
  fl_bf   = extr_value(coef_fl, 'disturbance'),
  fr_b0   = extr_value(coef_fr, 'b0'),
  fr_b1   = extr_value(coef_fr, 'logsize_t0'),
  fr_b2   = extr_value(coef_fr, 'logsize_t0_2'),
  fr_b3   = extr_value(coef_fr, 'logsize_t0_3'),
  fr_bf   = extr_value(coef_fr, 'disturbance'),
  fr2re_b0 = 0,
  fr2re_b1 = unname(fr2re_site_fixef["logfruits_site_t0"]),
  recr_sz = extr_value(coef_misc, 'rec_siz'),
  recr_sd = extr_value(coef_misc, 'rec_sd'),
  L       = extr_value(coef_misc, 'min_siz'),
  U       = extr_value(coef_misc, 'max_siz'),
  mat_siz = 200,
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index,
  mod_fl_index = v_mod_fl_index,
  mod_fr_index = v_mod_fr_index,
  mod_fr2re_site_index = v_mod_fr2re_site_index))


# Building the IPM -------------------------------------------------------------
# Function describing the invert logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

# Survival of x-sized individual to time t1
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  for (i in seq_len(num_pars)) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^(i)
    }
  }
  return(inv_logit(survival_value))
}

# Function describing standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- 0
  for (i in 0:num_pars) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}

# Function describing the transition kernel
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Flowering of x-sized individual at time t0
fl_x <- function(x, pars, num_pars = v_mod_fl_index) {
  val <- pars$fl_b0
  for (i in seq_len(num_pars)) {
    param <- paste0('fl_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  inv_logit(val)
}

# Fruiting of x-sized individuals at time t0
fr_x <- function(x, pars, num_pars = pars$mod_fr_index) {
  val <- pars$fr_b0
  for (i in seq_len(num_pars)) {
    param <- paste0('fr_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  exp(val)
}

# Expected fruits produced by an individual of size x
fruit_x <- function(x, pars) {
  fl_x(x, pars) * fr_x(x, pars)
}

# Make average annual observed size distribution for one site ------------------

make_initial_n_site <- function(pars, site_i, df_init = df) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  breaks <- seq(L, U, length.out = n + 1)
  
  years_i <- df_init %>%
    filter(site == site_i) %>%
    distinct(year) %>%
    arrange(year) %>%
    pull(year)
  
  if (length(years_i) == 0) {
    return(rep(0, n))
  }
  
  counts_by_year <- purrr::map(
    years_i,
    function(year_i) {
      
      x <- df_init %>%
        filter(
          site == site_i,
          year == year_i,
          !is.na(logsize_t0),
          is.finite(logsize_t0)
        ) %>%
        pull(logsize_t0)
      
      if (length(x) == 0) {
        return(rep(0, n))
      }
      
      hist(
        pmin(pmax(x, L), U),
        breaks = breaks,
        plot = FALSE,
        include.lowest = TRUE
      )$counts
    }
  )
  
  adult_counts_mean <- Reduce(`+`, counts_by_year) / length(counts_by_year)
  
  adult_density <- adult_counts_mean / h
  
  adult_density
}


# Expected total fruits at one site from its observed size distribution
expected_site_fruits <- function(pars, site_i, df_init = df) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  n_site <- make_initial_n_site(
    pars = pars, site_i = site_i, df_init = df_init)
  
  # n_site is a density, so n_site * h gives counts per size bin
  sum(fruit_x(y, pars) * n_site * h, na.rm = TRUE)
}


# Backtransform the site-level fruit-to-recruit model
predict_site_recruits_from_fruits <- function(site_fruits, pars) {
  
  if (is.na(site_fruits) || site_fruits <= 0) {
    return(0)
  }
  
  pred_log_recruits <- pars$fr2re_b1 * log1p(site_fruits)
  
  # Model was fit to log1p(recruits), so backtransform with expm1()
  pred_recruits <- expm1(pred_log_recruits)
  
  pmax(pred_recruits, 0)
}


# Convert site-level recruits back to recruits per fruit
site_recruits_per_fruit <- function(pars) {
  site_fruits <- pars$site_fruits_ref
  if (is.null(site_fruits) || is.na(site_fruits) || site_fruits <= 0) {
    return(0)
  }
  site_recruits <- predict_site_recruits_from_fruits(
    site_fruits = site_fruits,
    pars = pars)
  site_recruits / site_fruits
}

# Final site-specific parameter objects ----------------------------------------
pars_by_site <- pars_by_site %>%
  mutate(
    site_fruits_ref = purrr::map_dbl(
      site,
      ~ expected_site_fruits(
        pars = pars, site_i = .x, df_init = df)),
    
    site_recruits_ref = purrr::pmap_dbl(
      list(fr2re_b0_site, fr2re_b1_site, site_fruits_ref),
      function(fr2re_b0_site, fr2re_b1_site, site_fruits_ref) {
        p <- pars
        p$fr2re_b0 <- fr2re_b0_site
        p$fr2re_b1 <- fr2re_b1_site
        p$site_fruits_ref <- site_fruits_ref
        
        predict_site_recruits_from_fruits(
          site_fruits = site_fruits_ref,
          pars = p)
      }),
    
    recruits_per_fruit_ref = if_else(
      site_fruits_ref > 0, site_recruits_ref / site_fruits_ref, 0),
    
    pars_site = purrr::pmap(
      list(site, fr2re_b0_site, fr2re_b1_site, site_fruits_ref),
      function(site, fr2re_b0_site, fr2re_b1_site, site_fruits_ref) {
        p <- pars
        p$site <- site
        p$fr2re_b0 <- fr2re_b0_site
        p$fr2re_b1 <- fr2re_b1_site
        p$site_fruits_ref <- site_fruits_ref
        p
      }))

pars_by_site %>%
  select(
    site,
    fr2re_b0_site,
    fr2re_b1_site,
    site_fruits_ref,
    site_recruits_ref,
    recruits_per_fruit_ref)

# Recruitment size distribution at time t1 -------------------------------------

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


# F-kernel ---------------------------------------------------------------------

fyx <- function(y, x, pars) {
  
  fruit_x(x, pars) *
    site_recruits_per_fruit(pars) *
    re_y_dist(y, pars)
}


# Kernel -----------------------------------------------------------------------

kernel <- function(pars) {
  
  # number of bins over which to integrate
  n <- pars$mat_siz
  
  # lower and upper limits of integration
  L <- pars$L
  U <- pars$U
  
  # bin size
  h <- (U - L) / n
  
  # bin boundaries and mesh points
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  # Survival vector
  Smat <- sx(y, pars)
  
  # Growth matrix
  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Growth/survival transition matrix
  Tmat <- matrix(0, n, n)
  
  # Correct eviction of small offspring / small plants
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
  
  # Full kernel
  k_yx <- Fmat + Tmat
  
  return(list(
    k_yx    = k_yx,
    Fmat    = Fmat,
    Tmat    = Tmat,
    Gmat    = Gmat,
    meshpts = y
  ))
}


# Lambda function --------------------------------------------------------------

lambda_ipm <- function(pars) {
  Re(eigen(kernel(pars)$k_yx)$values[1])
}


# Site-specific IPM with disturbance regime -----------------------------------
# Goal:
# one lambda per site,
# with site-specific disturbance frequency
# and site-specific fruit-to-recruit transition.


# Site-specific disturbance frequency -----------------------------------------

p_disturbance_site <- df %>%
  group_by(site, year) %>%
  summarise(
    disturbance = max(disturbance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  summarise(
    p_disturbance = mean(disturbance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = factor(site, levels = levels(df$site)))

p_disturbance_site


# Apply one site's disturbance regime to the mean IPM parameters ---------------

apply_site_disturbance <- function(p, p_disturbance) {
  
  p2 <- p
  
  if (!is.null(p$surv_bf)) {
    p2$surv_b0 <- p$surv_b0 + p$surv_bf * p_disturbance
  }
  
  if (!is.null(p$grow_bf)) {
    p2$grow_b0 <- p$grow_b0 + p$grow_bf * p_disturbance
  }
  
  if (!is.null(p$fl_bf)) {
    p2$fl_b0 <- p$fl_b0 + p$fl_bf * p_disturbance
  }
  
  if (!is.null(p$fr_bf)) {
    p2$fr_b0 <- p$fr_b0 + p$fr_bf * p_disturbance
  }
  
  p2
}


# Build one parameter object per site -----------------------------------------

pars_by_site <- df_fr2re_site_coef %>%
  mutate(site = factor(site, levels = levels(df$site))) %>%
  left_join(p_disturbance_site, by = "site") %>%
  mutate(
    p_disturbance = replace_na(p_disturbance, 0),
    
    site_fruits_ref = purrr::map2_dbl(
      site,
      p_disturbance,
      function(site_i, p_dist_i) {
        
        p <- apply_site_disturbance(
          p = pars,
          p_disturbance = p_dist_i
        )
        
        expected_site_fruits(
          pars = p,
          site_i = site_i,
          df_init = df
        )
      }
    ),
    
    site_recruits_ref = purrr::pmap_dbl(
      list(fr2re_b0_site, fr2re_b1_site, site_fruits_ref, p_disturbance),
      function(fr2re_b0_site, fr2re_b1_site, site_fruits_ref, p_disturbance) {
        
        p <- apply_site_disturbance(
          p = pars,
          p_disturbance = p_disturbance
        )
        
        p$fr2re_b0 <- fr2re_b0_site
        p$fr2re_b1 <- fr2re_b1_site
        p$site_fruits_ref <- site_fruits_ref
        
        predict_site_recruits_from_fruits(
          site_fruits = site_fruits_ref,
          pars = p
        )
      }
    ),
    
    recruits_per_fruit_ref = if_else(
      site_fruits_ref > 0,
      site_recruits_ref / site_fruits_ref,
      0
    ),
    
    pars_site = purrr::pmap(
      list(site, fr2re_b0_site, fr2re_b1_site, site_fruits_ref, p_disturbance),
      function(site, fr2re_b0_site, fr2re_b1_site, site_fruits_ref, p_disturbance) {
        
        p <- apply_site_disturbance(
          p = pars,
          p_disturbance = p_disturbance
        )
        
        p$site <- site
        p$p_disturbance <- p_disturbance
        p$fr2re_b0 <- fr2re_b0_site
        p$fr2re_b1 <- fr2re_b1_site
        p$site_fruits_ref <- site_fruits_ref
        
        p
      }
    )
  )

pars_by_site %>%
  select(
    site,
    p_disturbance,
    fr2re_b0_site,
    fr2re_b1_site,
    site_fruits_ref,
    site_recruits_ref,
    recruits_per_fruit_ref
  ) %>%
  print(n = 100)


# Project one site-specific IPM ------------------------------------------------

project_site_ipm <- function(p, df_init = df) {
  
  n_obs <- make_initial_n_site(
    pars = p,
    site_i = p$site,
    df_init = df_init
  )
  
  K <- kernel(p)$k_yx
  
  h <- (p$U - p$L) / p$mat_siz
  
  n_proj <- K %*% n_obs
  
  tibble(
    lambda_asymptotic = Re(eigen(K)$values[1]),
    lambda_projected = as.numeric((sum(n_proj) * h) / (sum(n_obs) * h)),
    n_initial = sum(n_obs) * h,
    n_projected = sum(n_proj) * h,
    p_disturbance = p$p_disturbance,
    site_fruits_ref = p$site_fruits_ref,
    recruits_per_fruit_ref = site_recruits_per_fruit(p)
  )
}


# Modeled lambdas by site ------------------------------------------------------

df_lambda_model <- pars_by_site %>%
  transmute(
    site,
    lambda_data = purrr::map(pars_site, project_site_ipm)
  ) %>%
  tidyr::unnest(lambda_data)

df_lambda_model %>%
  print(n = 100)


# Observed site-level population growth ---------------------------------------

df_obs_counts_site_year <- df %>%
  filter(
    !is.na(logsize_t0),
    is.finite(logsize_t0)
  ) %>%
  group_by(site, year) %>%
  summarise(
    n_adults = n(),
    .groups = "drop"
  ) %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(
    n_adults_t1 = lead(n_adults),
    year_t1 = lead(year),
    year_gap = year_t1 - year,
    lambda_obs_year = n_adults_t1 / n_adults
  ) %>%
  ungroup() %>%
  filter(
    year_gap == 1,
    n_adults > 0,
    n_adults_t1 > 0
  )

df_lambda_obs <- df_obs_counts_site_year %>%
  group_by(site) %>%
  summarise(
    lambda_obs_arithmetic = mean(lambda_obs_year, na.rm = TRUE),
    lambda_obs_geometric = exp(mean(log(lambda_obs_year), na.rm = TRUE)),
    n_year_transitions = n(),
    mean_n_adults = mean(n_adults, na.rm = TRUE),
    .groups = "drop"
  )

df_lambda_obs %>%
  print(n = 100)


# Compare modeled and observed lambdas ----------------------------------------

df_lambda_compare <- df_lambda_model %>%
  left_join(df_lambda_obs, by = "site") %>%
  mutate(
    error_asymptotic_vs_obs = lambda_asymptotic - lambda_obs_geometric,
    error_projected_vs_obs = lambda_projected - lambda_obs_geometric
  )

df_lambda_compare %>%
  arrange(site) %>%
  print(n = 100)


# Plot observed versus modeled lambda -----------------------------------------

df_lambda_compare_plot <- df_lambda_compare %>%
  select(
    site,
    lambda_obs_geometric,
    lambda_asymptotic,
    lambda_projected
  ) %>%
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

fig_lambda_compare <- ggplot(
  df_lambda_compare_plot,
  aes(x = modeled_lambda, y = lambda_obs_geometric)
) +
  geom_point(size = 2.5) +
  geom_text(aes(label = site), nudge_y = 0.015, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~ modeled_lambda_type) +
  theme_bw() +
  labs(
    title = "Observed versus modeled population growth",
    subtitle = "Mean IPM with site-specific disturbance and fruit-to-recruit transition",
    x = expression("Modeled " * lambda),
    y = expression("Observed geometric " * lambda)
  )

fig_lambda_compare