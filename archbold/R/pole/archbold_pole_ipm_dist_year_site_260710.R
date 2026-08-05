# IPM year-specific with fire - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.10

# Study organism: Polygala lewtonii
# Citing publication: https://doi.org/10.1071/BT11271
# Time period: 2001-2024, through 2025


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
set.seed(100)
options(stringsAsFactors = F)


# Packages --------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(
  MASS,
  tidyverse,
  patchwork,
  skimr,
  ipmr,
  binom,
  bbmle,
  janitor,
  lme4,
  GGally)


# Specification ----------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Polygala lewtonii')
custom_delimiter <- c()
v_years_re <- c()

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Keep the manual model choices from the mean script unless changed here.
v_mod_set_su   <- c(6)
v_mod_set_gr   <- c(9)
v_mod_set_fl   <- c()
v_mod_set_fl_n <- c(2)
v_mod_set_re   <- c()

# Use site-specific recruitment random effects in projected lambdas? -----------
use_site_specific_recruitment <- TRUE


# Directory --------------------------------------------------------------------
dir_pub <- file.path(paste0(v_head))
dir_R <- file.path(dir_pub, 'R', v_sp_abb)
dir_data <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

if (!dir.exists(dir_R)) {dir.create(dir_R)}
if (!dir.exists(dir_data)) {dir.create(dir_data)}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Functions --------------------------------------------------------------------
source('helper_functions/plot_binned_prop.R')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df_og <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df_original_260519.csv')))

df_site <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df_sitehist_260519.csv')))

df_meta <- data.frame(
  variable = c(
    'unique_ID', 'quad', 'patch', 'year', 'month', 'date', 'qsurv',
    'stg', 'ht', 'mcd', 'st', 'flst'),
  description = c(
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
    'Number of flowering stems; Numeric'),
  stringsAsFactors = FALSE)


# Load workdata ----

df <- read.csv(
  file.path(
    dir_data,
    paste0('ab_', v_sp_abb, '_df_workdata_260519.csv')))



# Add disturbance from the exact previous year ----

df_dist_prev <- df %>%
  distinct(site, year, dist_transition) %>%
  transmute(
    site,
    year = year + 1,
    disturbance_prev = dist_transition)

df <- df %>%
  select(-any_of("disturbance_prev")) %>%
  left_join(df_dist_prev, by = c("site", "year")) %>%
  mutate(
    disturbance_prev = replace_na(disturbance_prev, 0),
    disturbance = factor(dist_transition, levels = c(0, 1)),
    disturbance_prev = factor(disturbance_prev, levels = c(0, 1)),
    year = as.integer(year),
    site = factor(site))

if (length(v_years_re) > 0) {
  df <- df %>% filter(!year %in% v_years_re)
}


# Controls ---------------------------------------------------------------------
ctrl_glmer <- glmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

ctrl_lmer <- lmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))


# Survival data ----------------------------------------------------------------
df_su <- df %>%
  filter(
    !is.na(survives),
    !is.na(logvol_t0),
    !is.na(logvol_t0_2),
    !is.na(logvol_t0_3),
    !is.na(recruits),
    size_t0 != 0) %>%
  mutate(
    year = factor(year),
    recruits = factor(recruits),
    disturbance = factor(disturbance, levels = c(0, 1))) %>%
  dplyr::select(
    id, site, year, size_t0, survives, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3,
    disturbance, stage, recruits)


# Survival model ---------------------------------------------------------------
mod_su_0 <- glmer(
  survives ~ disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_1 <- glmer(
  survives ~ logvol_t0 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_2 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_3 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance +
    (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_01 <- glmer(
  survives ~ recruits + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_11 <- glmer(
  survives ~ logvol_t0 + recruits + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_21 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + recruits + disturbance +
    (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_31 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits +
    disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_12 <- glmer(
  survives ~ logvol_t0 * recruits + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_22 <- glmer(
  survives ~ logvol_t0 * recruits + logvol_t0_2:recruits +
    disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_32 <- glmer(
  survives ~ logvol_t0 * recruits + logvol_t0_2:recruits +
    logvol_t0_3:recruits + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mods_su <- list(
  mod_su_0, mod_su_1, mod_su_2, mod_su_3,
  mod_su_01, mod_su_11, mod_su_21, mod_su_31,
  mod_su_12, mod_su_22, mod_su_32)
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


# Survival plots by year -------------------------------------------------------
make_su_year_plot <- function(year_i) {
  df_i <- df_su %>% filter(year == year_i)
  x <- seq(min(df_i$logvol_t0, na.rm = TRUE),
           max(df_i$logvol_t0, na.rm = TRUE), length.out = 100)
  pred_i <- tidyr::expand_grid(
    logvol_t0 = x,
    recruits = levels(df_su$recruits),
    disturbance = levels(df_su$disturbance)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      recruits = factor(recruits, levels = levels(df_su$recruits)),
      disturbance = factor(disturbance, levels = levels(df_su$disturbance)),
      year = factor(year_i, levels = levels(df_su$year)))
  
  pred_i <- pred_i %>%
    mutate(
      survives = predict(
        mod_su_best,
        newdata = pred_i,
        type = "response",
        re.form = NULL))
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, recruits, disturbance) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logvol_t0), n > 0)
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, survives, color = recruits, shape = disturbance),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, survives, color = recruits, linetype = disturbance),
      linewidth = 0.7) +
    scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
    scale_linetype_manual(values = c('0' = 'solid', '1' = 'dashed')) +
    labs(
      title = year_i,
      x = expression('log(volume)'[t0]),
      y = expression('Survival probability'[t1])) +
    ylim(0, 1) +
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


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(
    !is.na(size_t1), !is.na(size_t0),
    !is.na(logsize_t0), !is.na(logsize_t1),
    !is.na(logsize_t0_2), !is.na(logsize_t0_3),
    !is.na(volume_t0), !is.na(volume_t1),
    !is.na(logvol_t0), !is.na(logvol_t1),
    !is.na(logvol_t0_2), !is.na(logvol_t0_3),
    !is.na(recruits),
    size_t0 != 0,
    size_t1 != 0,
    is.finite(logvol_t0),
    is.finite(logvol_t1)) %>%
  mutate(
    year = factor(year),
    recruits = factor(recruits),
    disturbance = factor(disturbance, levels = c(0, 1))) %>%
  dplyr::select(
    id, site, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3,
    disturbance, stage, recruits)

fig_gr_raw <- ggplot(df_gr, aes(logvol_t0, logvol_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  labs(
    title = 'Growth',
    subtitle = v_ggp_suffix,
    x = expression('log(volume)'[t0]),
    y = expression('log(volume)'[t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_raw


# Growth model -----------------------------------------------------------------
mod_gr_0 <- lmer(
  logvol_t1 ~ disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_1 <- lmer(
  logvol_t1 ~ logvol_t0 + disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_2 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_3 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_01 <- lmer(
  logvol_t1 ~ recruits + disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_11 <- lmer(
  logvol_t1 ~ logvol_t0 + recruits + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_21 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + recruits + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_31 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    recruits + disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_12 <- lmer(
  logvol_t1 ~ logvol_t0 * recruits + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_22 <- lmer(
  logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits +
    disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_221 <- lmer(
  logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2 + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_32 <- lmer(
  logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits +
    logvol_t0_3:recruits + disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mods_gr <- list(
  mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3,
  mod_gr_01, mod_gr_11, mod_gr_21, mod_gr_31,
  mod_gr_12, mod_gr_22, mod_gr_221, mod_gr_32)
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


# Growth plots by year ---------------------------------------------------------
make_gr_year_plot <- function(year_i) {
  df_i <- df_gr %>% filter(year == year_i)
  x <- seq(min(df_i$logvol_t0, na.rm = TRUE),
           max(df_i$logvol_t0, na.rm = TRUE), length.out = 100)
  pred_i <- tidyr::expand_grid(
    logvol_t0 = x,
    recruits = levels(df_gr$recruits),
    disturbance = levels(df_gr$disturbance)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      recruits = factor(recruits, levels = levels(df_gr$recruits)),
      disturbance = factor(disturbance, levels = levels(df_gr$disturbance)),
      year = factor(year_i, levels = levels(df_gr$year)))
  
  pred_i <- pred_i %>%
    mutate(
      logvol_t1 = predict(
        mod_gr_best,
        newdata = pred_i,
        re.form = NULL))
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, recruits, disturbance) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      logvol_t1 = mean(logvol_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logvol_t0), n > 0)
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, logvol_t1, color = recruits, shape = disturbance),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, logvol_t1, color = recruits, linetype = disturbance),
      linewidth = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', lty = 2) +
    scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
    scale_linetype_manual(values = c('0' = 'solid', '1' = 'dashed')) +
    labs(
      title = year_i,
      x = expression('log(volume)'[t0]),
      y = expression('log(volume)'[t1])) +
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


# Flower data ------------------------------------------------------------------
# Recruits should theoretically be excluded from flowering probability and flower-number models.
# That exclusion is not supported by the data, because 35 recruits did flower, 
# including several fairly large individuals and some stage-3 plants.
df_fl <- df %>%
  filter(
    !is.na(flower),
    !is.na(logvol_t0),
    !is.na(logvol_t0_2),
    !is.na(logvol_t0_3)) %>%
  mutate(
    year = factor(year),
    disturbance_prev = factor(
      disturbance_prev, levels = c(0, 1))) %>%
  dplyr::select(
    id, site, year, size_t0, flower, fl_nr, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3,
    disturbance_prev, stage, recruits)


# Flower model -----------------------------------------------------------------
mod_fl_00 <- glmer(
  flower ~ 1 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_01 <- glmer(
  flower ~ disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_10 <- glmer(
  flower ~ logvol_t0 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_11 <- glmer(
  flower ~ logvol_t0 + disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_20 <- glmer(
  flower ~ logvol_t0 + logvol_t0_2 +
    (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_21 <- glmer(
  flower ~ logvol_t0 + logvol_t0_2 + disturbance_prev +
    (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_30 <- glmer(
  flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)
mod_fl_31 <- glmer(
  flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mods_fl <- list(mod_fl_00, mod_fl_01, mod_fl_10, mod_fl_11, 
                mod_fl_20, mod_fl_21, mod_fl_30, mod_fl_31)
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


# Flowering plots by year ------------------------------------------------------
make_fl_year_plot <- function(year_i) {
  df_i <- df_fl %>% filter(year == year_i)
  x <- seq(min(df_i$logvol_t0, na.rm = TRUE),
           max(df_i$logvol_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- tidyr::expand_grid(
    logvol_t0 = x,
    disturbance_prev = levels(df_fl$disturbance_prev)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      disturbance_prev = factor(
        disturbance_prev, levels = levels(df_fl$disturbance_prev)),
      year = factor(year_i, levels = levels(df_fl$year)))
  
  pred_i <- pred_i %>%
    mutate(
      flower = predict(
        mod_fl_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, disturbance_prev) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      flower = mean(flower, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logvol_t0), n > 0)
  
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, flower, color = disturbance_prev),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, flower, color = disturbance_prev),
      linewidth = 0.7) +
    scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(volume)'[t0]),
      y = 'Flowering probability') +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fl_yrs <- lapply(levels(df_fl$year), make_fl_year_plot)
fig_fl_years <- wrap_plots(fl_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Flowering - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fl_years


# Flower number conditional on flowering data ----------------------------------
df_fl_cond <- df_fl %>%
  filter(flower == 1, !is.na(fl_nr), fl_nr > 0, fl_nr == round(fl_nr)) %>%
  mutate(year = factor(year))


# Flower number model ----------------------------------------------------------
mod_fl_n_0 <- glmer.nb(
  fl_nr ~ disturbance_prev + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)
mod_fl_n_1 <- glmer.nb(
  fl_nr ~ logvol_t0 + disturbance_prev + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)
mod_fl_n_2 <- glmer.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2 + disturbance_prev +
    (1 | year),
  data = df_fl_cond, control = ctrl_glmer)
mod_fl_n_3 <- glmer.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance_prev +
    (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mods_fl_n <- list(mod_fl_n_0, mod_fl_n_1, mod_fl_n_2, mod_fl_n_3) # 
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

mod_fl_n_best
summary(mod_fl_n_best)
mods_fl_n_dAIC


# Flower number plots by year --------------------------------------------------
make_fl_n_year_plot <- function(year_i) {
  df_i <- df_fl_cond %>% filter(year == year_i)
  x <- seq(min(df_i$logvol_t0, na.rm = TRUE),
           max(df_i$logvol_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- tidyr::expand_grid(
    logvol_t0 = x,
    disturbance_prev = levels(df_fl_cond$disturbance_prev)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      disturbance_prev = factor(
        disturbance_prev, levels = levels(df_fl_cond$disturbance_prev)),
      year = factor(year_i, levels = levels(df_fl_cond$year)))
  
  pred_i <- pred_i %>%
    mutate(
      fl_nr = predict(
        mod_fl_n_best, newdata = pred_i, type = 'response',
        re.form = NULL))
  
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, disturbance_prev) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      fl_nr = mean(fl_nr, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logvol_t0), n > 0)
  
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, fl_nr, color = disturbance_prev),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, fl_nr, color = disturbance_prev),
      linewidth = 0.7) +
    scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
    labs(
      title = year_i,
      x = expression('log(volume)'[t0]),
      y = 'Number of flowering stems') +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')
}

fl_n_yrs <- lapply(levels(df_fl_cond$year), make_fl_n_year_plot)
fig_fl_n_years <- wrap_plots(fl_n_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Flower number - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_fl_n_years


# Flowering stems to recruits --------------------------------------------------
# Biological note:
# P. lewtonii has post-fire seedbank recruitment, so recruitment is modeled from
# site-year flowering stems at t0 to recruits at t1, with fire as a fixed effect.
df_fs2r <- df %>%
  group_by(site, year) %>%
  summarise(
    fs_t0 = sum(fl_nr, na.rm = TRUE),
    disturbance_t0 = factor(
      max(as.numeric(as.character(disturbance)), na.rm = TRUE),
      levels = c(0, 1)),
    .groups = "drop") %>%
  mutate(year_t1 = year + 1) %>%
  left_join(
    df %>%
      filter(recruits == 1) %>%
      count(site, year, name = "re_t1"),
    by = c("site", "year_t1" = "year")) %>%
  mutate(
    re_t1 = replace_na(re_t1, 0L),
    site = factor(site),
    year_t0 = factor(year),
    year_t1 = factor(year_t1),
    fs_t0_log = log1p(fs_t0),
    re_t1_log = log1p(re_t1))

# Keep the old mean-script recruitment plot as a quick check.
fig_fs2r_raw <- ggplot(
  df_fs2r,
  aes(x = fs_t0, y = re_t1, color = disturbance_t0)) +
  geom_jitter(height = 0.2, width = 0.5, alpha = 0.4) +
  scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
  theme_bw() +
  labs(
    title = 'Recruits t1 by flowering stems t0',
    subtitle = v_ggp_suffix,
    x = 'Total flowering stems at site in year t',
    y = 'Number of recruits at site in year t + 1',
    color = 'Fire')

fig_fs2r_raw


# Recruitment model ------------------------------------------------------------
ctrl_re <- lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))

mod_re_00 <- lm(
  re_t1_log ~ fs_t0_log + disturbance_t0,
  data = df_fs2r)

mod_re_10 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (1 | site),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_11 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (1 | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_12 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (1 | site) + (1 | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_20 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (fs_t0_log | site),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_21 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (fs_t0_log | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_22 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 +
    (fs_t0_log | site) + (fs_t0_log | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_30 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (0 + fs_t0_log | site),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_31 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 + (0 + fs_t0_log | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mod_re_32 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance_t0 +
    (0 + fs_t0_log | site) + (0 + fs_t0_log | year_t1),
  data = df_fs2r, REML = FALSE, control = ctrl_re)

mods_re <- list(
  mod_re_00,
  mod_re_10, mod_re_11, mod_re_12,
  mod_re_20, mod_re_21, mod_re_22,
  mod_re_30, mod_re_31, mod_re_32)
mods_re_dAIC <- bbmle::AICctab(mods_re, weights = TRUE, sort = FALSE)$dAIC
mods_re_sorted <- order(mods_re_dAIC)

if (length(v_mod_set_re) == 0) {
  mod_re_index_bestfit <- mods_re_sorted[1]
  v_mod_re_index <- mod_re_index_bestfit - 1
} else {
  mod_re_index_bestfit <- v_mod_set_re + 1
  v_mod_re_index <- v_mod_set_re
}

mod_re_best <- mods_re[[mod_re_index_bestfit]]

mod_re_best
summary(mod_re_best)
mods_re_dAIC


# Recruitment plots by recruit year -------------------------------------------
predict_re_safe <- function(model, newdata, re.form = NULL) {
  if (inherits(model, 'merMod')) {
    out <- predict(
      model, newdata = newdata, re.form = re.form, allow.new.levels = TRUE)
  } else {
    out <- predict(model, newdata = newdata)
  }
  return(out)
}

shape_values_re <- c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6, 9)
site_shapes_re <- setNames(
  shape_values_re[seq_along(levels(df_fs2r$site))],
  levels(df_fs2r$site))

make_re_year_plot <- function(year_i) {
  df_i <- df_fs2r %>% filter(year_t1 == year_i)
  x <- seq(min(df_i$fs_t0_log, na.rm = TRUE),
           max(df_i$fs_t0_log, na.rm = TRUE), length.out = 100)
  
  pred_fixed <- tidyr::expand_grid(
    fs_t0_log = x,
    disturbance_t0 = levels(df_fs2r$disturbance_t0)) %>%
    mutate(
      site = factor(
        levels(df_fs2r$site)[1], levels = levels(df_fs2r$site)),
      year_t1 = factor(year_i, levels = levels(df_fs2r$year_t1)),
      disturbance_t0 = factor(
        disturbance_t0, levels = levels(df_fs2r$disturbance_t0)))
  
  pred_fixed <- pred_fixed %>%
    mutate(
      pred_re = predict_re_safe(
        mod_re_best, pred_fixed, re.form = NA))
  
  pred_site <- tidyr::expand_grid(
    fs_t0_log = x,
    site = levels(df_fs2r$site),
    disturbance_t0 = levels(df_fs2r$disturbance_t0)) %>%
    mutate(
      site = factor(site, levels = levels(df_fs2r$site)),
      year_t1 = factor(year_i, levels = levels(df_fs2r$year_t1)),
      disturbance_t0 = factor(
        disturbance_t0, levels = levels(df_fs2r$disturbance_t0)))
  
  pred_site <- pred_site %>%
    mutate(
      pred_re = predict_re_safe(
        mod_re_best, pred_site, re.form = NULL))
  
  ggplot() +
    geom_point(
      data = df_i,
      aes(fs_t0_log, re_t1_log, color = site, shape = site),
      alpha = 0.75, size = 2) +
    geom_line(
      data = pred_site,
      aes(
        fs_t0_log, pred_re, color = site,
        group = interaction(site, disturbance_t0),
        linetype = disturbance_t0),
      linewidth = 0.55, alpha = 0.55) +
    geom_line(
      data = pred_fixed,
      aes(fs_t0_log, pred_re, linetype = disturbance_t0),
      inherit.aes = FALSE, color = 'black', linewidth = 1.1) +
    scale_color_viridis_d(option = 'turbo', end = 0.95) +
    scale_shape_manual(values = site_shapes_re) +
    scale_linetype_manual(values = c('0' = 'solid', '1' = 'dashed')) +
    theme_bw() +
    labs(
      title = year_i,
      x = expression('log(1 + flowering stems)'[t0]),
      y = expression('log(1 + recruits)'[t1]),
      color = 'Site',
      shape = 'Site') +
    theme(text = element_text(size = 5), legend.position = 'none')
}

re_yrs <- lapply(levels(df_fs2r$year_t1), make_re_year_plot)
fig_re_years <- wrap_plots(re_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Recruitment - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_re_years


# Recruit size distribution ----------------------------------------------------
df_recr <- df %>%
  filter(recruits == 1, volume_t0 > 0, is.finite(logvol_t0))

rc_sz <- data.frame(
  mean = mean(df_recr$logvol_t0, na.rm = TRUE),
  sd = sd(df_recr$logvol_t0, na.rm = TRUE))


# Parameter extraction helpers -------------------------------------------------
empty_coef_df <- function() {
  data.frame(coefficient = character(), value = numeric())
}

is_mixed_model <- function(model) {
  inherits(model, 'merMod')
}

get_fixef_safe <- function(model) {
  if (is_mixed_model(model)) {
    fixef(model)
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
    if (grepl('^regex:', alias)) {
      pattern <- sub('^regex:', '', alias)
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

check_duplicate_coefficients <- function(x, object_name = 'parameter table') {
  dup <- x$coefficient[duplicated(x$coefficient)]
  if (length(dup) > 0) {
    stop(
      object_name,
      ' contains duplicated coefficient names: ',
      paste(unique(dup), collapse = ', '))
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
    pivot_wider(names_from = coefficient, values_from = value) %>%
    as.list()
}

extract_fixed_pars <- function(model, term_map, fill_missing = TRUE) {
  model_coefs <- get_fixef_safe(model)
  model_terms <- names(model_coefs)
  out <- lapply(names(term_map), function(coef_name) {
    model_term <- find_model_term(model_terms, term_map[[coef_name]])
    if (is.na(model_term)) {
      if (fill_missing) {
        value <- 0
      } else {
        return(NULL)
      }
    } else {
      value <- unname(model_coefs[model_term])
    }
    data.frame(coefficient = coef_name, value = value)
  })
  bind_coef_rows(out)
}

extract_group_pars <- function(model, group_var, term_map, fill_missing = TRUE) {
  coef_matrix <- get_group_coef_safe(model, group_var)
  if (is.null(coef_matrix)) {
    return(empty_coef_df())
  }
  model_terms <- colnames(coef_matrix)
  out <- lapply(names(term_map), function(coef_prefix) {
    model_term <- find_model_term(model_terms, term_map[[coef_prefix]])
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
      coefficient = paste0(coef_prefix, '_', rownames(coef_matrix)),
      value = value)
  })
  bind_coef_rows(out)
}

extract_re_ranef_devs <- function(model, group_var, group_label) {
  if (!inherits(model, 'merMod')) {
    return(empty_coef_df())
  }
  ranef_list <- ranef(model)
  if (!(group_var %in% names(ranef_list))) {
    return(empty_coef_df())
  }
  ranef_matrix <- ranef_list[[group_var]]
  out <- list()
  if ('(Intercept)' %in% colnames(ranef_matrix)) {
    out[['intercept']] <- data.frame(
      coefficient = paste0('re_u0_', group_label, '_',
                           rownames(ranef_matrix)),
      value = ranef_matrix[, '(Intercept)'])
  }
  if ('fs_t0_log' %in% colnames(ranef_matrix)) {
    out[['slope']] <- data.frame(
      coefficient = paste0('re_u1_', group_label, '_',
                           rownames(ranef_matrix)),
      value = ranef_matrix[, 'fs_t0_log'])
  }
  bind_coef_rows(out)
}


# Term maps --------------------------------------------------------------------
size_term_map <- list(
  b0 = c('(Intercept)'),
  b1 = c('logvol_t0'),
  b2 = c('logvol_t0_2', 'I(logvol_t0^2)'),
  b3 = c('logvol_t0_3', 'I(logvol_t0^3)'),
  bf = c('disturbance1', 'regex:^disturbance'),
  br = c('recruits1', 'regex:^recruits1$'),
  b1_br = c('logvol_t0:recruits1', 'recruits1:logvol_t0'),
  b2_br0 = c('logvol_t0_2:recruits0', 'recruits0:logvol_t0_2'),
  b2_br1 = c('logvol_t0_2:recruits1', 'recruits1:logvol_t0_2'),
  b3_br0 = c('logvol_t0_3:recruits0', 'recruits0:logvol_t0_3'),
  b3_br1 = c('logvol_t0_3:recruits1', 'recruits1:logvol_t0_3'))

# Flowering disturbance term map
lagged_size_term_map <- list(
  b0 = c("(Intercept)"),
  b1 = c("logvol_t0"),
  b2 = c("logvol_t0_2", "I(logvol_t0^2)"),
  b3 = c("logvol_t0_3", "I(logvol_t0^3)"),
  bf = c(
    "disturbance_prev1",
    "regex:^disturbance_prev"))

make_term_map <- function(prefix, map) {
  out <- map
  names(out) <- paste0(prefix, names(out))
  out
}

su_term_map <- make_term_map('surv_', size_term_map)
gr_term_map <- make_term_map('grow_', size_term_map)
fl_term_map  <- make_term_map("fl_",  lagged_size_term_map)
fln_term_map <- make_term_map("fln_", lagged_size_term_map)
re_term_map <- list(
  re_b0 = c('(Intercept)'),
  re_b1 = c('fs_t0_log'),
  re_bf = c('disturbance_t01', 'regex:^disturbance_t0'))


# Fixed parameters -------------------------------------------------------------
su_fe <- extract_fixed_pars(mod_su_best, su_term_map, fill_missing = TRUE)
gr_fe <- extract_fixed_pars(mod_gr_best, gr_term_map, fill_missing = TRUE)
fl_fe <- extract_fixed_pars(mod_fl_best, fl_term_map, fill_missing = TRUE)
fln_fe <- extract_fixed_pars(mod_fl_n_best, fln_term_map, fill_missing = TRUE)
re_fe <- extract_fixed_pars(mod_re_best, re_term_map, fill_missing = TRUE)


# Constants --------------------------------------------------------------------
gr_var_coef <- coef(mod_gr_var)

constants <- tibble::tribble(
  ~coefficient, ~value,
  'recr_sz', rc_sz$mean,
  'recr_sd', rc_sz$sd,
  'a', as.numeric(gr_var_coef[1]),
  'b', as.numeric(gr_var_coef[2]),
  'L', min(c(df_gr$logvol_t0, df_fl$logvol_t0), na.rm = TRUE) - 0.1,
  'U', max(c(df_gr$logvol_t0, df_fl$logvol_t0), na.rm = TRUE) + 0.1,
  'mat_siz', 200,
  'fs_t0_ref', mean(df_fs2r$fs_t0, na.rm = TRUE),
  'fs_t0_ref_log', log1p(mean(df_fs2r$fs_t0, na.rm = TRUE)),
  're_sigma', sigma(mod_re_best),
  'mod_su_index', v_mod_su_index,
  'mod_gr_index', v_mod_gr_index,
  'mod_fl_index', v_mod_fl_index,
  'mod_fl_n_index', v_mod_fl_n_index,
  'mod_re_index', v_mod_re_index) %>%
  mutate(coefficient = as.character(coefficient), value = as.numeric(value))

pars_cons <- bind_coef_rows(list(su_fe, gr_fe, fl_fe, fln_fe, re_fe,
                                 constants))
check_duplicate_coefficients(pars_cons, object_name = 'pars_cons')
pars_all_mean <- coef_df_to_list(pars_cons)
pars_mean <- pars_all_mean


# Year-varying parameters ------------------------------------------------------
su_out_yr <- extract_group_pars(mod_su_best, 'year', su_term_map,
                                fill_missing = TRUE)
gr_out_yr <- extract_group_pars(mod_gr_best, 'year', gr_term_map,
                                fill_missing = TRUE)
fl_out_yr <- extract_group_pars(mod_fl_best, 'year', fl_term_map,
                                fill_missing = TRUE)
fln_out_yr <- extract_group_pars(mod_fl_n_best, 'year', fln_term_map,
                                 fill_missing = TRUE)

pars_var <- bind_coef_rows(list(su_out_yr, gr_out_yr, fl_out_yr,
                                fln_out_yr))
check_duplicate_coefficients(pars_var, object_name = 'pars_var')
pars_all_year <- coef_df_to_list(pars_var)


# Site and recruitment-year random-effect deviations ---------------------------
pars_re_site <- extract_re_ranef_devs(mod_re_best, 'site', 'site') %>%
  mutate(coefficient = as.character(coefficient))
check_duplicate_coefficients(pars_re_site, object_name = 'pars_re_site')
pars_all_site <- coef_df_to_list(pars_re_site)

pars_re_year_t1 <- extract_re_ranef_devs(
  mod_re_best, 'year_t1', 'year_t1') %>%
  mutate(coefficient = as.character(coefficient))
check_duplicate_coefficients(
  pars_re_year_t1, object_name = 'pars_re_year_t1')
pars_all_re_year_t1 <- coef_df_to_list(pars_re_year_t1)


# IPM helper functions ---------------------------------------------------------
inv_logit <- function(x) {
  plogis(x)
}

get_par <- function(pars, par_name, default = 0) {
  if (!is.null(pars[[par_name]])) {
    return(pars[[par_name]])
  }
  default
}

make_ipm_pars <- function(
    pars_mean,
    pars_year = NULL,
    pars_site = NULL,
    pars_re_year_t1 = NULL,
    year = NULL,
    site = NULL,
    re_year_t1 = NULL) {
  pars <- pars_mean
  if (!is.null(year) && !is.null(pars_year)) {
    year <- as.character(year)
    year_hits <- grep(paste0('_', year, '$'), names(pars_year), value = TRUE)
    for (nm in year_hits) {
      base_nm <- sub(paste0('_', year, '$'), '', nm)
      pars[[base_nm]] <- pars_year[[nm]]
    }
  }
  if (!is.null(site) && !is.null(pars_site)) {
    site <- as.character(site)
    pars$re_b0 <- get_par(pars, 're_b0', 0) +
      get_par(pars_site, paste0('re_u0_site_', site), 0)
    pars$re_b1 <- get_par(pars, 're_b1', 0) +
      get_par(pars_site, paste0('re_u1_site_', site), 0)
  }
  if (is.null(re_year_t1) && !is.null(year)) {
    re_year_t1 <- as.numeric(as.character(year)) + 1
  }
  if (!is.null(re_year_t1) && !is.null(pars_re_year_t1)) {
    re_year_t1 <- as.character(re_year_t1)
    pars$re_b0 <- get_par(pars, 're_b0', 0) +
      get_par(pars_re_year_t1, paste0('re_u0_year_t1_', re_year_t1), 0)
    pars$re_b1 <- get_par(pars, 're_b1', 0) +
      get_par(pars_re_year_t1, paste0('re_u1_year_t1_', re_year_t1), 0)
  }
  pars
}

surv_lp <- function(x, pars, stage = c('adult', 'recruit'), disturbance = 0) {
  stage <- match.arg(stage)
  eta <- get_par(pars, 'surv_b0', 0) +
    get_par(pars, 'surv_b1', 0) * x +
    get_par(pars, 'surv_b2', 0) * x^2 +
    get_par(pars, 'surv_b3', 0) * x^3 +
    get_par(pars, 'surv_bf', 0) * disturbance
  if (stage == 'recruit') {
    eta <- eta +
      get_par(pars, 'surv_br', 0) +
      get_par(pars, 'surv_b1_br', 0) * x +
      get_par(pars, 'surv_b2_br1', 0) * x^2 +
      get_par(pars, 'surv_b3_br1', 0) * x^3
  } else {
    eta <- eta +
      get_par(pars, 'surv_b2_br0', 0) * x^2 +
      get_par(pars, 'surv_b3_br0', 0) * x^3
  }
  eta
}

sx <- function(x, pars, stage = c('adult', 'recruit'), disturbance = 0) {
  inv_logit(surv_lp(x, pars, stage, disturbance))
}

grow_mu <- function(x, pars, stage = c('adult', 'recruit'), disturbance = 0) {
  stage <- match.arg(stage)
  mu <- get_par(pars, 'grow_b0', 0) +
    get_par(pars, 'grow_b1', 0) * x +
    get_par(pars, 'grow_b2', 0) * x^2 +
    get_par(pars, 'grow_b3', 0) * x^3 +
    get_par(pars, 'grow_bf', 0) * disturbance
  if (stage == 'recruit') {
    mu <- mu +
      get_par(pars, 'grow_br', 0) +
      get_par(pars, 'grow_b1_br', 0) * x +
      get_par(pars, 'grow_b2_br1', 0) * x^2 +
      get_par(pars, 'grow_b3_br1', 0) * x^3
  } else {
    mu <- mu +
      get_par(pars, 'grow_b2_br0', 0) * x^2 +
      get_par(pars, 'grow_b3_br0', 0) * x^3
  }
  mu
}

grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * x))
}

gxy <- function(x, y, pars, stage = c('adult', 'recruit'), disturbance = 0) {
  stage <- match.arg(stage)
  dnorm(
    y,
    mean = grow_mu(x, pars, stage = stage, disturbance = disturbance),
    sd = grow_sd(x, pars))
}

fl_x <- function(x, pars, disturbance_prev = 0) {
  eta <- get_par(pars, "fl_b0", 0) +
    get_par(pars, "fl_b1", 0) * x +
    get_par(pars, "fl_b2", 0) * x^2 +
    get_par(pars, "fl_b3", 0) * x^3 +
    get_par(pars, "fl_bf", 0) * disturbance_prev
  inv_logit(eta)
}

fl_n_x <- function(x, pars, disturbance_prev = 0) {
  eta <- get_par(pars, "fln_b0", 0) +
    get_par(pars, "fln_b1", 0) * x +
    get_par(pars, "fln_b2", 0) * x^2 +
    get_par(pars, "fln_b3", 0) * x^3 +
    get_par(pars, "fln_bf", 0) * disturbance_prev
  exp(eta)
}

fs_x <- function(x, pars, disturbance_prev = 0) {
  fl_x(x, pars, disturbance_prev) *
    fl_n_x(x, pars, disturbance_prev)
}

re_total_ref <- function(pars, disturbance = 0) {
  re_log <- get_par(pars, 're_b0', 0) +
    get_par(pars, 're_b1', 0) * get_par(pars, 'fs_t0_ref_log', 0) +
    get_par(pars, 're_bf', 0) * disturbance
  re_value <- expm1(re_log + 0.5 * get_par(pars, 're_sigma', 0)^2)
  pmax(re_value, 0)
}

rx <- function(
    x,
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  fs_value <- fs_x(
    x, pars, disturbance_prev = disturbance_prev)
  re_ref <- re_total_ref(
    pars, disturbance = disturbance)
  re_per_fs <- re_ref /
    pmax(get_par(pars, "fs_t0_ref", 0), .Machine$double.eps)
  fs_value * re_per_fs
}

re_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
  if (!is.null(h)) {
    dens <- dens / sum(dens * h)
  }
  dens
}

fyx <- function(
    y,
    x,
    pars,
    h,
    disturbance = 0,
    disturbance_prev = 0) {
  Pvec <- re_y_dist(y, pars, h = h)
  Rvec <- rx(
    x,
    pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  outer(Pvec, Rvec) * h
}


# Two-stage kernel -------------------------------------------------------------
kernel <- function(
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  S_rec <- sx(y, pars, stage = 'recruit', disturbance = disturbance)
  S_ad <- sx(y, pars, stage = 'adult', disturbance = disturbance)
  G_rec <- matrix(0, n, n)
  G_ad <- matrix(0, n, n)
  G_rec[] <- t(outer(
    y, y, Vectorize(function(x, y) {
      gxy(x, y, pars, stage = 'recruit', disturbance = disturbance)
    }))) * h
  G_ad[] <- t(outer(
    y, y, Vectorize(function(x, y) {
      gxy(x, y, pars, stage = 'adult', disturbance = disturbance)
    }))) * h
  T_rec <- matrix(0, n, n)
  T_ad <- matrix(0, n, n)
  for (i in 1:n) {
    if (i <= n / 2) {
      G_rec[1, i] <- G_rec[1, i] + 1 - sum(G_rec[, i])
      G_ad[1, i] <- G_ad[1, i] + 1 - sum(G_ad[, i])
    } else {
      G_rec[n, i] <- G_rec[n, i] + 1 - sum(G_rec[, i])
      G_ad[n, i] <- G_ad[n, i] + 1 - sum(G_ad[, i])
    }
    T_rec[, i] <- G_rec[, i] * S_rec[i]
    T_ad[, i] <- G_ad[, i] * S_ad[i]
  }
  
  F_ad <- fyx(
    y = y,
    x = y,
    pars = pars,
    h = h,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  K_rec_rec <- matrix(0, n, n)
  K_rec_ad <- F_ad
  K_ad_rec <- T_rec
  K_ad_ad <- T_ad
  k_yx <- rbind(cbind(K_rec_rec, K_rec_ad), cbind(K_ad_rec, K_ad_ad))
  list(
    k_yx = k_yx,
    K_rec_rec = K_rec_rec,
    K_rec_ad = K_rec_ad,
    K_ad_rec = K_ad_rec,
    K_ad_ad = K_ad_ad,
    F_ad = F_ad,
    T_rec = T_rec,
    T_ad = T_ad,
    G_rec = G_rec,
    G_ad = G_ad,
    S_rec = S_rec,
    S_ad = S_ad,
    meshpts = y,
    h = h,
    L = L,
    U = U)
}

lambda_ipm <- function(
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  Re(eigen(
    kernel(
      pars,
      disturbance = disturbance,
      disturbance_prev = disturbance_prev)$k_yx)$values[1])
}


# Mean IPMs --------------------------------------------------------------------

lambda_mean <- data.frame(
  scenario = c(
    'Undisturbed',
    'Current fire',
    'Previous fire',
    'Consecutive fire'),
  disturbance = c(0, 1, 0, 1),
  disturbance_prev = c(0, 0, 1, 1)) %>%
  mutate(
    lambda = purrr::map2_dbl(
      disturbance,
      disturbance_prev,
      ~ lambda_ipm(
        pars_mean,
        disturbance = .x,
        disturbance_prev = .y)))

lambda_mean


# Year-specific IPMs -----------------------------------------------------------

lambda_ipm_year <- function(
    year,
    disturbance = 0,
    disturbance_prev = 0,
    site = NULL,
    re_year_t1 = NULL) {
  pars_i <- make_ipm_pars(
    pars_mean = pars_all_mean,
    pars_year = pars_all_year,
    pars_site = pars_all_site,
    pars_re_year_t1 = pars_all_re_year_t1,
    year = year,
    site = site,
    re_year_t1 = re_year_t1)
  
  lambda_ipm(
    pars_i,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
}


# Year-specific disturbance scenarios -----------------------------------------

ipm_years <- sort(unique(as.integer(as.character(df_su$year))))

lambda_year <- tidyr::expand_grid(
  year = ipm_years,
  scenario = c(
    'Undisturbed',
    'Current fire',
    'Previous fire',
    'Consecutive fire')) %>%
  mutate(
    disturbance = case_when(
      scenario %in% c('Current fire', 'Consecutive fire') ~ 1,
      TRUE ~ 0),
    disturbance_prev = case_when(
      scenario %in% c('Previous fire', 'Consecutive fire') ~ 1,
      TRUE ~ 0),
    lambda = purrr::pmap_dbl(
      list(year, disturbance, disturbance_prev),
      function(year, disturbance, disturbance_prev) {
        lambda_ipm_year(
          year = year,
          disturbance = disturbance,
          disturbance_prev = disturbance_prev)
      }),
    scenario = factor(
      scenario,
      levels = c(
        'Undisturbed',
        'Current fire',
        'Previous fire',
        'Consecutive fire')))

lambda_year


# Plot year-specific scenarios -------------------------------------------------

fig_lambda_year <- ggplot(
  lambda_year,
  aes(x = year, y = lambda, color = scenario)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = c(
      'Undisturbed' = 'black',
      'Current fire' = 'red',
      'Previous fire' = 'blue',
      'Consecutive fire' = 'purple')) +
  theme_bw() +
  labs(
    title = 'Year-specific asymptotic lambda',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = expression(lambda),
    color = 'Disturbance scenario')

fig_lambda_year


# Observed and projected site-year population growth ---------------------------

# Site-year disturbance lookup -------------------------------------------------

fire_lookup_site <- df %>%
  group_by(site, year) %>%
  summarise(
    disturbance_num = max(
      as.numeric(as.character(disturbance)),
      na.rm = TRUE),
    disturbance_prev_num = max(
      as.numeric(as.character(disturbance_prev)),
      na.rm = TRUE),
    .groups = 'drop') %>%
  mutate(
    disturbance_num = if_else(
      is.finite(disturbance_num),
      disturbance_num,
      0),
    disturbance_prev_num = if_else(
      is.finite(disturbance_prev_num),
      disturbance_prev_num,
      0))


# Get site-year disturbance state ---------------------------------------------

get_disturbance_site_y <- function(site_i, yr, variable) {
  out <- fire_lookup_site %>%
    filter(site == site_i, year == yr) %>%
    pull({{ variable }})
  
  if (length(out) == 0 || is.na(out[1])) {
    return(0)
  }
  
  as.numeric(out[1])
}


# Observed site-year population counts -----------------------------------------

df_counts_site <- df %>%
  group_by(site, year) %>%
  summarise(
    n_sized = sum(is.finite(logvol_t0), na.rm = TRUE),
    n_unsized_recruits = sum(
      !is.finite(logvol_t0) & recruits == 1,
      na.rm = TRUE),
    n_ipm_state = n_sized + n_unsized_recruits,
    .groups = 'drop')


# Observed population growth ---------------------------------------------------

df_obs_pgr_site <- df_counts_site %>%
  rename(n_t0 = n_ipm_state) %>%
  left_join(
    df_counts_site %>%
      transmute(
        site,
        year = year - 1,
        n_t1 = n_ipm_state),
    by = c('site', 'year')) %>%
  filter(
    !is.na(n_t1),
    n_t0 > 0,
    n_t1 > 0) %>%
  left_join(
    fire_lookup_site,
    by = c('site', 'year')) %>%
  mutate(
    disturbance_num = replace_na(disturbance_num, 0),
    disturbance_prev_num = replace_na(disturbance_prev_num, 0),
    obs_pgr = n_t1 / n_t0)


# Construct observed initial size distribution --------------------------------

make_initial_n_site <- function(year0, site_i, pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  breaks <- seq(
    L,
    U,
    length.out = n + 1)
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  df0 <- df %>%
    filter(
      year == year0,
      site == site_i)
  
  rec_sized <- df0 %>%
    filter(
      recruits == 1,
      is.finite(logvol_t0)) %>%
    pull(logvol_t0)
  
  ad_sized <- df0 %>%
    filter(
      (is.na(recruits) | recruits != 1),
      is.finite(logvol_t0)) %>%
    pull(logvol_t0)
  
  rec_counts <- hist(
    pmin(pmax(rec_sized, L), U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE)$counts
  
  ad_counts <- hist(
    pmin(pmax(ad_sized, L), U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE)$counts
  
  rec_density <- rec_counts / h
  ad_density <- ad_counts / h
  
  n_unsized_recruits <- df0 %>%
    summarise(
      n = sum(
        !is.finite(logvol_t0) & recruits == 1,
        na.rm = TRUE)) %>%
    pull(n)
  
  if (n_unsized_recruits > 0) {
    rec_density <- rec_density +
      n_unsized_recruits * re_y_dist(
        y,
        pars,
        h = h)
  }
  
  c(rec_density, ad_density)
}


# Project one site-year --------------------------------------------------------

project_one_site_year <- function(yr, site_i) {
  site_re <- if (use_site_specific_recruitment) {
    site_i
  } else {
    NULL
  }
  
  disturbance_y <- get_disturbance_site_y(
    site_i,
    yr,
    disturbance_num)
  
  disturbance_prev_y <- get_disturbance_site_y(
    site_i,
    yr,
    disturbance_prev_num)
  
  pars_y <- make_ipm_pars(
    pars_mean = pars_all_mean,
    pars_year = pars_all_year,
    pars_site = pars_all_site,
    pars_re_year_t1 = pars_all_re_year_t1,
    year = yr,
    site = site_re,
    re_year_t1 = yr + 1)
  
  n_obs <- make_initial_n_site(
    year0 = yr,
    site_i = site_i,
    pars = pars_y)
  
  K <- kernel(
    pars = pars_y,
    disturbance = disturbance_y,
    disturbance_prev = disturbance_prev_y)$k_yx
  
  h <- (pars_y$U - pars_y$L) / pars_y$mat_siz
  n_initial <- sum(n_obs) * h
  
  if (n_initial <= 0) {
    return(data.frame(
      year = yr,
      site = site_i,
      disturbance_num = disturbance_y,
      disturbance_prev_num = disturbance_prev_y,
      n_obs_model = n_initial,
      n_proj_model = NA_real_,
      asym_lambda = NA_real_,
      proj_lambda = NA_real_))
  }
  
  n_proj <- K %*% n_obs
  n_projected <- sum(n_proj) * h
  
  data.frame(
    year = yr,
    site = site_i,
    disturbance_num = disturbance_y,
    disturbance_prev_num = disturbance_prev_y,
    n_obs_model = n_initial,
    n_proj_model = n_projected,
    asym_lambda = Re(eigen(K)$values[1]),
    proj_lambda = as.numeric(
      n_projected / n_initial))
}


# Run site-year projections ----------------------------------------------------

df_proj_site <- bind_rows(
  lapply(seq_len(nrow(df_obs_pgr_site)), function(i) {
    project_one_site_year(
      yr = df_obs_pgr_site$year[i],
      site_i = df_obs_pgr_site$site[i])
  }))


# Join observed and projected values -------------------------------------------

df_compare_site <- df_obs_pgr_site %>%
  left_join(
    df_proj_site,
    by = c(
      'year',
      'site',
      'disturbance_num',
      'disturbance_prev_num')) %>%
  mutate(
    error_asymptotic_vs_obs = asym_lambda - obs_pgr,
    error_projected_vs_obs = proj_lambda - obs_pgr)

df_compare_site %>%
  print(n = 100)


# Whole-population annual comparison ------------------------------------------

df_compare <- df_compare_site %>%
  group_by(year) %>%
  summarise(
    asym_lambda = weighted.mean(
      asym_lambda,
      w = n_obs_model,
      na.rm = TRUE),
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    obs_pgr = n_t1 / n_t0,
    n_obs_model = sum(n_obs_model, na.rm = TRUE),
    n_proj_model = sum(n_proj_model, na.rm = TRUE),
    proj_lambda = n_proj_model / n_obs_model,
    disturbance_num = as.integer(
      any(disturbance_num == 1, na.rm = TRUE)),
    disturbance_prev_num = as.integer(
      any(disturbance_prev_num == 1, na.rm = TRUE)),
    .groups = 'drop') %>%
  mutate(
    disturbance_scenario = case_when(
      disturbance_num == 0 & disturbance_prev_num == 0 ~
        'Undisturbed',
      disturbance_num == 1 & disturbance_prev_num == 0 ~
        'Current fire',
      disturbance_num == 0 & disturbance_prev_num == 1 ~
        'Previous fire',
      disturbance_num == 1 & disturbance_prev_num == 1 ~
        'Consecutive fire'),
    disturbance_scenario = factor(
      disturbance_scenario,
      levels = c(
        'Undisturbed',
        'Current fire',
        'Previous fire',
        'Consecutive fire')))


# Observed versus modeled data -------------------------------------------------

df_plot <- df_compare %>%
  select(
    year,
    obs_pgr,
    asym_lambda,
    proj_lambda,
    disturbance_scenario) %>%
  pivot_longer(
    cols = c(
      asym_lambda,
      proj_lambda),
    names_to = 'lambda_type',
    values_to = 'lambda') %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda =
        'Abundance-weighted asymptotic lambda',
      proj_lambda =
        'Projected lambda from observed site size distributions'))


# Observed versus modeled plot -------------------------------------------------

fig_mod_vs_obs <- ggplot(
  df_plot,
  aes(
    x = lambda,
    y = obs_pgr,
    color = disturbance_scenario)) +
  geom_point(size = 3) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2) +
  facet_wrap(
    ~ lambda_type,
    scales = 'free_x') +
  scale_color_manual(
    values = c(
      'Undisturbed' = 'black',
      'Current fire' = 'red',
      'Previous fire' = 'blue',
      'Consecutive fire' = 'purple')) +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = v_ggp_suffix,
    x = expression('Modeled ' * lambda),
    y = 'Observed population growth rate',
    color = 'Disturbance') +
  theme_classic()

fig_mod_vs_obs


# Log-transformed comparison ---------------------------------------------------

df_plot_log <- df_plot %>%
  filter(
    obs_pgr > 0,
    lambda > 0) %>%
  mutate(
    log_obs_pgr = log(obs_pgr),
    log_lambda = log(lambda))

fig_mod_vs_obs_log <- ggplot(
  df_plot_log,
  aes(
    x = log_lambda,
    y = log_obs_pgr,
    color = disturbance_scenario)) +
  geom_point(size = 3) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2) +
  facet_wrap(
    ~ lambda_type,
    scales = 'free_x') +
  scale_color_manual(
    values = c(
      'Undisturbed' = 'black',
      'Current fire' = 'red',
      'Previous fire' = 'blue',
      'Consecutive fire' = 'purple')) +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = paste(
      v_ggp_suffix,
      '- log-transformed lambda'),
    x = expression('log modeled ' * lambda),
    y = 'log observed population growth rate',
    color = 'Disturbance') +
  theme_classic()

fig_mod_vs_obs_log


# Summary statistics -----------------------------------------------------------

df_compare_summary <- df_compare %>%
  filter(
    obs_pgr > 0,
    asym_lambda > 0,
    proj_lambda > 0) %>%
  summarise(
    n_years = n(),
    arithmetic_mean_obs_pgr = mean(
      obs_pgr,
      na.rm = TRUE),
    geometric_mean_obs_pgr = exp(
      mean(log(obs_pgr), na.rm = TRUE)),
    arithmetic_mean_asym_lambda = mean(
      asym_lambda,
      na.rm = TRUE),
    geometric_mean_asym_lambda = exp(
      mean(log(asym_lambda), na.rm = TRUE)),
    arithmetic_mean_proj_lambda = mean(
      proj_lambda,
      na.rm = TRUE),
    geometric_mean_proj_lambda = exp(
      mean(log(proj_lambda), na.rm = TRUE)),
    mean_error_asymptotic_vs_obs = mean(
      asym_lambda - obs_pgr,
      na.rm = TRUE),
    mean_error_projected_vs_obs = mean(
      proj_lambda - obs_pgr,
      na.rm = TRUE),
    rmse_asymptotic_vs_obs = sqrt(
      mean(
        (asym_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    rmse_projected_vs_obs = sqrt(
      mean(
        (proj_lambda - obs_pgr)^2,
        na.rm = TRUE)))

df_compare_summary

df_compare %>%
  print(n = 100)


# Site-level summary -----------------------------------------------------------

df_compare_site_summary <- df_compare_site %>%
  filter(
    !is.na(obs_pgr),
    !is.na(asym_lambda),
    !is.na(proj_lambda),
    obs_pgr > 0,
    asym_lambda > 0,
    proj_lambda > 0) %>%
  group_by(site) %>%
  summarise(
    n_year_transitions = n(),
    lambda_obs_geometric = exp(
      mean(log(obs_pgr), na.rm = TRUE)),
    lambda_obs_arithmetic = mean(
      obs_pgr,
      na.rm = TRUE),
    lambda_asymptotic_geometric = exp(
      mean(log(asym_lambda), na.rm = TRUE)),
    lambda_asymptotic_arithmetic = mean(
      asym_lambda,
      na.rm = TRUE),
    lambda_projected_geometric = exp(
      mean(log(proj_lambda), na.rm = TRUE)),
    lambda_projected_arithmetic = mean(
      proj_lambda,
      na.rm = TRUE),
    error_projected_geo_vs_obs_geo =
      lambda_projected_geometric -
      lambda_obs_geometric,
    rmse_projected_vs_obs = sqrt(
      mean(
        (proj_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    mean_n_initial = mean(
      n_obs_model,
      na.rm = TRUE),
    mean_current_disturbance = mean(
      disturbance_num,
      na.rm = TRUE),
    mean_previous_disturbance = mean(
      disturbance_prev_num,
      na.rm = TRUE),
    .groups = 'drop') %>%
  arrange(
    as.numeric(as.character(site)))

df_compare_site_summary %>%
  print(
    n = 100,
    width = Inf)
 


df_plot_site <- df_compare_site %>%
  mutate(
    disturbance_scenario = case_when(
      disturbance_num == 0 & disturbance_prev_num == 0 ~ "Undisturbed",
      disturbance_num == 1 & disturbance_prev_num == 0 ~ "Current fire",
      disturbance_num == 0 & disturbance_prev_num == 1 ~ "Previous fire",
      disturbance_num == 1 & disturbance_prev_num == 1 ~ "Consecutive fire"),
    disturbance_scenario = factor(
      disturbance_scenario,
      levels = c(
        "Undisturbed", "Current fire",
        "Previous fire", "Consecutive fire"))) %>%
  select(
    year, site, obs_pgr, asym_lambda, proj_lambda,
    disturbance_scenario) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = "lambda_type",
    values_to = "lambda") %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = "Site-year asymptotic lambda",
      proj_lambda = "Site-year projected lambda"))

fig_mod_vs_obs_site <- ggplot(
  df_plot_site,
  aes(
    x = lambda,
    y = obs_pgr,
    color = disturbance_scenario)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~ lambda_type, scales = "free_x") +
  scale_color_manual(
    values = c(
      "Undisturbed" = "black",
      "Current fire" = "red",
      "Previous fire" = "blue",
      "Consecutive fire" = "purple")) +
  labs(
    title = "Observed population growth vs site-year modeled lambda",
    subtitle = v_ggp_suffix,
    x = expression("Modeled " * lambda),
    y = "Observed population growth rate",
    color = "Disturbance",
    shape = "Site") +
  theme_classic()

fig_mod_vs_obs_site



fig_mod_vs_obs_site_log <- ggplot(
  df_plot_site %>% filter(obs_pgr > 0, lambda > 0),
  aes(
    x = lambda,
    y = obs_pgr,
    color = disturbance_scenario)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~ lambda_type) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Undisturbed" = "black",
      "Current fire" = "red",
      "Previous fire" = "blue",
      "Consecutive fire" = "purple")) +
  labs(
    title = "Observed population growth vs site-year modeled lambda",
    subtitle = v_ggp_suffix,
    x = expression("Modeled " * lambda * " (log scale)"),
    y = "Observed population growth rate (log scale)",
    color = "Disturbance",
    shape = "Site") +
  theme_classic()

fig_mod_vs_obs_site_log

# # Save key outputs -------------------------------------------------------------
# # Uncomment if needed.
# # saveRDS(pars_all_mean,
# #         file.path(dir_result, 'pole_year_fire_pars_mean.rds'))
# # saveRDS(pars_all_year,
# #         file.path(dir_result, 'pole_year_fire_pars_year.rds'))
# # saveRDS(df_compare,
# #         file.path(dir_result, 'pole_year_fire_lambda_compare.rds'))
# # saveRDS(df_compare_site,
# #         file.path(dir_result, 'pole_year_fire_lambda_compare_site.rds'))
