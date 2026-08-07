# IPM year-specific with fire - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.10

# Study organism: Polygala lewtonii
# Citing publication: https://doi.org/10.1071/BT11271
# Time period: 2001-2024, through 2025


################################################################################
################################################################################

# ATTENTION: This runs over multiple hours

################################################################################
################################################################################


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
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
v_mod_set_fl   <- c()
v_mod_set_fl_n <- c()
v_mod_set_re   <- c()


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



# Load workdata ---------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0('ab_', v_sp_abb, '_df_workdata_260519.csv'))) %>%
  mutate(
    disturbance = as.numeric(dist_transition),
    disturbance_prev = as.numeric(disturbance_prev),
    year = as.integer(year),
    site = factor(site)) %>%
  filter(
    !is.na(year),
    !(year %in% v_years_re))



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
    size_t0 != 0,
    is.finite(logvol_t0),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = factor(disturbance, levels = c(0, 1))) %>%
  dplyr::select(
    id, site, year, size_t0, survives, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3,
    disturbance, stage)


# Survival model ---------------------------------------------------------------
mod_su_0 <- glmer(
  survives ~ disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_10 <- glmer(
  survives ~ logvol_t0 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_11 <- glmer(
  survives ~ logvol_t0 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_20 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_21 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_30 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)
mod_su_31 <- glmer(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance +
    (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)


mods_su <- list(
  mod_su_0, mod_su_10, mod_su_11, 
  mod_su_20, mod_su_21, mod_su_30, mod_su_31)
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
  
  x <- seq(
    min(df_i$logvol_t0, na.rm = TRUE),
    max(df_i$logvol_t0, na.rm = TRUE), length.out = 100)
  
  pred_i <- tidyr::expand_grid(
    logvol_t0 = x,
    disturbance = levels(df_su$disturbance)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      disturbance = factor(
        disturbance, levels = levels(df_su$disturbance)),
      year = factor(year_i, levels = levels(df_su$year)))
  
  pred_i$survives <- predict(
    mod_su_best, newdata = pred_i,
    type = "response", re.form = NULL)
  
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, disturbance) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(), .groups = "drop") %>%
    filter(!is.na(logvol_t0), n > 0)
  
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, survives, color = disturbance), size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, survives, color = disturbance,
          linetype = disturbance),
      linewidth = 0.7) +
    scale_color_manual(values = c("0" = "black", "1" = "red")) +
    labs(
      title = year_i,
      x = expression("log(volume)"[t0]),
      y = expression("Survival probability"[t1])) +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = "none")
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
    size_t0 != 0,
    size_t1 != 0,
    is.finite(logvol_t0),
    is.finite(logvol_t1),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = factor(
      disturbance,
      levels = c(0, 1))) %>%
  dplyr::select(
    id, site, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3,
    disturbance, stage)

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
mod_gr_10 <- lmer(
  logvol_t1 ~ logvol_t0 + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_11 <- lmer(
  logvol_t1 ~ logvol_t0 + disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_20 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_21 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + disturbance +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_30 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)
mod_gr_31 <- lmer(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 +
    disturbance + (logvol_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)


mods_gr <- list(
  mod_gr_0, mod_gr_10, mod_gr_11, mod_gr_20, mod_gr_21, 
  mod_gr_30, mod_gr_31)
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
    disturbance = levels(df_gr$disturbance)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      disturbance = factor(
        disturbance, levels = levels(df_gr$disturbance)),
      year = factor(year_i, levels = levels(df_gr$year)))
  
  pred_i <- pred_i %>%
    mutate(
      logvol_t1 = predict(
        mod_gr_best,
        newdata = pred_i,
        re.form = NULL))
  pts_i <- df_i %>%
    mutate(bin = cut(logvol_t0, breaks = 8)) %>%
    group_by(bin, disturbance) %>%
    summarise(
      logvol_t0 = mean(logvol_t0, na.rm = TRUE),
      logvol_t1 = mean(logvol_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logvol_t0), n > 0)
  ggplot() +
    geom_point(
      data = pts_i,
      aes(logvol_t0, logvol_t1, color = disturbance),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logvol_t0, logvol_t1, color = disturbance,
          linetype = disturbance),
      linewidth = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', lty = 2) +
    scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
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
    is.finite(logvol_t0),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance_prev)) %>%
  mutate(
    year = factor(year),
    disturbance_prev = factor(
      disturbance_prev,
      levels = c(0, 1))) %>%
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
# Total flowering stems in quadrat q during year t predict recruits in the
# same quadrat during year t + 1. Current fire is retained at the quadrat level.
df_quad_year_sampled <- df %>%
  distinct(
    site,
    quad,
    year)

df_fs2r <- df %>%
  group_by(
    site,
    quad,
    year) %>%
  summarise(
    fs_t0 = sum(
      fl_nr,
      na.rm = TRUE),
    disturbance = first(
      disturbance[!is.na(disturbance)],
      default = NA_real_),
    .groups = 'drop') %>%
  mutate(
    year_t1 = year + 1) %>%
  semi_join(
    df_quad_year_sampled %>%
      transmute(
        site,
        quad,
        year_t1 = year),
    by = c(
      'site',
      'quad',
      'year_t1')) %>%
  left_join(
    df %>%
      filter(recruits == 1) %>%
      count(
        site,
        quad,
        year,
        name = 're_t1'),
    by = c(
      'site',
      'quad',
      'year_t1' = 'year')) %>%
  mutate(
    re_t1 = replace_na(re_t1, 0L),
    fs_t0_log = log1p(fs_t0),
    re_t1_log = log1p(re_t1)) %>%
  filter(
    !is.na(year),
    !is.na(year_t1),
    !is.na(disturbance),
    !(year %in% v_years_re),
    !(year_t1 %in% v_years_re)) %>%
  mutate(
    site = factor(site),
    year_t0 = factor(year),
    year_t1 = factor(year_t1),
    disturbance_t0 = factor(
      disturbance,
      levels = c(0, 1)))

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
    x = 'Total flowering stems in quadrat in year t',
    y = 'Number of recruits in quadrat in year t + 1',
    color = 'Fire')

fig_fs2r_raw


# Flowering-stem reference by previous-year fire -------------------------------
df_fs2r_previous_fire <- df_fs2r %>%
  left_join(
    df %>%
      distinct(
        site,
        quad,
        year,
        disturbance_prev),
    by = c(
      'site',
      'quad',
      'year')) %>%
  filter(!is.na(disturbance_prev))

df_stem_reference <- df_fs2r_previous_fire %>%
  group_by(disturbance_prev) %>%
  summarise(
    fs_t0_ref = mean(fs_t0),
    n_quadrat_transitions = n(),
    .groups = 'drop')

df_stem_reference


# Flower 2 Recruit model -------------------------------------------------------
ctrl_re <- lmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

df_fs2r <- df_fs2r %>%
  mutate(
    disturbance =
      as.numeric(as.character(disturbance_t0)))

# Intercept only
mod_re_00 <- lmer(
  re_t1_log ~ 1 +
    (1 | year_t1),
  data = df_fs2r,
  REML = FALSE,
  control = ctrl_re)

# Disturbance only
mod_re_0 <- lmer(
  re_t1_log ~ disturbance +
    (1 | year_t1),
  data = df_fs2r,
  REML = FALSE,
  control = ctrl_re)

# Flowering stems only
mod_re_10 <- lmer(
  re_t1_log ~ fs_t0_log +
    (1 | year_t1),
  data = df_fs2r,
  REML = FALSE,
  control = ctrl_re)

# Flowering stems + disturbance
mod_re_1 <- lmer(
  re_t1_log ~ fs_t0_log + disturbance +
    (1 | year_t1),
  data = df_fs2r,
  REML = FALSE,
  control = ctrl_re)

# Flowering stems x disturbance
mod_re_20 <- lmer(
  re_t1_log ~ fs_t0_log * disturbance +
    (1 | year_t1),
  data = df_fs2r,
  REML = FALSE,
  control = ctrl_re)

mods_re <- list(
  mod_re_00,
  mod_re_0,
  mod_re_10,
  mod_re_1,
  mod_re_20)

mods_re_dAIC <- bbmle::AICctab(
  mods_re,
  weights = TRUE,
  sort = FALSE)$dAIC

mods_re_sorted <- order(mods_re_dAIC)

if (length(v_mod_set_re) == 0) {
  mod_re_index_bestfit <- mods_re_sorted[1]
} else {
  mod_re_index_bestfit <- v_mod_set_re + 1
}

v_mod_re_index <- mod_re_index_bestfit - 1
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
  
  x <- seq(
    0, max(df_i$fs_t0_log, na.rm = TRUE), length.out = 100)
  
  pred_i <- tidyr::expand_grid(
    fs_t0_log = x,
    disturbance_t0 = levels(df_fs2r$disturbance_t0)) %>%
    mutate(
      disturbance =
        as.numeric(as.character(disturbance_t0)),
      year_t1 = factor(year_i, levels = levels(df_fs2r$year_t1)),
      disturbance_t0 = factor(
        disturbance_t0, levels = levels(df_fs2r$disturbance_t0)))
  
  pred_i <- pred_i %>%
    mutate(
      pred_re = predict_re_safe(
        mod_re_best, pred_i, re.form = NULL))
  
  ggplot() +
    geom_point(
      data = df_i,
      aes(fs_t0_log, re_t1_log, color = site, shape = site),
      alpha = 0.75, size = 2) +
    geom_line(
      data = pred_i,
      aes(
        fs_t0_log, pred_re,
        group = disturbance_t0,
        linetype = disturbance_t0),
      inherit.aes = FALSE, color = 'black', linewidth = 0.8) +
    scale_color_viridis_d(option = 'turbo', end = 0.95) +
    scale_shape_manual(values = site_shapes_re) +
    scale_linetype_manual(values = c('0' = 'solid', '1' = 'dashed')) +
    theme_bw() +
    labs(
      title = year_i,
      x = expression('log(1 + flowering stems)'[t0]),
      y = expression('log(1 + recruits)'[t1]),
      color = 'Site',
      shape = 'Site',
      linetype = 'Fire') +
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
  bf = c('disturbance1', 'regex:^disturbance'))

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
  re_bd = c('disturbance'),
  re_bf = c(
    'fs_t0_log:disturbance',
    'disturbance:fs_t0_log'))


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
  'fs_t0_ref_0', df_stem_reference %>% filter(disturbance_prev == 0) %>% pull(fs_t0_ref),
  'fs_t0_ref_1', df_stem_reference %>% filter(disturbance_prev == 1) %>% pull(fs_t0_ref),
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


# Recruitment-year random-effect deviations ---------------------------
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
    pars_re_year_t1 = NULL,
    year = NULL,
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
  
  if (is.null(re_year_t1) && !is.null(year)) {
    re_year_t1 <- as.numeric(as.character(year)) + 1
  }
  
  if (!is.null(re_year_t1) && !is.null(pars_re_year_t1)) {
    re_year_t1 <- as.character(re_year_t1)
    pars$re_b0 <- get_par(pars, 're_b0', 0) +
      get_par(
        pars_re_year_t1,
        paste0('re_u0_year_t1_', re_year_t1), 0)
    pars$re_b1 <- get_par(pars, 're_b1', 0) +
      get_par(
        pars_re_year_t1,
        paste0('re_u1_year_t1_', re_year_t1), 0)
  }
  
  pars
}

surv_lp <- function(x, pars, disturbance = 0) {
  get_par(pars, 'surv_b0', 0) +
    get_par(pars, 'surv_b1', 0) * x +
    get_par(pars, 'surv_b2', 0) * x^2 +
    get_par(pars, 'surv_b3', 0) * x^3 +
    get_par(pars, 'surv_bf', 0) * disturbance
}

sx <- function(x, pars, disturbance = 0) {
  inv_logit(surv_lp(x, pars, disturbance))
}

grow_mu <- function(x, pars, disturbance = 0) {
  get_par(pars, 'grow_b0', 0) +
    get_par(pars, 'grow_b1', 0) * x +
    get_par(pars, 'grow_b2', 0) * x^2 +
    get_par(pars, 'grow_b3', 0) * x^3 +
    get_par(pars, 'grow_bf', 0) * disturbance
}

grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * x))
}

gxy <- function(x, y, pars, disturbance = 0) {
  dnorm(
    y,
    mean = grow_mu(x, pars, disturbance),
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

get_fs_ref <- function(
    pars,
    disturbance_prev = 0) {
  
  if (disturbance_prev == 1) {
    get_par(pars, 'fs_t0_ref_1')
  } else {
    get_par(pars, 'fs_t0_ref_0')
  }
}

re_total_ref <- function(
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  
  re_intercept <- get_par(pars, 're_b0', 0) +
    get_par(pars, 're_bd', 0) * disturbance
  
  re_slope <- get_par(pars, 're_b1', 0) +
    get_par(pars, 're_bf', 0) * disturbance
  
  fs_ref <- get_fs_ref(
    pars,
    disturbance_prev = disturbance_prev)
  
  re_log <- re_intercept +
    re_slope * log1p(fs_ref)
  
  re_value <- expm1(
    re_log +
      0.5 * get_par(pars, 're_sigma', 0)^2)
  
  pmax(re_value, 0)
}

rx <- function(
    x,
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  
  fs_value <- fs_x(
    x,
    pars,
    disturbance_prev = disturbance_prev)
  
  fs_ref <- get_fs_ref(
    pars,
    disturbance_prev = disturbance_prev)
  
  re_ref <- re_total_ref(
    pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  re_per_fs <- re_ref /
    pmax(
      fs_ref,
      .Machine$double.eps)
  
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


# Single-state size-structured kernel ------------------------------------------
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
  
  S <- sx(y, pars, disturbance = disturbance)
  
  G <- t(outer(
    y, y, Vectorize(function(x, y) {
      gxy(x, y, pars, disturbance = disturbance)
    }))) * h
  
  for (i in 1:n) {
    if (i <= n / 2) {
      G[1, i] <- G[1, i] + 1 - sum(G[, i])
    } else {
      G[n, i] <- G[n, i] + 1 - sum(G[, i])
    }
  }
  
  T <- sweep(G, 2, S, '*')
  
  F <- fyx(
    y = y,
    x = y,
    pars = pars,
    h = h,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  k_yx <- T + F
  
  list(
    k_yx = k_yx,
    T = T,
    F = F,
    G = G,
    S = S,
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
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 0)
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 1)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 1)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 0)



# Year-specific IPMs -----------------------------------------------------------
lambda_ipm_year <- function(
    year,
    disturbance = 0,
    disturbance_prev = 0,
    re_year_t1 = NULL) {
  pars_i <- make_ipm_pars(
    pars_mean = pars_all_mean,
    pars_year = pars_all_year,
    pars_re_year_t1 = pars_all_re_year_t1,
    year = year,
    re_year_t1 = re_year_t1)
  
  lambda_ipm(
    pars_i,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
}

# Year-specific disturbance scenarios ----

ipm_years <- sort(unique(as.integer(as.character(df_su$year))))

lambda_year <- data.frame(
  year = ipm_years,
  lambda_undisturbed = sapply(ipm_years, function(i) {
    lambda_ipm_year(
      year = i,
      disturbance = 0,
      disturbance_prev = 0)}),
  lambda_current_fire = sapply(ipm_years, function(i) {
    lambda_ipm_year(
      year = i,
      disturbance = 1,
      disturbance_prev = 0)}),
  lambda_previous_fire = sapply(ipm_years, function(i) {
    lambda_ipm_year(
      year = i,
      disturbance = 0,
      disturbance_prev = 1)}),
  lambda_consecutive_fire = sapply(ipm_years, function(i) {
    lambda_ipm_year(
      year = i,
      disturbance = 1,
      disturbance_prev = 1)})
)

fig_lambda_year <- lambda_year %>%
  pivot_longer(
    cols = -year,
    names_to = 'disturbance',
    values_to = 'lambda') %>%
  mutate(
    disturbance = recode(
      disturbance,
      lambda_undisturbed = 'Undisturbed',
      lambda_current_fire = 'Current-year fire',
      lambda_previous_fire = 'Previous-year fire',
      lambda_consecutive_fire = 'Fire in both years')) %>%
  ggplot(aes(year, lambda, color = disturbance)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(
    title = 'Year-specific asymptotic lambda',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = expression(lambda),
    color = 'Disturbance')

fig_lambda_year


# Observed and projected population growth ------------------------------------
# Observed growth is calculated from the same quadrats represented in both
# consecutive years. Each quadrat is projected with its own current- and
# previous-year disturbance state.


# Quadrat-year disturbance lookup ---------------------------------------------
fire_lookup_quad <- df %>%
  group_by(
    site,
    quad,
    year) %>%
  summarise(
    disturbance_num = first(
      disturbance[!is.na(disturbance)],
      default = NA_real_),
    disturbance_prev_num = first(
      disturbance_prev[!is.na(disturbance_prev)],
      default = NA_real_),
    .groups = 'drop')


# Quadrat-year abundance -------------------------------------------------------
# Counts include sized plants and recruits without a usable size measurement.
df_counts_quad_year <- df %>%
  filter(
    !is.na(site),
    !is.na(quad),
    !is.na(year)) %>%
  group_by(
    site,
    quad,
    year) %>%
  summarise(
    n_sized = sum(
      is.finite(logvol_t0),
      na.rm = TRUE),
    n_unsized_recruits = sum(
      !is.finite(logvol_t0) & recruits == 1,
      na.rm = TRUE),
    n_ipm_state = n_sized + n_unsized_recruits,
    .groups = 'drop')


# Matched consecutive quadrat transitions -------------------------------------
df_obs_pgr_quad <- df_counts_quad_year %>%
  arrange(
    site,
    quad,
    year) %>%
  group_by(
    site,
    quad) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n_ipm_state),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(
    year_gap == 1,
    n_ipm_state > 0) %>%
  rename(
    n_t0 = n_ipm_state) %>%
  left_join(
    fire_lookup_quad,
    by = c(
      'site',
      'quad',
      'year')) %>%
  filter(
    !is.na(disturbance_num),
    !is.na(disturbance_prev_num)) %>%
  mutate(
    obs_pgr = n_t1 / n_t0)


# Initial quadrat size distribution -------------------------------------------
make_initial_n_quad <- function(
    year0,
    site_i,
    quad_i,
    pars) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  breaks <- seq(
    L,
    U,
    length.out = n + 1)
  
  b <- L + c(0:n) * h
  
  y <- 0.5 *
    (b[1:n] + b[2:(n + 1)])
  
  df0 <- df %>%
    filter(
      year == year0,
      site == site_i,
      quad == quad_i)
  
  sizes <- df0 %>%
    filter(
      is.finite(logvol_t0)) %>%
    pull(logvol_t0)
  
  size_counts <- hist(
    pmin(
      pmax(sizes, L),
      U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE)$counts
  
  n_density <- size_counts / h
  
  n_unsized_recruits <- df0 %>%
    summarise(
      n = sum(
        !is.finite(logvol_t0) &
          recruits == 1,
        na.rm = TRUE)) %>%
    pull(n)
  
  if (n_unsized_recruits > 0) {
    n_density <- n_density +
      n_unsized_recruits *
      re_y_dist(
        y,
        pars,
        h = h)
  }
  
  n_density
}


# Project one quadrat-year -----------------------------------------------------
project_one_quad_year <- function(
    yr,
    site_i,
    quad_i,
    disturbance_y,
    disturbance_prev_y) {
  
  pars_y <- make_ipm_pars(
    pars_mean = pars_all_mean,
    pars_year = pars_all_year,
    pars_re_year_t1 = pars_all_re_year_t1,
    year = yr,
    re_year_t1 = yr + 1)
  
  n_obs <- make_initial_n_quad(
    year0 = yr,
    site_i = site_i,
    quad_i = quad_i,
    pars = pars_y)
  
  K <- kernel(
    pars = pars_y,
    disturbance = disturbance_y,
    disturbance_prev = disturbance_prev_y)$k_yx
  
  h <- (pars_y$U - pars_y$L) /
    pars_y$mat_siz
  
  n_initial <- sum(n_obs) * h
  
  if (n_initial <= 0) {
    return(
      data.frame(
        year = yr,
        site = site_i,
        quad = quad_i,
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
    quad = quad_i,
    disturbance_num = disturbance_y,
    disturbance_prev_num = disturbance_prev_y,
    n_obs_model = n_initial,
    n_proj_model = n_projected,
    asym_lambda = Re(
      eigen(K)$values[1]),
    proj_lambda = as.numeric(
      n_projected / n_initial))
}


# Project all matched quadrat transitions -------------------------------------
df_proj_quad <- bind_rows(
  lapply(
    seq_len(nrow(df_obs_pgr_quad)),
    function(i) {
      
      project_one_quad_year(
        yr = df_obs_pgr_quad$year[i],
        site_i = df_obs_pgr_quad$site[i],
        quad_i = df_obs_pgr_quad$quad[i],
        disturbance_y =
          df_obs_pgr_quad$disturbance_num[i],
        disturbance_prev_y =
          df_obs_pgr_quad$disturbance_prev_num[i])
    }))


# Quadrat-level observed and modeled growth -----------------------------------
df_compare_quad <- df_obs_pgr_quad %>%
  left_join(
    df_proj_quad,
    by = c(
      'year',
      'site',
      'quad',
      'disturbance_num',
      'disturbance_prev_num')) %>%
  mutate(
    error_asymptotic_vs_obs =
      asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr)


# Site-year comparison from matched quadrats ----------------------------------
df_compare_site <- df_compare_quad %>%
  group_by(
    site,
    year) %>%
  summarise(
    asym_lambda = weighted.mean(
      asym_lambda,
      w = n_obs_model,
      na.rm = TRUE),
    n_t0 = sum(
      n_t0,
      na.rm = TRUE),
    n_t1 = sum(
      n_t1,
      na.rm = TRUE),
    n_obs_model = sum(
      n_obs_model,
      na.rm = TRUE),
    n_proj_model = sum(
      n_proj_model,
      na.rm = TRUE),
    n_quadrats = n(),
    n_burned_quadrats = sum(
      disturbance_num == 1,
      na.rm = TRUE),
    disturbance_num = as.numeric(
      any(
        disturbance_num == 1,
        na.rm = TRUE)),
    disturbance_prev_num = as.numeric(
      any(
        disturbance_prev_num == 1,
        na.rm = TRUE)),
    .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 / n_t0,
    proj_lambda = n_proj_model / n_obs_model,
    error_asymptotic_vs_obs =
      asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr)

df_compare_site %>%
  print(
    n = 100,
    width = Inf)


# Whole-population annual comparison ------------------------------------------
# Annual totals contain only quadrats represented in both t and t + 1.
df_compare <- df_compare_quad %>%
  group_by(year) %>%
  summarise(
    asym_lambda = weighted.mean(
      asym_lambda,
      w = n_obs_model,
      na.rm = TRUE),
    n_t0 = sum(
      n_t0,
      na.rm = TRUE),
    n_t1 = sum(
      n_t1,
      na.rm = TRUE),
    n_obs_model = sum(
      n_obs_model,
      na.rm = TRUE),
    n_proj_model = sum(
      n_proj_model,
      na.rm = TRUE),
    n_quadrats = n(),
    n_sites = n_distinct(site),
    n_burned_quadrats = sum(
      disturbance_num == 1,
      na.rm = TRUE),
    p_burned_quadrats = mean(
      disturbance_num == 1,
      na.rm = TRUE),
    disturbance = if_else(
      any(
        disturbance_num == 1,
        na.rm = TRUE),
      'Fire',
      'No fire'),
    .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 / n_t0,
    proj_lambda = n_proj_model / n_obs_model,
    disturbance = factor(
      disturbance,
      levels = c(
        'No fire',
        'Fire')))

# Observed vs modeled plot

df_plot <- df_compare %>%
  select(year, obs_pgr, asym_lambda, proj_lambda, disturbance) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = 'lambda_type',
    values_to = 'lambda') %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = 'Abundance-weighted asymptotic lambda',
      proj_lambda = 'Projected lambda from observed quadrat size distributions'))

fig_mod_vs_obs <- ggplot(
  df_plot,
  aes(x = lambda, y = obs_pgr, color = disturbance)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = 'free_x') +
  scale_color_manual(values = c('No fire' = 'black', 'Fire' = 'red')) +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = v_ggp_suffix,
    x = expression('Modeled ' * lambda),
    y = 'Observed population growth rate',
    color = 'Fire') +
  theme_classic()

fig_mod_vs_obs

# Log-transformed observed vs modeled plot

df_plot_log <- df_plot %>%
  filter(obs_pgr > 0, lambda > 0) %>%
  mutate(
    log_obs_pgr = log(obs_pgr),
    log_lambda = log(lambda))

fig_mod_vs_obs_log <- ggplot(
  df_plot_log,
  aes(x = log_lambda, y = log_obs_pgr, color = disturbance)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = 'free_x') +
  scale_color_manual(values = c('No fire' = 'black', 'Fire' = 'red')) +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = paste(v_ggp_suffix, '- log-transformed lambda'),
    x = expression('log modeled ' * lambda),
    y = 'log observed population growth rate',
    color = 'Fire') +
  theme_classic()

fig_mod_vs_obs_log


# Summary statistics -----------------------------------------------------------
df_compare_summary <- df_compare %>%
  summarise(
    n_years = n(),
    arithmetic_mean_obs_pgr = mean(obs_pgr, na.rm = TRUE),
    geometric_mean_obs_pgr = exp(mean(log(obs_pgr), na.rm = TRUE)),
    arithmetic_mean_asym_lambda = mean(asym_lambda, na.rm = TRUE),
    geometric_mean_asym_lambda = exp(mean(log(asym_lambda), na.rm = TRUE)),
    arithmetic_mean_proj_lambda = mean(proj_lambda, na.rm = TRUE),
    geometric_mean_proj_lambda = exp(mean(log(proj_lambda), na.rm = TRUE)),
    mean_error_asymptotic_vs_obs = mean(
      asym_lambda - obs_pgr, na.rm = TRUE),
    mean_error_projected_vs_obs = mean(
      proj_lambda - obs_pgr, na.rm = TRUE),
    percent_bias_asymptotic_vs_obs = 100 * sum(
      asym_lambda - obs_pgr, na.rm = TRUE) /
      sum(obs_pgr, na.rm = TRUE),
    percent_bias_projected_vs_obs = 100 * sum(
      proj_lambda - obs_pgr, na.rm = TRUE) /
      sum(obs_pgr, na.rm = TRUE),
    rmse_asymptotic_vs_obs = sqrt(mean(
      (asym_lambda - obs_pgr)^2, na.rm = TRUE)),
    rmse_projected_vs_obs = sqrt(mean(
      (proj_lambda - obs_pgr)^2, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "statistic",
    values_to = "value")

df_compare_summary

df_compare %>% print(n = 100)


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
    lambda_obs_geometric = exp(mean(log(obs_pgr), na.rm = TRUE)),
    lambda_obs_arithmetic = mean(obs_pgr, na.rm = TRUE),
    lambda_asymptotic_geometric = exp(mean(log(asym_lambda),
                                           na.rm = TRUE)),
    lambda_asymptotic_arithmetic = mean(asym_lambda, na.rm = TRUE),
    lambda_projected_geometric = exp(mean(log(proj_lambda),
                                          na.rm = TRUE)),
    lambda_projected_arithmetic = mean(proj_lambda, na.rm = TRUE),
    error_projected_geo_vs_obs_geo =
      lambda_projected_geometric - lambda_obs_geometric,
    rmse_projected_vs_obs = sqrt(mean((proj_lambda - obs_pgr)^2,
                                      na.rm = TRUE)),
    mean_n_initial = mean(n_obs_model, na.rm = TRUE),
    mean_disturbance = mean(disturbance_num, na.rm = TRUE),
    .groups = 'drop') %>%
  arrange(as.numeric(as.character(site)))

df_compare_site_summary %>% print(n = 100, width = Inf)

# 
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
