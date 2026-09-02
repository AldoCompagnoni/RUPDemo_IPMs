# IPM year-specific, environmental and historical, with fire
# Archbold - Eriogonum longifolium

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.08.31

# Study organism: Eriogonum longifolium var. gnaphalifolium
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.226.1
# Meta data link:
# https://portal.edirepository.org/nis/metadataviewer?packageid=edi.226.1
# Citing publication: Satterthwaite et al. 2002, Ecological Applications
# Time period: 1990-2013


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


# Specification ---------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Eriogonum longifolium')
custom_delimiter <- c()
v_years_re <- c()

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Keep manual model choices empty for AICc selection.
# For survival, growth, dormancy, flowering, and scape number, values 0:3
# restrict selection to the corresponding polynomial degree.
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
v_mod_set_do   <- c()
v_mod_set_fl   <- c()
v_mod_set_fl_n <- c()

# Recruitment candidates are indexed 0:8 below.
v_mod_set_re <- c()


# Directory -------------------------------------------------------------------
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


# Functions -------------------------------------------------------------------
source('helper_functions/plot_binned_prop.R')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data ------------------------------------------------------------------------
disturbance_case_levels <- c(
  'Undisturbed',
  'Current-year fire',
  'Previous-year fire',
  'Fire in both years')

disturbance_case_cols <- c(
  'Undisturbed' = 'black',
  'Current-year fire' = 'red',
  'Previous-year fire' = 'purple',
  'Fire in both years' = 'magenta')

df <- read.csv(
  file.path(
    dir_data,
    paste0('ab_', v_sp_abb, '_df_workdata_260820.csv'))) %>%
  mutate(
    disturbance = as.numeric(dist_transition),
    disturbance_prev = as.numeric(disturbance_prev),
    year = as.integer(year),
    row_type = as.character(row_type),
    disturbance_case = factor(
      paste0(disturbance, disturbance_prev),
      levels = c('00', '10', '01', '11'),
      labels = disturbance_case_levels)) %>%
  filter(
    !is.na(year),
    !(year %in% v_years_re))

df_ind <- df %>%
  filter(row_type == 'individual')

df_re_all <- df %>%
  filter(row_type == 'recruitment')

df_re <- df_re_all %>%
  filter(
    recruitment_complete,
    !is.na(recruits_simple))


# Controls --------------------------------------------------------------------
ctrl_glmer <- glmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

ctrl_lmer <- lmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

ctrl_re <- glmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))


# Survival data ---------------------------------------------------------------
df_su <- df %>%
  filter(
    state == 'active',
    !is.na(survives),
    size_t0 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t0_2),
    is.finite(logsize_t0_3),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    id, year, size_t0, survives, disturbance,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Survival model --------------------------------------------------------------
mod_su_00 <- glmer(
  survives ~ 1 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_0 <- glmer(
  survives ~ disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_10 <- glmer(
  survives ~ logsize_t0 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_1 <- glmer(
  survives ~ logsize_t0 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_20 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_30 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    disturbance + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mods_su <- list(
  mod_su_00, mod_su_0, mod_su_10, mod_su_1,
  mod_su_20, mod_su_2, mod_su_30, mod_su_3)

mods_su_dAICc <- bbmle::AICctab(
  mods_su, weights = TRUE, sort = FALSE)$dAICc

if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- which.min(mods_su_dAICc)
} else {
  keep <- 2 * v_mod_set_su + c(1, 2)
  mod_su_index_bestfit <- keep[which.min(mods_su_dAICc[keep])]
}

v_mod_su_index <- floor((mod_su_index_bestfit - 1) / 2)
mod_su_bestfit <- mods_su[[mod_su_index_bestfit]]

mod_su_bestfit
summary(mod_su_bestfit)
mods_su_dAICc


# Survival plots by year ------------------------------------------------------
make_su_year_plot <- function(year_i) {
  df_i <- df_su %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tidyr::expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_su$year)))

  pred_i <- pred_i %>%
    mutate(
      survives = predict(
        mod_su_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, disturbance) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, survives, color = factor(disturbance)),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, survives, color = factor(disturbance),
          linetype = factor(disturbance)), linewidth = 0.7) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = expression('log(diameter)'[t0]),
      y = 'Survival probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

su_yrs <- lapply(levels(df_su$year), make_su_year_plot)
fig_su_years <- wrap_plots(su_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Survival - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_su_years


# Growth data -----------------------------------------------------------------
df_gr <- df %>%
  filter(
    state == 'active',
    state_t1 == 'active',
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    is.finite(logsize_t0_2),
    is.finite(logsize_t0_3),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    id, year, size_t0, size_t1,
    logsize_t0, logsize_t1,
    logsize_t0_2, logsize_t0_3, disturbance)

fig_gr_raw <- ggplot(df_gr, aes(logsize_t0, logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  labs(
    title = 'Growth',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = expression('log(diameter)'[t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_raw


# Growth model ----------------------------------------------------------------
mod_gr_00 <- lmer(
  logsize_t1 ~ 1 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_0 <- lmer(
  logsize_t1 ~ disturbance + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_10 <- lmer(
  logsize_t1 ~ logsize_t0 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + disturbance + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_20 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 +
    (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance +
    (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_30 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    disturbance + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mods_gr <- list(
  mod_gr_00, mod_gr_0, mod_gr_10, mod_gr_1,
  mod_gr_20, mod_gr_2, mod_gr_30, mod_gr_3)

mods_gr_dAICc <- bbmle::AICctab(
  mods_gr, weights = TRUE, sort = FALSE)$dAICc

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- which.min(mods_gr_dAICc)
} else {
  keep <- 2 * v_mod_set_gr + c(1, 2)
  mod_gr_index_bestfit <- keep[which.min(mods_gr_dAICc[keep])]
}

v_mod_gr_index <- floor((mod_gr_index_bestfit - 1) / 2)
mod_gr_bestfit <- mods_gr[[mod_gr_index_bestfit]]

mod_gr_bestfit
summary(mod_gr_bestfit)
mods_gr_dAICc


# Growth plots by year --------------------------------------------------------
make_gr_year_plot <- function(year_i) {
  df_i <- df_gr %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tidyr::expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_gr$year)))

  pred_i <- pred_i %>%
    mutate(
      logsize_t1 = predict(
        mod_gr_bestfit, newdata = pred_i, re.form = NULL,
        allow.new.levels = TRUE))

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, disturbance) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      logsize_t1 = mean(logsize_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, logsize_t1, color = factor(disturbance)),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, logsize_t1, color = factor(disturbance),
          linetype = factor(disturbance)), linewidth = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', lty = 2) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = expression('log(diameter)'[t0]),
      y = expression('log(diameter)'[t1])) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

gr_yrs <- lapply(levels(df_gr$year), make_gr_year_plot)
fig_gr_years <- wrap_plots(gr_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Growth - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_gr_years


# Growth variance -------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_bestfit)
mod_gr_y <- resid(mod_gr_bestfit)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Dormancy entry data ---------------------------------------------------------
df_do <- df %>%
  filter(
    !is.na(enter_dormancy),
    size_t0 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t0_2),
    is.finite(logsize_t0_3),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance)) %>%
  dplyr::select(
    id, year, size_t0, enter_dormancy, disturbance,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Dormancy entry model --------------------------------------------------------
mod_do_00 <- glmer(
  enter_dormancy ~ 1 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_0 <- glmer(
  enter_dormancy ~ disturbance + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_10 <- glmer(
  enter_dormancy ~ logsize_t0 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_1 <- glmer(
  enter_dormancy ~ logsize_t0 + disturbance + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_20 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_2 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + disturbance +
    (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_30 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_3 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    disturbance + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mods_do <- list(
  mod_do_00, mod_do_0, mod_do_10, mod_do_1,
  mod_do_20, mod_do_2, mod_do_30, mod_do_3)

mods_do_dAICc <- bbmle::AICctab(
  mods_do, weights = TRUE, sort = FALSE)$dAICc

if (length(v_mod_set_do) == 0) {
  mod_do_index_bestfit <- which.min(mods_do_dAICc)
} else {
  keep <- 2 * v_mod_set_do + c(1, 2)
  mod_do_index_bestfit <- keep[which.min(mods_do_dAICc[keep])]
}

v_mod_do_index <- floor((mod_do_index_bestfit - 1) / 2)
mod_do_bestfit <- mods_do[[mod_do_index_bestfit]]

mod_do_bestfit
summary(mod_do_bestfit)
mods_do_dAICc


# Dormancy plots by year ------------------------------------------------------
make_do_year_plot <- function(year_i) {
  df_i <- df_do %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tidyr::expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_do$year)))

  pred_i <- pred_i %>%
    mutate(
      enter_dormancy = predict(
        mod_do_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, disturbance) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      enter_dormancy = mean(enter_dormancy, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, enter_dormancy, color = factor(disturbance)),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, enter_dormancy, color = factor(disturbance),
          linetype = factor(disturbance)), linewidth = 0.7) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = expression('log(diameter)'[t0]),
      y = 'Dormancy probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

do_yrs <- lapply(levels(df_do$year), make_do_year_plot)
fig_do_years <- wrap_plots(do_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Dormancy entry - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_do_years


# Reactivation data -----------------------------------------------------------
df_ra <- df %>%
  filter(
    !is.na(reactivate),
    !is.na(disturbance)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance))


# Reactivation model ----------------------------------------------------------
mod_ra_0 <- glmer(
  reactivate ~ 1 + (1 | year),
  data = df_ra, family = binomial, control = ctrl_glmer)

mod_ra_1 <- glmer(
  reactivate ~ disturbance + (1 | year),
  data = df_ra, family = binomial, control = ctrl_glmer)

mods_ra <- list(mod_ra_0, mod_ra_1)

mods_ra_dAICc <- bbmle::AICctab(
  mods_ra, weights = TRUE, sort = FALSE)$dAICc

mod_ra_index_bestfit <- which.min(mods_ra_dAICc)
mod_ra_bestfit <- mods_ra[[mod_ra_index_bestfit]]

mod_ra_bestfit
summary(mod_ra_bestfit)
mods_ra_dAICc


# Reactivation plots by year --------------------------------------------------
make_ra_year_plot <- function(year_i) {
  df_i <- df_ra %>%
    filter(year == year_i)

  pred_i <- tibble(
    disturbance = c(0, 1),
    year = factor(year_i, levels = levels(df_ra$year)))

  pred_i <- pred_i %>%
    mutate(
      reactivate = predict(
        mod_ra_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  pts_i <- df_i %>%
    group_by(disturbance) %>%
    summarise(
      reactivate = mean(reactivate, na.rm = TRUE),
      n = n(), .groups = 'drop')

  ggplot() +
    geom_point(
      data = pts_i,
      aes(factor(disturbance), reactivate,
          color = factor(disturbance)), size = 1.8) +
    geom_point(
      data = pred_i,
      aes(factor(disturbance), reactivate,
          color = factor(disturbance)), shape = 95, size = 6) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = 'Current-year fire',
      y = 'Reactivation probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

ra_yrs <- lapply(levels(df_ra$year), make_ra_year_plot)
fig_ra_years <- wrap_plots(ra_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Reactivation - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_ra_years


# Reactivation size distribution ---------------------------------------------
df_ra_size <- df %>%
  filter(
    reactivate == 1,
    size_reactivate_t1 > 0) %>%
  mutate(logsize_reactivate = log(size_reactivate_t1))

react_sz <- mean(
  df_ra_size$logsize_reactivate,
  na.rm = TRUE)

react_sd <- sd(
  df_ra_size$logsize_reactivate,
  na.rm = TRUE)


# Flower data -----------------------------------------------------------------
df_fl <- df %>%
  filter(
    state == 'active',
    !is.na(flower),
    size_t0 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t0_2),
    is.finite(logsize_t0_3),
    !is.na(disturbance_prev)) %>%
  mutate(
    year = factor(year),
    disturbance_prev = as.numeric(disturbance_prev))


# Flower model ----------------------------------------------------------------
mod_fl_00 <- glmer(
  flower ~ 1 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_0 <- glmer(
  flower ~ disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_10 <- glmer(
  flower ~ logsize_t0 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_1 <- glmer(
  flower ~ logsize_t0 + disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_20 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + disturbance_prev +
    (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_30 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    disturbance_prev + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mods_fl <- list(
  mod_fl_00, mod_fl_0, mod_fl_10, mod_fl_1,
  mod_fl_20, mod_fl_2, mod_fl_30, mod_fl_3)

mods_fl_dAICc <- bbmle::AICctab(
  mods_fl, weights = TRUE, sort = FALSE)$dAICc

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- which.min(mods_fl_dAICc)
} else {
  keep <- 2 * v_mod_set_fl + c(1, 2)
  mod_fl_index_bestfit <- keep[which.min(mods_fl_dAICc[keep])]
}

v_mod_fl_index <- floor((mod_fl_index_bestfit - 1) / 2)
mod_fl_bestfit <- mods_fl[[mod_fl_index_bestfit]]

mod_fl_bestfit
summary(mod_fl_bestfit)
mods_fl_dAICc


# Flowering plots by year -----------------------------------------------------
make_fl_year_plot <- function(year_i) {
  df_i <- df_fl %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tidyr::expand_grid(
    logsize_t0 = x,
    disturbance_prev = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fl$year)))

  pred_i <- pred_i %>%
    mutate(
      flower = predict(
        mod_fl_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, disturbance_prev) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      flower = mean(flower, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, flower, color = factor(disturbance_prev)),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, flower, color = factor(disturbance_prev),
          linetype = factor(disturbance_prev)), linewidth = 0.7) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = expression('log(diameter)'[t0]),
      y = 'Flowering probability') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

fl_yrs <- lapply(levels(df_fl$year), make_fl_year_plot)
fig_fl_years <- wrap_plots(fl_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Flowering - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_fl_years


# Number of scapes conditional on flowering data -----------------------------
df_fl_cond <- df_fl %>%
  filter(
    flower == 1,
    !is.na(fl_nr),
    fl_nr > 0,
    fl_nr == round(fl_nr),
    !is.na(disturbance_prev)) %>%
  mutate(year = factor(year))


# Number of scapes model ------------------------------------------------------
mod_fl_n_00 <- glmer.nb(
  fl_nr ~ 1 + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_0 <- glmer.nb(
  fl_nr ~ disturbance_prev + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_10 <- glmer.nb(
  fl_nr ~ logsize_t0 + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_1 <- glmer.nb(
  fl_nr ~ logsize_t0 + disturbance_prev + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_20 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_2 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + disturbance_prev +
    (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_30 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mod_fl_n_3 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    disturbance_prev + (1 | year),
  data = df_fl_cond, control = ctrl_glmer)

mods_fl_n <- list(
  mod_fl_n_00, mod_fl_n_0, mod_fl_n_10, mod_fl_n_1,
  mod_fl_n_20, mod_fl_n_2, mod_fl_n_30, mod_fl_n_3)

mods_fl_n_dAICc <- bbmle::AICctab(
  mods_fl_n, weights = TRUE, sort = FALSE)$dAICc

if (length(v_mod_set_fl_n) == 0) {
  mod_fl_n_index_bestfit <- which.min(mods_fl_n_dAICc)
} else {
  keep <- 2 * v_mod_set_fl_n + c(1, 2)
  mod_fl_n_index_bestfit <- keep[which.min(mods_fl_n_dAICc[keep])]
}

v_mod_fl_n_index <- floor((mod_fl_n_index_bestfit - 1) / 2)
mod_fl_n_bestfit <- mods_fl_n[[mod_fl_n_index_bestfit]]

mod_fl_n_bestfit
summary(mod_fl_n_bestfit)
mods_fl_n_dAICc


# Scape-number plots by year --------------------------------------------------
make_fl_n_year_plot <- function(year_i) {
  df_i <- df_fl_cond %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tidyr::expand_grid(
    logsize_t0 = x,
    disturbance_prev = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fl_cond$year)))

  pred_i <- pred_i %>%
    mutate(
      fl_nr = predict(
        mod_fl_n_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, disturbance_prev) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      fl_nr = mean(fl_nr, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, fl_nr, color = factor(disturbance_prev)),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, fl_nr, color = factor(disturbance_prev),
          linetype = factor(disturbance_prev)), linewidth = 0.7) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
      name = 'Fire') +
    labs(
      title = year_i,
      x = expression('log(diameter)'[t0]),
      y = 'Number of flowering scapes') +
    theme_bw() +
    theme(text = element_text(size = 5))
}

fl_n_yrs <- lapply(levels(df_fl_cond$year), make_fl_n_year_plot)
fig_fl_n_years <- wrap_plots(fl_n_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Scape number - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_fl_n_years


# Recruitment data ------------------------------------------------------------
df_sc2re <- df_re %>%
  filter(
    !is.na(nr_scapes),
    !is.na(recruits_simple)) %>%
  group_by(site, pop, year) %>%
  summarise(
    total_scapes = sum(nr_scapes),
    recruit_count = sum(recruits_simple),
    disturbance = case_when(
      any(disturbance == 1, na.rm = TRUE) ~ 1,
      any(is.na(disturbance)) ~ NA_real_,
      TRUE ~ 0),
    disturbance_prev = case_when(
      any(disturbance_prev == 1, na.rm = TRUE) ~ 1,
      any(is.na(disturbance_prev)) ~ NA_real_,
      TRUE ~ 0),
    .groups = 'drop')

df_sc2re_mod <- df_sc2re %>%
  filter(
    !is.na(total_scapes),
    !is.na(recruit_count),
    !is.na(disturbance),
    !(year %in% v_years_re)) %>%
  mutate(
    year = factor(year),
    disturbance = as.numeric(disturbance),
    disturbance_prev = as.numeric(disturbance_prev))


# Recruitment model -----------------------------------------------------------
mod_re_00 <- glmer.nb(
  recruit_count ~ 1 + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_0 <- glmer.nb(
  recruit_count ~ disturbance + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_sc0 <- glmer.nb(
  recruit_count ~ total_scapes + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_sc1 <- glmer.nb(
  recruit_count ~ total_scapes + disturbance + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_sc2 <- glmer.nb(
  recruit_count ~ total_scapes * disturbance + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_logsc0 <- glmer.nb(
  recruit_count ~ log1p(total_scapes) + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_logsc1 <- glmer.nb(
  recruit_count ~ log1p(total_scapes) + disturbance + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

mod_re_logsc2 <- glmer.nb(
  recruit_count ~ log1p(total_scapes) * disturbance + (1 | year),
  data = df_sc2re_mod, control = ctrl_re)

# Current fire effect varies among recruitment years, as in the POLE yearly IPM.
# Previous fire is not added here; it already acts through flowering and scapes.
mod_re_yrfire <- glmer.nb(
  recruit_count ~ disturbance + (1 + disturbance || year),
  data = df_sc2re_mod, control = ctrl_re)

mods_re <- list(
  mod_re_00,
  mod_re_0,
  mod_re_sc0,
  mod_re_sc1,
  mod_re_sc2,
  mod_re_logsc0,
  mod_re_logsc1,
  mod_re_logsc2,
  mod_re_yrfire)

mods_re_dAICc <- bbmle::AICctab(
  mods_re, weights = TRUE, sort = FALSE)$dAICc

mods_re_sorted <- order(mods_re_dAICc)

if (length(v_mod_set_re) == 0) {
  mod_re_index_bestfit <- mods_re_sorted[1]
} else {
  mod_re_index_bestfit <- v_mod_set_re + 1
}

v_mod_re_index <- mod_re_index_bestfit - 1
mod_re_bestfit <- mods_re[[mod_re_index_bestfit]]

mod_re_bestfit
summary(mod_re_bestfit)
mods_re_dAICc


# Recruitment plots by year --------------------------------------------------
make_re_year_plot <- function(year_i) {
  df_i <- df_sc2re_mod %>%
    filter(year == year_i)

  x <- seq(
    0,
    max(df_i$total_scapes, na.rm = TRUE),
    length.out = 100)

  pred_i <- tidyr::expand_grid(
    total_scapes = x,
    disturbance = c(0, 1)) %>%
    mutate(
      year = factor(year_i, levels = levels(df_sc2re_mod$year)))

  pred_i <- pred_i %>%
    mutate(
      recruit_count = predict(
        mod_re_bestfit, newdata = pred_i, type = 'response',
        re.form = NULL, allow.new.levels = TRUE))

  ggplot() +
    geom_point(
      data = df_i,
      aes(total_scapes, recruit_count, color = factor(disturbance)),
      alpha = 0.75, size = 1.5) +
    geom_line(
      data = pred_i,
      aes(total_scapes, recruit_count, color = factor(disturbance),
          linetype = factor(disturbance)), linewidth = 0.8) +
    scale_color_manual(
      values = c('0' = 'black', '1' = 'red'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    scale_linetype_manual(
      values = c('0' = 'solid', '1' = 'dashed'),
      labels = c('0' = 'No fire', '1' = 'Current-year fire'),
      name = 'Fire') +
    theme_bw() +
    labs(
      title = year_i,
      x = expression('Total flowering scapes '[t0]),
      y = expression('Number of recruits '[t1])) +
    theme(text = element_text(size = 5))
}

re_yrs <- lapply(levels(df_sc2re_mod$year), make_re_year_plot)
fig_re_years <- wrap_plots(re_yrs, guides = 'collect') +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Recruitment - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9))) &
  theme(legend.position = 'bottom')

fig_re_years


# Mean/fixed-effect vital-rate diagnostic plots -------------------------------
# These retain the diagnostic plots from the ERLO mean IPM. Because the models
# are now mixed, lines use fixed effects only (re.form = NA); yearly plots above
# include the corresponding year random effects.

# Mean survival plot ----------------------------------------------------------
df_su_mean_pred <- tidyr::expand_grid(
  logsize_t0 = seq(
    min(df_su$logsize_t0, na.rm = TRUE),
    max(df_su$logsize_t0, na.rm = TRUE),
    length.out = 100),
  disturbance = c(0, 1)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    year = factor(levels(df_su$year)[1], levels = levels(df_su$year)))

df_su_mean_pred <- df_su_mean_pred %>%
  mutate(
    survives = predict(
      mod_su_bestfit, newdata = df_su_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

df_su_mean_binned <- df_su %>%
  mutate(bin = cut(logsize_t0, breaks = 10)) %>%
  group_by(bin, disturbance) %>%
  summarise(
    logsize_t0 = mean(logsize_t0, na.rm = TRUE),
    survives = mean(survives, na.rm = TRUE),
    n = n(), .groups = 'drop') %>%
  filter(!is.na(logsize_t0), n > 0)

fig_su <- ggplot() +
  geom_point(
    data = df_su_mean_binned,
    aes(logsize_t0, survives, color = factor(disturbance))) +
  geom_line(
    data = df_su_mean_pred,
    aes(logsize_t0, survives, color = factor(disturbance),
        linetype = factor(disturbance)), linewidth = 1) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    name = NULL) +
  scale_linetype_manual(
    values = c('0' = 'solid', '1' = 'dashed'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    name = NULL) +
  labs(
    title = 'Survival probability by size',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = 'Probability of survival') +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

fig_su


# Mean growth plot ------------------------------------------------------------
df_gr_mean_pred <- tidyr::expand_grid(
  logsize_t0 = seq(
    min(df_gr$logsize_t0, na.rm = TRUE),
    max(df_gr$logsize_t0, na.rm = TRUE),
    length.out = 100),
  disturbance = c(0, 1)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    year = factor(levels(df_gr$year)[1], levels = levels(df_gr$year)))

df_gr_mean_pred <- df_gr_mean_pred %>%
  mutate(
    logsize_t1 = predict(
      mod_gr_bestfit, newdata = df_gr_mean_pred, re.form = NA,
      allow.new.levels = TRUE))

fig_gr <- ggplot(
  df_gr,
  aes(logsize_t0, logsize_t1, color = factor(disturbance))) +
  geom_point(alpha = 0.25) +
  geom_line(
    data = df_gr_mean_pred,
    aes(logsize_t0, logsize_t1, color = factor(disturbance)),
    linewidth = 1) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    name = NULL) +
  labs(
    title = 'Growth prediction',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = expression('log(diameter)'[t1])) +
  theme_bw()

fig_gr


# Mean dormancy plot ----------------------------------------------------------
df_do_mean_pred <- tidyr::expand_grid(
  logsize_t0 = seq(
    min(df_do$logsize_t0, na.rm = TRUE),
    max(df_do$logsize_t0, na.rm = TRUE),
    length.out = 100),
  disturbance = c(0, 1)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    year = factor(levels(df_do$year)[1], levels = levels(df_do$year)))

df_do_mean_pred <- df_do_mean_pred %>%
  mutate(
    enter_dormancy = predict(
      mod_do_bestfit, newdata = df_do_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

fig_do <- ggplot(
  df_do,
  aes(logsize_t0, enter_dormancy, color = factor(disturbance))) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.15) +
  geom_line(
    data = df_do_mean_pred,
    aes(logsize_t0, enter_dormancy, color = factor(disturbance)),
    linewidth = 1) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    name = NULL) +
  labs(
    title = 'Dormancy probability by size',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = 'Probability of entering dormancy') +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

fig_do


# Mean reactivation plot ------------------------------------------------------
df_ra_mean_plot <- df_ra %>%
  group_by(disturbance) %>%
  summarise(
    n = n(),
    nr_reactivate = sum(reactivate == 1),
    reactivation = mean(reactivate), .groups = 'drop') %>%
  mutate(
    lwr = binom.confint(
      nr_reactivate, n, methods = 'wilson')$lower,
    upr = binom.confint(
      nr_reactivate, n, methods = 'wilson')$upper)

df_ra_mean_pred <- tibble(
  disturbance = c(0, 1),
  year = factor(levels(df_ra$year)[1], levels = levels(df_ra$year)))

df_ra_mean_pred <- df_ra_mean_pred %>%
  mutate(
    reactivation = predict(
      mod_ra_bestfit, newdata = df_ra_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

fig_ra <- ggplot(
  df_ra_mean_plot,
  aes(factor(disturbance), reactivation,
      color = factor(disturbance))) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
  geom_point(
    data = df_ra_mean_pred,
    aes(factor(disturbance), reactivation,
        color = factor(disturbance)), shape = 95, size = 10) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    guide = 'none') +
  labs(
    title = 'Reactivation probability',
    subtitle = v_ggp_suffix,
    x = '',
    y = 'Probability of reactivation') +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

fig_ra


# Mean flowering plot ---------------------------------------------------------
df_fl_mean_pred <- tidyr::expand_grid(
  logsize_t0 = seq(
    min(df_fl$logsize_t0, na.rm = TRUE),
    max(df_fl$logsize_t0, na.rm = TRUE),
    length.out = 100),
  disturbance_prev = c(0, 1)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    year = factor(levels(df_fl$year)[1], levels = levels(df_fl$year)))

df_fl_mean_pred <- df_fl_mean_pred %>%
  mutate(
    flower = predict(
      mod_fl_bestfit, newdata = df_fl_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

fig_fl <- ggplot(
  df_fl,
  aes(logsize_t0, flower, color = factor(disturbance_prev))) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.15) +
  geom_line(
    data = df_fl_mean_pred,
    aes(logsize_t0, flower, color = factor(disturbance_prev)),
    linewidth = 1) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
    name = NULL) +
  labs(
    title = 'Flowering probability by size',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = 'Probability of flowering') +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

fig_fl


# Mean scape-number plot ------------------------------------------------------
df_fl_n_mean_pred <- tidyr::expand_grid(
  logsize_t0 = seq(
    min(df_fl_cond$logsize_t0, na.rm = TRUE),
    max(df_fl_cond$logsize_t0, na.rm = TRUE),
    length.out = 100),
  disturbance_prev = c(0, 1)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    year = factor(
      levels(df_fl_cond$year)[1],
      levels = levels(df_fl_cond$year)))

df_fl_n_mean_pred <- df_fl_n_mean_pred %>%
  mutate(
    fl_nr = predict(
      mod_fl_n_bestfit, newdata = df_fl_n_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

fig_fl_n <- ggplot(
  df_fl_cond,
  aes(logsize_t0, fl_nr, color = factor(disturbance_prev))) +
  geom_jitter(alpha = 0.2, width = 0.08, height = 0.3) +
  geom_line(
    data = df_fl_n_mean_pred,
    aes(logsize_t0, fl_nr, color = factor(disturbance_prev)),
    linewidth = 1) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Previous-year fire'),
    name = NULL) +
  labs(
    title = 'Flower number',
    subtitle = v_ggp_suffix,
    x = expression('log(diameter)'[t0]),
    y = 'Number of flowering scapes') +
  theme_bw()

fig_fl_n


# Mean scapes-to-recruits plot -----------------------------------------------
df_sc2re_mean_pred <- tidyr::expand_grid(
  total_scapes = seq(
    0,
    max(df_sc2re_mod$total_scapes, na.rm = TRUE),
    length.out = 100),
  disturbance = c(0, 1)) %>%
  mutate(
    year = factor(
      levels(df_sc2re_mod$year)[1],
      levels = levels(df_sc2re_mod$year)))

df_sc2re_mean_pred <- df_sc2re_mean_pred %>%
  mutate(
    recruit_count = predict(
      mod_re_bestfit, newdata = df_sc2re_mean_pred, type = 'response',
      re.form = NA, allow.new.levels = TRUE))

fig_sc2re <- ggplot(
  df_sc2re_mod,
  aes(total_scapes, recruit_count, color = factor(disturbance))) +
  geom_jitter(height = 0.2, width = 0.5, alpha = 0.4) +
  geom_line(
    data = df_sc2re_mean_pred,
    aes(total_scapes, recruit_count, color = factor(disturbance)),
    linewidth = 1.2) +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    labels = c('0' = 'No fire', '1' = 'Current-year fire'),
    name = NULL) +
  labs(
    title = 'Recruits t1 by flowering scapes t0',
    subtitle = v_ggp_suffix,
    x = expression('Total flowering scapes '[t0]),
    y = expression('Number of recruits '[t1])) +
  theme_bw()

fig_sc2re


# Recruit size distribution ---------------------------------------------------
df_re_size <- df_ind %>%
  filter(
    recruit_type == 'seedling',
    size_t0 > 0,
    is.finite(logsize_t0))

recr_sz <- mean(
  df_re_size$logsize_t0,
  na.rm = TRUE)

recr_sd <- sd(
  df_re_size$logsize_t0,
  na.rm = TRUE)


# Recruitment conversion helpers ---------------------------------------------
# The ERLO data/model describe recruits from total scapes at the site-pop-year
# level. The IPM uses the same observed scape distribution to obtain a recruits
# per scape conversion, as in the mean ERLO IPM. For a yearly IPM, the year
# random effect is included and that year's observed scape distribution is used.
get_recr_per_scape <- function(disturbance_i, year_i = NULL) {
  if (is.null(year_i)) {
    df_i <- df_sc2re_mod %>%
      mutate(disturbance = disturbance_i)

    pred_i <- predict(
      mod_re_bestfit,
      newdata = df_i,
      type = 'response',
      re.form = NA,
      allow.new.levels = TRUE)
  } else {
    df_i <- df_sc2re_mod %>%
      filter(as.character(year) == as.character(year_i)) %>%
      mutate(disturbance = disturbance_i)

    if (nrow(df_i) == 0 ||
        sum(df_i$total_scapes, na.rm = TRUE) <= 0) {
      return(get_recr_per_scape(disturbance_i))
    }

    pred_i <- predict(
      mod_re_bestfit,
      newdata = df_i,
      type = 'response',
      re.form = NULL,
      allow.new.levels = TRUE)
  }

  total_scapes_i <- sum(df_i$total_scapes, na.rm = TRUE)

  if (!is.finite(total_scapes_i) || total_scapes_i <= 0) {
    return(0)
  }

  sum(pred_i, na.rm = TRUE) / total_scapes_i
}

recr_per_scape_0 <- get_recr_per_scape(0)
recr_per_scape_1 <- get_recr_per_scape(1)

recr_per_scape_0
recr_per_scape_1


# Parameter extraction helpers -----------------------------------------------
get_fixef_safe <- function(model) {
  if (inherits(model, 'merMod')) {
    fixef(model)
  } else {
    coef(model)
  }
}

get_year_coef <- function(model, year_i = NULL, group_var = 'year') {
  fixed <- get_fixef_safe(model)

  if (is.null(year_i) || !inherits(model, 'merMod')) {
    return(fixed)
  }

  coef_list <- coef(model)

  if (!(group_var %in% names(coef_list))) {
    return(fixed)
  }

  coef_matrix <- coef_list[[group_var]]
  year_i <- as.character(year_i)

  if (!(year_i %in% rownames(coef_matrix))) {
    return(fixed)
  }

  out <- unlist(coef_matrix[year_i, , drop = TRUE])
  names(out) <- colnames(coef_matrix)
  out
}

get_coef <- function(x, term, default = 0) {
  if (term %in% names(x)) {
    unname(x[term])
  } else {
    default
  }
}


# IPM constants ---------------------------------------------------------------
mesh_limits <- range(
  c(
    df_gr$logsize_t0,
    df_gr$logsize_t1,
    df_re_size$logsize_t0,
    df_ra_size$logsize_reactivate),
  na.rm = TRUE,
  finite = TRUE)

gr_var_coef <- coef(mod_gr_var)


# Build mean or year-specific IPM parameters ---------------------------------
make_ipm_pars <- function(year_i = NULL) {
  su_cf <- get_year_coef(mod_su_bestfit, year_i)
  gr_cf <- get_year_coef(mod_gr_bestfit, year_i)
  do_cf <- get_year_coef(mod_do_bestfit, year_i)
  ra_cf <- get_year_coef(mod_ra_bestfit, year_i)
  fl_cf <- get_year_coef(mod_fl_bestfit, year_i)
  fln_cf <- get_year_coef(mod_fl_n_bestfit, year_i)

  recr_0_i <- get_recr_per_scape(0, year_i)
  recr_1_i <- get_recr_per_scape(1, year_i)

  list(
    prefix = v_script_prefix,
    species = v_species,

    surv_b0 = get_coef(su_cf, '(Intercept)'),
    surv_b1 = get_coef(su_cf, 'logsize_t0'),
    surv_b2 = get_coef(su_cf, 'logsize_t0_2'),
    surv_b3 = get_coef(su_cf, 'logsize_t0_3'),
    surv_bd = get_coef(su_cf, 'disturbance'),

    grow_b0 = get_coef(gr_cf, '(Intercept)'),
    grow_b1 = get_coef(gr_cf, 'logsize_t0'),
    grow_b2 = get_coef(gr_cf, 'logsize_t0_2'),
    grow_b3 = get_coef(gr_cf, 'logsize_t0_3'),
    grow_bd = get_coef(gr_cf, 'disturbance'),
    a = as.numeric(gr_var_coef['a']),
    b = as.numeric(gr_var_coef['b']),

    dorm_b0 = get_coef(do_cf, '(Intercept)'),
    dorm_b1 = get_coef(do_cf, 'logsize_t0'),
    dorm_b2 = get_coef(do_cf, 'logsize_t0_2'),
    dorm_b3 = get_coef(do_cf, 'logsize_t0_3'),
    dorm_bd = get_coef(do_cf, 'disturbance'),

    react_b0 = get_coef(ra_cf, '(Intercept)'),
    react_bd = get_coef(ra_cf, 'disturbance'),

    fl_b0 = get_coef(fl_cf, '(Intercept)'),
    fl_b1 = get_coef(fl_cf, 'logsize_t0'),
    fl_b2 = get_coef(fl_cf, 'logsize_t0_2'),
    fl_b3 = get_coef(fl_cf, 'logsize_t0_3'),
    fl_bd = get_coef(fl_cf, 'disturbance_prev'),

    fln_b0 = get_coef(fln_cf, '(Intercept)'),
    fln_b1 = get_coef(fln_cf, 'logsize_t0'),
    fln_b2 = get_coef(fln_cf, 'logsize_t0_2'),
    fln_b3 = get_coef(fln_cf, 'logsize_t0_3'),
    fln_bd = get_coef(fln_cf, 'disturbance_prev'),

    recr_0 = recr_0_i,
    recr_1 = recr_1_i,
    recr_sz = recr_sz,
    recr_sd = recr_sd,
    react_sz = react_sz,
    react_sd = react_sd,

    L = mesh_limits[1],
    U = mesh_limits[2],
    mat_siz = 200,

    mod_su_index = v_mod_su_index,
    mod_gr_index = v_mod_gr_index,
    mod_do_index = v_mod_do_index,
    mod_fl_index = v_mod_fl_index,
    mod_fl_n_index = v_mod_fl_n_index,
    mod_re_index = v_mod_re_index)
}

pars_mean <- make_ipm_pars()
pars <- pars_mean


# Disturbance transition probabilities ---------------------------------------
# Current disturbance refers to t -> t + 1 and previous disturbance to
# t - 1 -> t. Recruitment rows provide the monitored quadrat-year inventory.
df_disturbance_regime <- df_re_all %>%
  filter(
    !is.na(site),
    !is.na(pop),
    !is.na(qu),
    !is.na(year),
    !is.na(disturbance),
    !is.na(disturbance_prev)) %>%
  distinct(
    site, pop, qu, year,
    disturbance, disturbance_prev) %>%
  count(
    disturbance,
    disturbance_prev,
    name = 'n_quadrat_years') %>%
  complete(
    disturbance = c(0, 1),
    disturbance_prev = c(0, 1),
    fill = list(n_quadrat_years = 0)) %>%
  arrange(
    disturbance_prev,
    disturbance)

df_disturbance_regime


df_disturbance_transition <- df_disturbance_regime %>%
  group_by(disturbance_prev) %>%
  mutate(
    p_transition = n_quadrat_years /
      sum(n_quadrat_years)) %>%
  ungroup()

df_disturbance_transition


# Site-specific disturbance transition probabilities -------------------------
df_disturbance_regime_site <- df_re_all %>%
  filter(
    !is.na(site),
    !is.na(pop),
    !is.na(qu),
    !is.na(year),
    !is.na(disturbance),
    !is.na(disturbance_prev)) %>%
  distinct(
    site, pop, qu, year,
    disturbance, disturbance_prev) %>%
  count(
    site, disturbance, disturbance_prev,
    name = 'n_quadrat_years') %>%
  complete(
    site,
    disturbance = c(0, 1),
    disturbance_prev = c(0, 1),
    fill = list(n_quadrat_years = 0)) %>%
  arrange(site, disturbance_prev, disturbance)

df_disturbance_transition_site <- df_disturbance_regime_site %>%
  group_by(site, disturbance_prev) %>%
  mutate(
    n_from_state = sum(n_quadrat_years),
    p_transition = if_else(
      n_from_state > 0,
      n_quadrat_years / n_from_state,
      NA_real_)) %>%
  ungroup()

df_transition_site <- df_disturbance_transition_site %>%
  mutate(
    transition = case_when(
      disturbance == 0 & disturbance_prev == 0 ~ 'p00',
      disturbance == 1 & disturbance_prev == 0 ~ 'p10',
      disturbance == 0 & disturbance_prev == 1 ~ 'p01',
      disturbance == 1 & disturbance_prev == 1 ~ 'p11')) %>%
  select(site, transition, p_transition) %>%
  pivot_wider(names_from = transition, values_from = p_transition)

df_disturbance_transition_site
df_transition_site

df_transition_site %>%
  filter(if_any(c(p00, p10, p01, p11), is.na))

get_transition_probability <- function(
    disturbance_i,
    disturbance_prev_i) {

  df_disturbance_transition %>%
    filter(
      disturbance == disturbance_i,
      disturbance_prev == disturbance_prev_i) %>%
    pull(p_transition)
}

p_0_given_0 <- get_transition_probability(0, 0)
p_1_given_0 <- get_transition_probability(1, 0)
p_0_given_1 <- get_transition_probability(0, 1)
p_1_given_1 <- get_transition_probability(1, 1)


# Previous-disturbance probabilities -----------------------------------------
df_previous_disturbance <- df_disturbance_regime %>%
  group_by(disturbance_prev) %>%
  summarise(
    n_quadrat_years = sum(n_quadrat_years),
    .groups = 'drop') %>%
  mutate(
    p_previous = n_quadrat_years /
      sum(n_quadrat_years))

df_previous_disturbance

p_prev_0 <- df_previous_disturbance %>%
  filter(disturbance_prev == 0) %>%
  pull(p_previous)

p_prev_1 <- df_previous_disturbance %>%
  filter(disturbance_prev == 1) %>%
  pull(p_previous)


# IPM helper functions --------------------------------------------------------
inv_logit <- function(x) {
  plogis(x)
}


# Survival --------------------------------------------------------------------
sx <- function(x, pars, disturbance = 0) {
  inv_logit(
    pars$surv_b0 +
      pars$surv_b1 * x +
      pars$surv_b2 * x^2 +
      pars$surv_b3 * x^3 +
      pars$surv_bd * disturbance)
}


# Growth ----------------------------------------------------------------------
grow_mu <- function(x, pars, disturbance = 0) {
  pars$grow_b0 +
    pars$grow_b1 * x +
    pars$grow_b2 * x^2 +
    pars$grow_b3 * x^3 +
    pars$grow_bd * disturbance
}


# Growth variation ------------------------------------------------------------
grow_sd <- function(mu, pars) {
  sqrt(pars$a * exp(pars$b * mu))
}


# Growth transition -----------------------------------------------------------
gxy <- function(x, y, pars, disturbance = 0) {
  mu <- grow_mu(
    x, pars,
    disturbance = disturbance)

  dnorm(
    y,
    mean = mu,
    sd = grow_sd(mu, pars))
}


# Dormancy entry --------------------------------------------------------------
dorm_x <- function(x, pars, disturbance = 0) {
  inv_logit(
    pars$dorm_b0 +
      pars$dorm_b1 * x +
      pars$dorm_b2 * x^2 +
      pars$dorm_b3 * x^3 +
      pars$dorm_bd * disturbance)
}


# Reactivation ----------------------------------------------------------------
react_p <- function(pars, disturbance = 0) {
  inv_logit(
    pars$react_b0 +
      pars$react_bd * disturbance)
}


# Flowering probability -------------------------------------------------------
fl_x <- function(x, pars, disturbance_prev = 0) {
  inv_logit(
    pars$fl_b0 +
      pars$fl_b1 * x +
      pars$fl_b2 * x^2 +
      pars$fl_b3 * x^3 +
      pars$fl_bd * disturbance_prev)
}


# Scape number conditional on flowering --------------------------------------
fl_n_x <- function(x, pars, disturbance_prev = 0) {
  exp(
    pars$fln_b0 +
      pars$fln_b1 * x +
      pars$fln_b2 * x^2 +
      pars$fln_b3 * x^3 +
      pars$fln_bd * disturbance_prev)
}


# Recruitment conversion ------------------------------------------------------
recruits_per_scape <- function(pars, disturbance = 0) {
  if (disturbance == 0) {
    pars$recr_0
  } else {
    pars$recr_1
  }
}


# Recruitment size distribution ---------------------------------------------
re_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(
    y,
    mean = pars$recr_sz,
    sd = pars$recr_sd)

  if (!is.null(h)) {
    dens <- dens / sum(dens * h)
  } else {
    norm <- pnorm(
      pars$U, mean = pars$recr_sz, sd = pars$recr_sd) -
      pnorm(
        pars$L, mean = pars$recr_sz, sd = pars$recr_sd)
    dens <- dens / norm
  }

  dens
}


# Reactivation size distribution ---------------------------------------------
react_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(
    y,
    mean = pars$react_sz,
    sd = pars$react_sd)

  if (!is.null(h)) {
    dens <- dens / sum(dens * h)
  } else {
    norm <- pnorm(
      pars$U, mean = pars$react_sz, sd = pars$react_sd) -
      pnorm(
        pars$L, mean = pars$react_sz, sd = pars$react_sd)
    dens <- dens / norm
  }

  dens
}


# F-kernel --------------------------------------------------------------------
fyx <- function(
    y, x, pars,
    disturbance = 0,
    disturbance_prev = 0) {

  fl_x(
    x, pars,
    disturbance_prev = disturbance_prev) *
    fl_n_x(
      x, pars,
      disturbance_prev = disturbance_prev) *
    recruits_per_scape(
      pars,
      disturbance = disturbance) *
    re_y_dist(y, pars)
}


# Kernel ----------------------------------------------------------------------
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

  # Active survival and dormancy
  Smat <- sx(
    y, pars,
    disturbance = disturbance)

  Dmat <- dorm_x(
    y, pars,
    disturbance = disturbance)

  # Growth
  Gmat <- matrix(0, n, n)

  Gmat[] <- t(outer(
    y, y,
    Vectorize(function(x, y) {
      gxy(
        x, y, pars,
        disturbance = disturbance)
    }))) * h

  # Correct growth eviction
  for (i in seq_len(n)) {
    if (i <= n / 2) {
      Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    } else {
      Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    }
  }

  # Active -> active
  Tmat <- sweep(
    Gmat,
    2,
    Smat * (1 - Dmat),
    '*')

  # Active -> dormant
  Dorm_row <- matrix(
    Smat * Dmat,
    nrow = 1)

  # Dormant -> active
  React_vec <- react_y_dist(
    y, pars,
    h = h) * h

  React_vec <- React_vec / sum(React_vec)

  p_react <- react_p(
    pars,
    disturbance = disturbance)

  React_col <- matrix(
    p_react * React_vec,
    ncol = 1)

  # Dormant -> dormant
  Dorm_stasis <- 1 - p_react

  # Fertility
  Fmat <- outer(
    y, y,
    Vectorize(function(y, x) {
      fyx(
        y, x, pars,
        disturbance = disturbance,
        disturbance_prev = disturbance_prev)
    })) * h

  # Full 201 x 201 demographic kernel
  K <- rbind(
    cbind(Tmat + Fmat, React_col),
    cbind(Dorm_row, Dorm_stasis))

  list(
    k_yx = K,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    Dorm_row = Dorm_row,
    React_col = React_col,
    Smat = Smat,
    Dmat = Dmat,
    p_reactivate = p_react,
    meshpts = y,
    h = h,
    L = L,
    U = U)
}


# Population growth rate ------------------------------------------------------
lambda_ipm <- function(
    pars,
    disturbance = 0,
    disturbance_prev = 0) {

  Re(eigen(
    kernel(
      pars,
      disturbance = disturbance,
      disturbance_prev = disturbance_prev)$k_yx,
    only.values = TRUE)$values[1])
}


# Environmental-history kernel -----------------------------------------------
kernel_env_hist <- function(
    pars,
    p00 = p_0_given_0,
    p10 = p_1_given_0,
    p01 = p_0_given_1,
    p11 = p_1_given_1) {

  K_00 <- kernel(
    pars = pars,
    disturbance = 0,
    disturbance_prev = 0)$k_yx

  K_10 <- kernel(
    pars = pars,
    disturbance = 1,
    disturbance_prev = 0)$k_yx

  K_01 <- kernel(
    pars = pars,
    disturbance = 0,
    disturbance_prev = 1)$k_yx

  K_11 <- kernel(
    pars = pars,
    disturbance = 1,
    disturbance_prev = 1)$k_yx

  K_env <- rbind(
    cbind(
      p00 * K_00,
      p01 * K_01),
    cbind(
      p10 * K_10,
      p11 * K_11))

  list(
    K_env = K_env,
    K_00 = K_00,
    K_10 = K_10,
    K_01 = K_01,
    K_11 = K_11)
}


# Environmental-history asymptotic lambda ------------------------------------
lambda_env_hist <- function(
    pars,
    p00 = p_0_given_0,
    p10 = p_1_given_0,
    p01 = p_0_given_1,
    p11 = p_1_given_1) {

  K <- kernel_env_hist(
    pars = pars,
    p00 = p00,
    p10 = p10,
    p01 = p01,
    p11 = p11)$K_env

  Re(eigen(
    K,
    only.values = TRUE)$values[1])
}


# Mean IPMs -------------------------------------------------------------------
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 0)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 0)
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 1)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 1)

lambda_env_hist(pars_mean)


# Mean observed population structure -----------------------------------------
make_initial_n_mean <- function(pars, df_init = df_ind) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U

  breaks <- seq(
    L, U,
    length.out = n + 1)

  years_i <- df_init %>%
    filter(state %in% c('active', 'dormant')) %>%
    distinct(year) %>%
    arrange(year) %>%
    pull(year)

  counts_by_year <- map(years_i, function(year_i) {
    df_i <- df_init %>%
      filter(year == year_i)

    active_size <- df_i %>%
      filter(
        state == 'active',
        is.finite(logsize_t0)) %>%
      pull(logsize_t0)

    active_counts <- hist(
      pmin(pmax(active_size, L), U),
      breaks = breaks,
      plot = FALSE,
      include.lowest = TRUE)$counts

    dormant_count <- sum(
      df_i$state == 'dormant',
      na.rm = TRUE)

    c(active_counts, dormant_count)
  })

  Reduce(`+`, counts_by_year) / length(counts_by_year)
}

n_obs <- make_initial_n_mean(pars_mean)


# Initial environmental-history population ----------------------------------
n_env_initial <- c(
  p_prev_0 * n_obs,
  p_prev_1 * n_obs)


# One-step projection ---------------------------------------------------------
project_kernel <- function(K, n_initial, label) {
  n_initial_total <- sum(n_initial)

  n_proj <- as.numeric(
    K %*% n_initial)

  n_projected <- sum(n_proj)

  tibble(
    lambda_type = label,
    asym_lambda = Re(eigen(
      K,
      only.values = TRUE)$values[1]),
    proj_lambda = n_projected / n_initial_total,
    n_initial = n_initial_total,
    n_projected = n_projected)
}


# Mean disturbance-state lambdas ---------------------------------------------
K_mean <- kernel_env_hist(pars_mean)

K_00 <- K_mean$K_00
K_10 <- K_mean$K_10
K_01 <- K_mean$K_01
K_11 <- K_mean$K_11
K_env_hist <- K_mean$K_env

df_lambda_states <- bind_rows(
  project_kernel(
    K_00,
    n_obs,
    'Undisturbed') %>%
    mutate(
      disturbance = 0,
      disturbance_prev = 0),

  project_kernel(
    K_10,
    n_obs,
    'Current-year fire') %>%
    mutate(
      disturbance = 1,
      disturbance_prev = 0),

  project_kernel(
    K_01,
    n_obs,
    'Previous-year fire') %>%
    mutate(
      disturbance = 0,
      disturbance_prev = 1),

  project_kernel(
    K_11,
    n_obs,
    'Fire in both years') %>%
    mutate(
      disturbance = 1,
      disturbance_prev = 1))

df_lambda_states


df_lambda_env_hist <- project_kernel(
  K_env_hist,
  n_env_initial,
  'Environmental historical regime') %>%
  mutate(
    disturbance = NA_real_,
    disturbance_prev = NA_real_)

df_lambda_env_hist


df_lambda_env <- bind_rows(
  df_lambda_states,
  df_lambda_env_hist)

df_lambda_env


# Year-specific IPMs ----------------------------------------------------------
ipm_years <- sort(unique(as.integer(as.character(df_su$year))))

lambda_ipm_year <- function(
    year,
    disturbance = 0,
    disturbance_prev = 0) {

  pars_i <- make_ipm_pars(year)

  lambda_ipm(
    pars_i,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
}


# Year-specific environmental-history lambda ---------------------------------
lambda_env_hist_year <- function(year) {
  pars_i <- make_ipm_pars(year)
  lambda_env_hist(pars_i)
}


lambda_env_hist_year_site <- function(year, site_i) {
  p_site <- df_transition_site %>%
    filter(as.character(site) == as.character(site_i))

  if (nrow(p_site) != 1 ||
      any(!is.finite(c(
        p_site$p00, p_site$p10,
        p_site$p01, p_site$p11)))) {
    return(NA_real_)
  }

  pars_i <- make_ipm_pars(year)

  lambda_env_hist(
    pars = pars_i,
    p00 = p_site$p00,
    p10 = p_site$p10,
    p01 = p_site$p01,
    p11 = p_site$p11)
}


# Year-specific disturbance scenarios ----------------------------------------
lambda_year <- tibble(
  year = ipm_years,
  lambda_undisturbed = map_dbl(
    year, ~ lambda_ipm_year(.x, 0, 0)),
  lambda_current_fire = map_dbl(
    year, ~ lambda_ipm_year(.x, 1, 0)),
  lambda_previous_fire = map_dbl(
    year, ~ lambda_ipm_year(.x, 0, 1)),
  lambda_consecutive_fire = map_dbl(
    year, ~ lambda_ipm_year(.x, 1, 1)),
  lambda_environmental_history = map_dbl(
    year, lambda_env_hist_year))

lambda_year


# Year-specific lambda plot ---------------------------------------------------
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
      lambda_consecutive_fire = 'Fire in both years',
      lambda_environmental_history = 'Environmental-history regime')) %>%
  ggplot(aes(year, lambda, color = disturbance)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c(
    disturbance_case_cols,
    'Environmental-history regime' = 'blue')) +
  theme_bw() +
  labs(
    title = 'Year-specific asymptotic lambda',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = expression(lambda),
    color = 'Disturbance')

fig_lambda_year


# Observed and projected population growth -----------------------------------
# Observed growth follows the ERLO mean IPM: known living individuals are
# counted, and new adults are backfilled one year. One-step IPM projection uses
# the model-resolved active + dormant state distribution in each quadrat.


# Quadrat-year disturbance lookup --------------------------------------------
fire_lookup_quad <- df_re_all %>%
  filter(
    !is.na(site),
    !is.na(pop),
    !is.na(qu),
    !is.na(year)) %>%
  group_by(site, pop, qu, year) %>%
  summarise(
    disturbance_num = first(
      disturbance[!is.na(disturbance)],
      default = NA_real_),
    disturbance_prev_num = first(
      disturbance_prev[!is.na(disturbance_prev)],
      default = NA_real_),
    .groups = 'drop')


# Quadrat-year abundance ------------------------------------------------------
df_counts_alive <- df_ind %>%
  filter(state %in% c('active', 'dormant', 'alive_or_dormant')) %>%
  count(site, pop, qu, year, name = 'n_alive')

df_newadult_backfill <- df_ind %>%
  filter(recruit_type == 'new_adult') %>%
  transmute(site, pop, qu, year = year - 1) %>%
  count(site, pop, qu, year, name = 'n_newadult_backfill')

# Recruitment rows define the monitored quadrat-year inventory.
df_counts_quad <- df_re_all %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_counts_alive,
    by = c('site', 'pop', 'qu', 'year')) %>%
  left_join(
    df_newadult_backfill,
    by = c('site', 'pop', 'qu', 'year')) %>%
  mutate(
    n_alive = replace_na(n_alive, 0L),
    n_newadult_backfill = replace_na(n_newadult_backfill, 0L),
    n = n_alive + n_newadult_backfill) %>%
  group_by(site, pop, qu) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n),
    annual_transition = year_t1 == year + 1) %>%
  ungroup() %>%
  filter(annual_transition)

# The final three transitions cannot fully resolve terminal dormancy.
dormancy_cutoff <- 3
last_complete_transition <- max(df_ind$year, na.rm = TRUE) -
  dormancy_cutoff - 1


# Matched consecutive quadrat transitions ------------------------------------
df_obs_pgr_quad <- df_counts_quad %>%
  filter(
    year <= last_complete_transition,
    n > 0) %>%
  rename(n_t0 = n) %>%
  left_join(
    fire_lookup_quad,
    by = c('site', 'pop', 'qu', 'year')) %>%
  filter(
    !is.na(disturbance_num),
    !is.na(disturbance_prev_num)) %>%
  mutate(obs_pgr = n_t1 / n_t0)


# Initial quadrat state distribution -----------------------------------------
make_initial_n_quad <- function(
    year0,
    site_i,
    pop_i,
    qu_i,
    pars) {

  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U

  breaks <- seq(
    L, U,
    length.out = n + 1)

  df0 <- df_ind %>%
    filter(
      year == year0,
      site == site_i,
      pop == pop_i,
      qu == qu_i)

  active_size <- df0 %>%
    filter(
      state == 'active',
      is.finite(logsize_t0)) %>%
    pull(logsize_t0)

  active_counts <- hist(
    pmin(pmax(active_size, L), U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE)$counts

  dormant_count <- sum(
    df0$state == 'dormant',
    na.rm = TRUE)

  c(active_counts, dormant_count)
}


# Project one quadrat-year ----------------------------------------------------
project_one_quad_year <- function(
    yr,
    site_i,
    pop_i,
    qu_i,
    disturbance_y,
    disturbance_prev_y) {

  pars_y <- make_ipm_pars(yr)

  n_obs_model <- make_initial_n_quad(
    year0 = yr,
    site_i = site_i,
    pop_i = pop_i,
    qu_i = qu_i,
    pars = pars_y)

  n_initial <- sum(n_obs_model)

  if (!is.finite(n_initial) || n_initial <= 0) {
    return(tibble(
      year = yr,
      site = site_i,
      pop = pop_i,
      qu = qu_i,
      disturbance_num = disturbance_y,
      disturbance_prev_num = disturbance_prev_y,
      n_obs_model = n_initial,
      n_proj_model = NA_real_,
      state_asym_lambda = NA_real_,
      proj_lambda = NA_real_))
  }

  K <- kernel(
    pars = pars_y,
    disturbance = disturbance_y,
    disturbance_prev = disturbance_prev_y)$k_yx

  n_proj <- as.numeric(K %*% n_obs_model)
  n_projected <- sum(n_proj)

  tibble(
    year = yr,
    site = site_i,
    pop = pop_i,
    qu = qu_i,
    disturbance_num = disturbance_y,
    disturbance_prev_num = disturbance_prev_y,
    n_obs_model = n_initial,
    n_proj_model = n_projected,
    state_asym_lambda = Re(eigen(
      K,
      only.values = TRUE)$values[1]),
    proj_lambda = n_projected / n_initial)
}


# Project all matched quadrat transitions ------------------------------------
df_proj_quad <- bind_rows(
  lapply(
    seq_len(nrow(df_obs_pgr_quad)),
    function(i) {
      project_one_quad_year(
        yr = df_obs_pgr_quad$year[i],
        site_i = df_obs_pgr_quad$site[i],
        pop_i = df_obs_pgr_quad$pop[i],
        qu_i = df_obs_pgr_quad$qu[i],
        disturbance_y = df_obs_pgr_quad$disturbance_num[i],
        disturbance_prev_y =
          df_obs_pgr_quad$disturbance_prev_num[i])
    }))


# Quadrat-level observed and modeled growth ----------------------------------
df_compare_quad <- df_obs_pgr_quad %>%
  left_join(
    df_proj_quad,
    by = c(
      'year',
      'site',
      'pop',
      'qu',
      'disturbance_num',
      'disturbance_prev_num')) %>%
  mutate(
    error_state_asymptotic_vs_obs =
      state_asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr,
    kernel_state = case_when(
      disturbance_num == 0 & disturbance_prev_num == 0 ~ 'K_00',
      disturbance_num == 1 & disturbance_prev_num == 0 ~ 'K_10',
      disturbance_num == 0 & disturbance_prev_num == 1 ~ 'K_01',
      disturbance_num == 1 & disturbance_prev_num == 1 ~ 'K_11',
      TRUE ~ NA_character_))


# Disturbance-state kernel used for each quadrat-year -------------------------
df_kernel_state_use <- df_compare_quad %>%
  count(
    year,
    kernel_state,
    name = 'n_quadrat_transitions') %>%
  arrange(
    year,
    kernel_state)

df_kernel_state_use


# Site-year environmental-history lambda -------------------------------------
lambda_site_year <- df_compare_quad %>%
  distinct(site, year) %>%
  mutate(
    lambda_environmental_history_site = map2_dbl(
      year,
      site,
      lambda_env_hist_year_site))

lambda_site_year


# Site-year comparison from matched quadrats ---------------------------------
df_compare_site <- df_compare_quad %>%
  group_by(site, year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    n_obs_model = sum(n_obs_model, na.rm = TRUE),
    n_proj_model = sum(n_proj_model, na.rm = TRUE),
    n_quadrats = n(),
    n_burned_quadrats = sum(
      disturbance_num == 1,
      na.rm = TRUE),
    disturbance_num = as.numeric(any(
      disturbance_num == 1,
      na.rm = TRUE)),
    disturbance_prev_num = as.numeric(any(
      disturbance_prev_num == 1,
      na.rm = TRUE)),
    .groups = 'drop') %>%
  left_join(
    lambda_site_year,
    by = c('site', 'year')) %>%
  rename(
    asym_lambda = lambda_environmental_history_site) %>%
  mutate(
    obs_pgr = n_t1 / n_t0,
    proj_lambda = n_proj_model / n_obs_model,
    error_asymptotic_vs_obs =
      asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr)

df_compare_site %>%
  print(n = 100, width = Inf)


# Whole-population annual comparison -----------------------------------------
# Annual totals contain only quadrats represented in both t and t + 1.
df_compare <- df_compare_quad %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    n_obs_model = sum(n_obs_model, na.rm = TRUE),
    n_proj_model = sum(n_proj_model, na.rm = TRUE),
    n_quadrats = n(),
    n_sites = n_distinct(site),
    n_burned_quadrats = sum(
      disturbance_num == 1,
      na.rm = TRUE),
    p_burned_quadrats = mean(
      disturbance_num == 1,
      na.rm = TRUE),
    disturbance = if_else(
      any(disturbance_num == 1, na.rm = TRUE),
      'Fire',
      'No fire'),
    .groups = 'drop') %>%
  left_join(
    lambda_year %>%
      select(
        year,
        lambda_environmental_history),
    by = 'year') %>%
  rename(
    asym_lambda = lambda_environmental_history) %>%
  mutate(
    obs_pgr = n_t1 / n_t0,
    proj_lambda = n_proj_model / n_obs_model,
    disturbance = factor(
      disturbance,
      levels = c('No fire', 'Fire')))


df_compare %>%
  print(n = 100, width = Inf)


# Mean observed population growth --------------------------------------------
lam_obs_y <- df_compare$obs_pgr

lam_obs_mean <- mean(lam_obs_y, na.rm = TRUE)
lam_obs_geo <- exp(mean(log(lam_obs_y), na.rm = TRUE))

lam_obs_mean
lam_obs_geo


# Modeled versus observed mean population growth -----------------------------
c(
  asym_lambda_env = lambda_env_hist(pars_mean),
  projected_lambda_env = df_lambda_env_hist$proj_lambda,
  lambda_obs_arithmetic = lam_obs_mean,
  lambda_obs_geometric = lam_obs_geo)


df_lambda_mean <- tibble(
  estimate = c(
    'IPM asymptotic',
    'IPM one-step projected',
    'Observed arithmetic',
    'Observed geometric'),
  lambda = c(
    lambda_env_hist(pars_mean),
    df_lambda_env_hist$proj_lambda,
    lam_obs_mean,
    lam_obs_geo))

df_lambda_mean


# Observed vs modeled plot ----------------------------------------------------
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
      asym_lambda = 'Environmental-history asymptotic lambda',
      proj_lambda =
        'Projected lambda from observed quadrat state distributions'))

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


# Log-transformed observed vs modeled plot -----------------------------------
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


# Summary statistics ----------------------------------------------------------
df_compare_summary <- df_compare %>%
  summarise(
    n_years = n(),
    arithmetic_mean_obs_pgr =
      mean(obs_pgr, na.rm = TRUE),
    geometric_mean_obs_pgr =
      exp(mean(log(obs_pgr), na.rm = TRUE)),
    arithmetic_mean_env_hist_lambda =
      mean(asym_lambda, na.rm = TRUE),
    geometric_mean_env_hist_lambda =
      exp(mean(log(asym_lambda), na.rm = TRUE)),
    arithmetic_mean_proj_lambda =
      mean(proj_lambda, na.rm = TRUE),
    geometric_mean_proj_lambda =
      exp(mean(log(proj_lambda), na.rm = TRUE)),
    mean_error_env_hist_vs_obs =
      mean(
        asym_lambda - obs_pgr,
        na.rm = TRUE),
    mean_error_projected_vs_obs =
      mean(
        proj_lambda - obs_pgr,
        na.rm = TRUE),
    percent_bias_env_hist_vs_obs =
      100 * sum(
        asym_lambda - obs_pgr,
        na.rm = TRUE) /
      sum(
        obs_pgr,
        na.rm = TRUE),
    percent_bias_projected_vs_obs =
      100 * sum(
        proj_lambda - obs_pgr,
        na.rm = TRUE) /
      sum(
        obs_pgr,
        na.rm = TRUE),
    rmse_env_hist_vs_obs =
      sqrt(mean(
        (asym_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    rmse_projected_vs_obs =
      sqrt(mean(
        (proj_lambda - obs_pgr)^2,
        na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'statistic',
    values_to = 'value')

df_compare_summary


# Log-scale summary statistics ------------------------------------------------
df_compare_log_summary <- df_compare %>%
  filter(
    obs_pgr > 0,
    asym_lambda > 0,
    proj_lambda > 0) %>%
  mutate(
    log_obs_pgr = log(obs_pgr),
    log_asym_lambda = log(asym_lambda),
    log_proj_lambda = log(proj_lambda),
    log_error_asym = log_asym_lambda - log_obs_pgr,
    log_error_proj = log_proj_lambda - log_obs_pgr) %>%
  summarise(
    n_years = n(),
    mean_log_obs_pgr = mean(log_obs_pgr, na.rm = TRUE),
    mean_log_env_hist_lambda = mean(log_asym_lambda, na.rm = TRUE),
    mean_log_proj_lambda = mean(log_proj_lambda, na.rm = TRUE),
    mean_log_error_env_hist = mean(log_error_asym, na.rm = TRUE),
    mean_log_error_projected = mean(log_error_proj, na.rm = TRUE),
    mean_absolute_log_error_env_hist =
      mean(abs(log_error_asym), na.rm = TRUE),
    mean_absolute_log_error_projected =
      mean(abs(log_error_proj), na.rm = TRUE),
    rmse_log_env_hist =
      sqrt(mean(log_error_asym^2, na.rm = TRUE)),
    rmse_log_projected =
      sqrt(mean(log_error_proj^2, na.rm = TRUE)),
    multiplicative_bias_env_hist =
      exp(mean(log_error_asym, na.rm = TRUE)),
    multiplicative_bias_projected =
      exp(mean(log_error_proj, na.rm = TRUE)),
    percent_multiplicative_bias_env_hist =
      100 * (exp(mean(log_error_asym, na.rm = TRUE)) - 1),
    percent_multiplicative_bias_projected =
      100 * (exp(mean(log_error_proj, na.rm = TRUE)) - 1)) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'statistic',
    values_to = 'value')

df_compare_log_summary


# Site-level summary ----------------------------------------------------------
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
    lambda_obs_geometric =
      exp(mean(log(obs_pgr), na.rm = TRUE)),
    lambda_obs_arithmetic =
      mean(obs_pgr, na.rm = TRUE),
    lambda_env_hist_geometric =
      exp(mean(log(asym_lambda), na.rm = TRUE)),
    lambda_env_hist_arithmetic =
      mean(asym_lambda, na.rm = TRUE),
    lambda_projected_geometric =
      exp(mean(log(proj_lambda), na.rm = TRUE)),
    lambda_projected_arithmetic =
      mean(proj_lambda, na.rm = TRUE),
    error_env_hist_geo_vs_obs_geo =
      lambda_env_hist_geometric -
      lambda_obs_geometric,
    error_projected_geo_vs_obs_geo =
      lambda_projected_geometric -
      lambda_obs_geometric,
    rmse_env_hist_vs_obs =
      sqrt(mean(
        (asym_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    rmse_projected_vs_obs =
      sqrt(mean(
        (proj_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    mean_n_initial =
      mean(n_obs_model, na.rm = TRUE),
    mean_disturbance =
      mean(disturbance_num, na.rm = TRUE),
    .groups = 'drop') %>%
  arrange(as.numeric(as.character(site)))

df_compare_site_summary %>%
  print(n = 100, width = Inf)


# Projection-state coverage diagnostic ---------------------------------------
# Observed counts include alive_or_dormant individuals and backfilled new
# adults. The initial IPM vector contains only model-resolved active/dormant
# states. This table makes any difference explicit instead of hiding it.
df_projection_state_coverage <- df_compare_quad %>%
  mutate(
    n_not_in_initial_ipm = n_t0 - n_obs_model,
    p_initial_ipm_covered = if_else(
      n_t0 > 0,
      n_obs_model / n_t0,
      NA_real_)) %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_obs_model = sum(n_obs_model, na.rm = TRUE),
    n_not_in_initial_ipm = sum(n_not_in_initial_ipm, na.rm = TRUE),
    p_initial_ipm_covered = n_obs_model / n_t0,
    .groups = 'drop')

df_projection_state_coverage


# Yearly lambda comparison diagnostic ----------------------------------------
df_lambda_year_check <- lambda_year %>%
  mutate(
    min_state_lambda = pmin(
      lambda_undisturbed,
      lambda_current_fire,
      lambda_previous_fire,
      lambda_consecutive_fire),
    max_state_lambda = pmax(
      lambda_undisturbed,
      lambda_current_fire,
      lambda_previous_fire,
      lambda_consecutive_fire),
    env_hist_minus_undisturbed =
      lambda_environmental_history -
      lambda_undisturbed) %>%
  select(
    year,
    lambda_undisturbed,
    lambda_current_fire,
    lambda_previous_fire,
    lambda_consecutive_fire,
    lambda_environmental_history,
    min_state_lambda,
    max_state_lambda,
    env_hist_minus_undisturbed)

df_lambda_year_check %>%
  as_tibble() %>%
  print(n = 100, width = Inf)


# Fire-effect kernel diagnostics ----------------------------------------------
# These checks separate current-fire and previous-fire effects on the mean
# demographic kernels. Survival/growth/dormancy/reactivation use current fire;
# flowering/scapes use previous fire; recruitment conversion uses current fire.
c(
  current_fire_no_previous =
    max(abs(K_00 - K_10)),
  previous_fire_no_current =
    max(abs(K_00 - K_01)),
  both_vs_undisturbed =
    max(abs(K_00 - K_11)))

c(
  current_fire_total =
    max(abs(colSums(K_00) - colSums(K_10))),
  previous_fire_total =
    max(abs(colSums(K_00) - colSums(K_01))))


# # Save key outputs ----------------------------------------------------------
# # Uncomment if needed.
# saveRDS(
#   pars_mean,
#   file.path(dir_result, 'erlo_year_fire_pars_mean.rds'))
# saveRDS(
#   lambda_year,
#   file.path(dir_result, 'erlo_year_fire_lambda_year.rds'))
# saveRDS(
#   df_compare,
#   file.path(dir_result, 'erlo_year_fire_lambda_compare.rds'))
# saveRDS(
#   df_compare_site,
#   file.path(dir_result, 'erlo_year_fire_lambda_compare_site.rds'))
