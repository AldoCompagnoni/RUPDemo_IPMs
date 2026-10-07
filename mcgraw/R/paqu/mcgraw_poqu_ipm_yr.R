# IPM year-specific; dormancy - McGraw - Panax quinquefolius

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.10.07

# Study organism: Panax quinquefolius L.
# Time period: 1998-2016


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
set.seed(100)
options(stringsAsFactors = FALSE)


# Packages --------------------------------------------------------------------
source('C:/code/RUPDemo_IPMs/helper_functions/load_packages.R')
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
v_head <- c('mcgraw')
v_species <- c('Panax quinquefolius')
custom_delimiter <- c()
v_years_re <- c()
v_states <- c('WV')

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Empty = AICc selection; 0:3 = fixed polynomial degree.
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
v_mod_set_do   <- c()
v_mod_set_fl   <- c()
v_mod_set_se   <- c()
v_mod_set_se_n <- c()


# Directory -------------------------------------------------------------------
dir_wd <- file.path('C:/code/RUPDemo_IPMs')
dir_pub <- file.path(dir_wd, v_head)
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
source(file.path(dir_wd, 'helper_functions/plot_binned_prop.R'))
source(file.path(dir_wd, 'helper_functions/plot_binned_prop_year.R'))
source(file.path(dir_wd, 'helper_functions/line_color_pred_fun.R'))
source(file.path(dir_wd, 'helper_functions/predictor_fun.R'))


# Load workdata ---------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0(v_head, '_', v_sp_abb, '_df_workdata.csv'))) %>%
  mutate(year = as.integer(year), year_t1 = as.integer(year_t1)) %>%
  filter(state %in% v_states, !is.na(year), !(year %in% v_years_re)) %>%
  # Remove the inspected duplicate empty 2001 -> 2001 row only.
  filter(!(population == 30 & id == 354 & year == 2001 &
             year_t1 %in% 2001L))


# Controls --------------------------------------------------------------------
ctrl_glmer <- glmerControl(
  optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))

ctrl_lmer <- lmerControl(
  optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5))


# Survival data ---------------------------------------------------------------
df_su <- df %>%
  filter(
    reliable_demography, dormant_t0 == 0, !is.na(survives),
    size_t0 > 0, is.finite(logsize_t0),
    is.finite(logsize_t0_2), is.finite(logsize_t0_3),
    !(year_t1 %in% v_years_re)) %>%
  mutate(year = factor(year)) %>%
  dplyr::select(
    state, population, id, year, size_t0, survives,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Survival model ---------------------------------------------------------------
mod_su_0 <- glmer(
  survives ~ 1 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_1 <- glmer(
  survives ~ logsize_t0 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
  data = df_su, family = binomial, control = ctrl_glmer)

mods_su <- list(
  mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAICc <- bbmle::AICctab(
  mods_su, weights = TRUE, sort = FALSE)$dAICc
mods_su_sorted <- order(mods_su_dAICc)

if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index <- mod_su_index_bestfit - 1
} else {
  mod_su_index_bestfit <- v_mod_set_su + 1
  v_mod_su_index <- v_mod_set_su}

mod_su_best <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef <- coef(mod_su_best)$year

mod_su_best
summary(mod_su_best)
mods_su_dAICc


# Survival plots by year -------------------------------------------------------
make_su_year_plot <- function(year_i) {
  df_i <- df_su %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_su$year)))

  pred_i$survives <- predict(
    mod_su_best, newdata = pred_i, type = 'response', re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, survives), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, survives), linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'Survival probability') +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

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


# Growth data -----------------------------------------------------------------
df_gr <- df %>%
  filter(
    reliable_demography, survives == 1,
    dormant_t0 == 0, dormant_t1 == 0,
    size_t0 > 0, size_t1 > 0,
    is.finite(logsize_t0), is.finite(logsize_t1),
    is.finite(logsize_t0_2), is.finite(logsize_t0_3),
    !(year_t1 %in% v_years_re)) %>%
  mutate(year = factor(year)) %>%
  dplyr::select(
    state, population, id, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

fig_gr_raw <- ggplot(df_gr, aes(logsize_t0, logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  labs(
    title = 'Growth', subtitle = v_ggp_suffix,
    x = expression('log(leaf area)'[t0]),
    y = expression('log(leaf area)'[t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_raw


# Growth model -----------------------------------------------------------------
mod_gr_0 <- lmer(
  logsize_t1 ~ 1 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (logsize_t0 | year),
  data = df_gr, REML = FALSE, control = ctrl_lmer)

mods_gr <- list(
  mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAICc <- bbmle::AICctab(
  mods_gr, weights = TRUE, sort = FALSE)$dAICc
mods_gr_sorted <- order(mods_gr_dAICc)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index <- mod_gr_index_bestfit - 1
} else {
  mod_gr_index_bestfit <- v_mod_set_gr + 1
  v_mod_gr_index <- v_mod_set_gr}

mod_gr_best <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef <- coef(mod_gr_best)$year

mod_gr_best
summary(mod_gr_best)
mods_gr_dAICc


# Growth plots by year ---------------------------------------------------------
make_gr_year_plot <- function(year_i) {
  df_i <- df_gr %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_gr$year)))

  pred_i$logsize_t1 <- predict(
    mod_gr_best, newdata = pred_i, re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      logsize_t1 = mean(logsize_t1, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, logsize_t1), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, logsize_t1), linewidth = 0.7) +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'Next log leaf area') +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

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


# Growth variance -------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_best)
mod_gr_y <- resid(mod_gr_best)^2
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Dormancy entry data ---------------------------------------------------------
df_do <- df %>%
  filter(
    reliable_demography, dormant_t0 == 0, survives == 1,
    !is.na(dormant_t1), size_t0 > 0, is.finite(logsize_t0),
    is.finite(logsize_t0_2), is.finite(logsize_t0_3),
    !(year_t1 %in% v_years_re)) %>%
  mutate(year = factor(year), enter_dormancy = dormant_t1) %>%
  dplyr::select(
    state, population, id, year, size_t0, enter_dormancy,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Dormancy entry model ---------------------------------------------------------
mod_do_0 <- glmer(
  enter_dormancy ~ 1 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_1 <- glmer(
  enter_dormancy ~ logsize_t0 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_2 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mod_do_3 <- glmer(
  enter_dormancy ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
  data = df_do, family = binomial, control = ctrl_glmer)

mods_do <- list(
  mod_do_0, mod_do_1, mod_do_2, mod_do_3)
mods_do_dAICc <- bbmle::AICctab(
  mods_do, weights = TRUE, sort = FALSE)$dAICc
mods_do_sorted <- order(mods_do_dAICc)

if (length(v_mod_set_do) == 0) {
  mod_do_index_bestfit <- mods_do_sorted[1]
  v_mod_do_index <- mod_do_index_bestfit - 1
} else {
  mod_do_index_bestfit <- v_mod_set_do + 1
  v_mod_do_index <- v_mod_set_do}

mod_do_best <- mods_do[[mod_do_index_bestfit]]
mod_do_ranef <- coef(mod_do_best)$year

mod_do_best
summary(mod_do_best)
mods_do_dAICc


# Dormancy entry plots by year -------------------------------------------------
make_do_year_plot <- function(year_i) {
  df_i <- df_do %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_do$year)))

  pred_i$enter_dormancy <- predict(
    mod_do_best, newdata = pred_i, type = 'response', re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      enter_dormancy = mean(enter_dormancy, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, enter_dormancy), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, enter_dormancy), linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'Dormancy entry given survival') +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

do_yrs <- lapply(levels(df_do$year), make_do_year_plot)
fig_do_years <- wrap_plots(do_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Dormancy entry - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_do_years


# Dormant survival and reactivation -------------------------------------------
df_dorm_su <- df %>%
  filter(
    reliable_demography, dormant_t0 == 1, !is.na(survives),
    !(year_t1 %in% v_years_re)) %>%
  mutate(year = factor(year))

# A constant response cannot estimate a year random effect.
if (n_distinct(df_dorm_su$survives) > 1) {
  mod_dorm_su <- glmer(
    survives ~ 1 + (1 | year),
    data = df_dorm_su, family = binomial, control = ctrl_glmer)
} else {
  mod_dorm_su <- glm(
    survives ~ 1, data = df_dorm_su, family = binomial)}

mod_dorm_su
summary(mod_dorm_su)


df_ra <- df_dorm_su %>%
  filter(survives == 1, !is.na(dormant_t1)) %>%
  mutate(
    year = factor(as.character(year)),
    reactivate = as.integer(dormant_t1 == 0))

if (n_distinct(df_ra$reactivate) > 1) {
  mod_ra <- glmer(
    reactivate ~ 1 + (1 | year),
    data = df_ra, family = binomial, control = ctrl_glmer)
} else {
  mod_ra <- glm(
    reactivate ~ 1, data = df_ra, family = binomial)}

mod_ra
summary(mod_ra)


df_ra_size <- df %>%
  filter(
    reliable_demography, dormant_t0 == 1, survives == 1,
    dormant_t1 == 0, size_t1 > 0, is.finite(logsize_t1),
    !(year_t1 %in% v_years_re)) %>%
  mutate(logsize_reactivate = logsize_t1)

react_sz <- mean(df_ra_size$logsize_reactivate, na.rm = TRUE)
react_sd <- sd(df_ra_size$logsize_reactivate, na.rm = TRUE)


# Flower data -----------------------------------------------------------------
df_fl <- df %>%
  filter(
    !is.na(flower), size_t0 > 0, is.finite(logsize_t0),
    is.finite(logsize_t0_2), is.finite(logsize_t0_3)) %>%
  mutate(year = factor(year)) %>%
  dplyr::select(
    state, population, id, year, size_t0, flower,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Flower model -----------------------------------------------------------------
mod_fl_0 <- glmer(
  flower ~ 1 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_1 <- glmer(
  flower ~ logsize_t0 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
  data = df_fl, family = binomial, control = ctrl_glmer)

mods_fl <- list(
  mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAICc <- bbmle::AICctab(
  mods_fl, weights = TRUE, sort = FALSE)$dAICc
mods_fl_sorted <- order(mods_fl_dAICc)

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index <- mod_fl_index_bestfit - 1
} else {
  mod_fl_index_bestfit <- v_mod_set_fl + 1
  v_mod_fl_index <- v_mod_set_fl}

mod_fl_best <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef <- coef(mod_fl_best)$year

mod_fl_best
summary(mod_fl_best)
mods_fl_dAICc


# Flowering plots by year ------------------------------------------------------
make_fl_year_plot <- function(year_i) {
  df_i <- df_fl %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_fl$year)))

  pred_i$flower <- predict(
    mod_fl_best, newdata = pred_i, type = 'response', re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      flower = mean(flower, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, flower), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, flower), linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'Flowering probability') +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

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


# Seed production data --------------------------------------------------------
# Probability of positive seed production, conditional on flowering.
df_se <- df %>%
  filter(
    flower == 1, !is.na(seed_nr), size_t0 > 0,
    is.finite(logsize_t0), is.finite(logsize_t0_2),
    is.finite(logsize_t0_3)) %>%
  mutate(year = factor(year), seed_prod = as.integer(seed_nr > 0)) %>%
  dplyr::select(
    state, population, id, year, size_t0, flower, seed_prod, seed_nr,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Seed production probability model --------------------------------------------
mod_se_0 <- glmer(
  seed_prod ~ 1 + (1 | year),
  data = df_se, family = binomial, control = ctrl_glmer)

mod_se_1 <- glmer(
  seed_prod ~ logsize_t0 + (1 | year),
  data = df_se, family = binomial, control = ctrl_glmer)

mod_se_2 <- glmer(
  seed_prod ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_se, family = binomial, control = ctrl_glmer)

mod_se_3 <- glmer(
  seed_prod ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
  data = df_se, family = binomial, control = ctrl_glmer)

mods_se <- list(
  mod_se_0, mod_se_1, mod_se_2, mod_se_3)
mods_se_dAICc <- bbmle::AICctab(
  mods_se, weights = TRUE, sort = FALSE)$dAICc
mods_se_sorted <- order(mods_se_dAICc)

if (length(v_mod_set_se) == 0) {
  mod_se_index_bestfit <- mods_se_sorted[1]
  v_mod_se_index <- mod_se_index_bestfit - 1
} else {
  mod_se_index_bestfit <- v_mod_set_se + 1
  v_mod_se_index <- v_mod_set_se}

mod_se_best <- mods_se[[mod_se_index_bestfit]]
mod_se_ranef <- coef(mod_se_best)$year

mod_se_best
summary(mod_se_best)
mods_se_dAICc


# Seed production plots by year ------------------------------------------------
make_se_year_plot <- function(year_i) {
  df_i <- df_se %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_se$year)))

  pred_i$seed_prod <- predict(
    mod_se_best, newdata = pred_i, type = 'response', re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      seed_prod = mean(seed_prod, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, seed_prod), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, seed_prod), linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'P(positive seeds | flowering)') +
    ylim(0, 1) +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

se_yrs <- lapply(levels(df_se$year), make_se_year_plot)
fig_se_years <- wrap_plots(se_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Seed production - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_se_years


# Seed number conditional on seed production ----------------------------------
df_se_n <- df_se %>%
  filter(seed_prod == 1, seed_nr > 0) %>%
  mutate(year = factor(as.character(year)))

# Retain the ordinary negative-binomial model from PAQU mean.


# Seed number model ------------------------------------------------------------
mod_se_n_0 <- glmer.nb(
  seed_nr ~ 1 + (1 | year),
  data = df_se_n, control = ctrl_glmer)

mod_se_n_1 <- glmer.nb(
  seed_nr ~ logsize_t0 + (1 | year),
  data = df_se_n, control = ctrl_glmer)

mod_se_n_2 <- glmer.nb(
  seed_nr ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_se_n, control = ctrl_glmer)

# mod_se_n_3 <- glmer.nb(
#   seed_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
#   data = df_se_n, control = ctrl_glmer)

mods_se_n <- list(
  mod_se_n_0, mod_se_n_1, mod_se_n_2) # , mod_se_n_3
mods_se_n_dAICc <- bbmle::AICctab(
  mods_se_n, weights = TRUE, sort = FALSE)$dAICc
mods_se_n_sorted <- order(mods_se_n_dAICc)

if (length(v_mod_set_se_n) == 0) {
  mod_se_n_index_bestfit <- mods_se_n_sorted[1]
  v_mod_se_n_index <- mod_se_n_index_bestfit - 1
} else {
  mod_se_n_index_bestfit <- v_mod_set_se_n + 1
  v_mod_se_n_index <- v_mod_set_se_n}

mod_se_n_best <- mods_se_n[[mod_se_n_index_bestfit]]
mod_se_n_ranef <- coef(mod_se_n_best)$year

mod_se_n_best
summary(mod_se_n_best)
mods_se_n_dAICc


# Seed number plots by year ----------------------------------------------------
make_se_n_year_plot <- function(year_i) {
  df_i <- df_se_n %>% filter(year == year_i)
  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE), length.out = 100)

  pred_i <- tibble(logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(year_i, levels = levels(df_se_n$year)))

  pred_i$seed_nr <- predict(
    mod_se_n_best, newdata = pred_i, type = 'response', re.form = NULL)

  pts_i <- df_i %>%
    mutate(bin = cut(logsize_t0, breaks = 8, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      seed_nr = mean(seed_nr, na.rm = TRUE),
      n = n(), .groups = 'drop') %>%
    filter(!is.na(logsize_t0), n > 0)

  ggplot() +
    geom_point(
      data = pts_i, aes(logsize_t0, seed_nr), size = 1.1) +
    geom_line(
      data = pred_i, aes(logsize_t0, seed_nr), linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = 'Seeds given positive production') +
    theme_bw() +
    theme(text = element_text(size = 5), legend.position = 'none')}

se_n_yrs <- lapply(levels(df_se_n$year), make_se_n_year_plot)
fig_se_n_years <- wrap_plots(se_n_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Seed number - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_se_n_years


# Seed to recruit transition --------------------------------------------------
# Seeds produced in year t are related to recruits observed in year t + 2.
df_re_seed <- df %>%
  group_by(population, year) %>%
  summarise(
    nr_seeds = if (all(is.na(seed_nr))) NA_real_ else
      sum(seed_nr, na.rm = TRUE),
    nr_seed_obs = sum(!is.na(seed_nr)), .groups = 'drop') %>%
  mutate(year_t2 = year + 2L)


df_re_sampled <- df %>%
  distinct(population, year)

df_re_count <- df %>%
  filter(recruit == 1) %>%
  count(population, year, name = 'nr_recruits')

df_re_count <- df_re_sampled %>%
  left_join(df_re_count, by = c('population', 'year')) %>%
  mutate(nr_recruits = replace_na(nr_recruits, 0L))


df_re <- df_re_seed %>%
  inner_join(df_re_count, by = c('population', 'year_t2' = 'year')) %>%
  filter(nr_seed_obs > 0) %>%
  mutate(recruit_year = factor(year_t2))


# Seed to recruit model -------------------------------------------------------
# The offset gives recruits per seed; the year effect is for recruit arrival.
df_re_mod <- df_re %>%
  filter(nr_seeds > 0) %>%
  mutate(recruit_year = factor(as.character(recruit_year)))

mod_re_best <- glmer.nb(
  nr_recruits ~ 1 + offset(log(nr_seeds)) + (1 | recruit_year),
  data = df_re_mod, control = ctrl_glmer)

mod_re_ranef <- coef(mod_re_best)$recruit_year
recr_per_seed <- exp(unname(fixef(mod_re_best)['(Intercept)']))

mod_re_best
summary(mod_re_best)
recr_per_seed


# Recruitment plots by recruit year -------------------------------------------
make_re_year_plot <- function(year_i) {
  df_i <- df_re %>% filter(recruit_year == year_i)
  pred_i <- tibble(
    nr_seeds = seq(1, max(df_i$nr_seeds), length.out = 100),
    recruit_year = factor(year_i, levels = levels(df_re_mod$recruit_year)))

  pred_i$pred_recruits <- predict(
    mod_re_best, newdata = pred_i, type = 'response', re.form = NULL)

  ggplot() +
    geom_point(
      data = df_i,
      aes(nr_seeds, nr_recruits, color = factor(population)),
      alpha = 0.75, size = 2) +
    geom_line(
      data = pred_i, aes(nr_seeds, pred_recruits), linewidth = 0.8) +
    theme_bw() +
    labs(
      title = year_i,
      x = 'Seeds two years earlier', y = 'Recruits', color = 'Population') +
    theme(text = element_text(size = 5), legend.position = 'none')}

re_yrs <- lapply(levels(df_re_mod$recruit_year), make_re_year_plot)
fig_re_years <- wrap_plots(re_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Seed to recruit transition - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = 'bold'),
      plot.subtitle = element_text(size = 9)))

fig_re_years


# Recruit size distribution ---------------------------------------------------
df_re_size <- df %>%
  filter(recruit == 1, size_t0 > 0, is.finite(logsize_t0))

recr_sz <- mean(df_re_size$logsize_t0, na.rm = TRUE)
recr_sd <- sd(df_re_size$logsize_t0, na.rm = TRUE)


# Parameter extraction helpers ------------------------------------------------
empty_coef_df <- function() {
  data.frame(coefficient = character(), value = numeric())}

is_mixed_model <- function(model) {
  inherits(model, 'merMod')}

get_fixef_safe <- function(model) {
  if (is_mixed_model(model)) fixef(model) else coef(model)}

get_group_coef_safe <- function(model, group_var) {
  if (!is_mixed_model(model)) return(NULL)
  coef(model)[[group_var]]}

find_model_term <- function(model_terms, aliases) {
  for (alias in unlist(aliases)) {
    if (grepl('^regex:', alias)) {
      hit <- grep(sub('^regex:', '', alias), model_terms, value = TRUE)
    } else {
      hit <- model_terms[model_terms == alias]}
    if (length(hit) > 0) return(hit[1])}
  NA_character_}

bind_coef_rows <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0) return(empty_coef_df())
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out}

coef_df_to_list <- function(x) {
  if (nrow(x) == 0) return(list())
  x %>%
    mutate(coefficient = as.character(coefficient)) %>%
    pivot_wider(names_from = coefficient, values_from = value) %>%
    as.list()}

extract_fixed_pars <- function(model, term_map, fill_missing = TRUE) {
  model_coefs <- get_fixef_safe(model)
  out <- lapply(names(term_map), function(coef_name) {
    model_term <- find_model_term(names(model_coefs), term_map[[coef_name]])
    if (is.na(model_term)) {
      if (!fill_missing) return(NULL)
      value <- 0
    } else {
      value <- unname(model_coefs[model_term])}
    data.frame(coefficient = coef_name, value = value)})
  bind_coef_rows(out)}

extract_group_pars <- function(
    model, group_var, term_map, fill_missing = TRUE) {
  coef_matrix <- get_group_coef_safe(model, group_var)
  if (is.null(coef_matrix)) return(empty_coef_df())
  out <- lapply(names(term_map), function(coef_prefix) {
    model_term <- find_model_term(
      colnames(coef_matrix), term_map[[coef_prefix]])
    if (is.na(model_term)) {
      if (!fill_missing) return(NULL)
      value <- rep(0, nrow(coef_matrix))
    } else {
      value <- coef_matrix[, model_term]}
    data.frame(
      coefficient = paste0(coef_prefix, '_', rownames(coef_matrix)),
      value = value)})
  bind_coef_rows(out)}

extract_re_ranef_devs <- function(model, group_var, group_label) {
  ranef_matrix <- ranef(model)[[group_var]]
  data.frame(
    coefficient = paste0('re_u0_', group_label, '_', rownames(ranef_matrix)),
    value = ranef_matrix[, '(Intercept)'])}


# Term maps -------------------------------------------------------------------
size_term_map <- list(
  b0 = c('(Intercept)'),
  b1 = c('logsize_t0'),
  b2 = c('logsize_t0_2', 'I(logsize_t0^2)'),
  b3 = c('logsize_t0_3', 'I(logsize_t0^3)'))

make_term_map <- function(prefix, map) {
  out <- map
  names(out) <- paste0(prefix, names(out))
  out}

su_term_map <- make_term_map('surv_', size_term_map)
gr_term_map <- make_term_map('grow_', size_term_map)
do_term_map <- make_term_map('dorm_', size_term_map)
fl_term_map <- make_term_map('fl_', size_term_map)
se_term_map <- make_term_map('se_', size_term_map)
se_n_term_map <- make_term_map('se_n_', size_term_map)
re_term_map <- list(re_b0 = c('(Intercept)'))


# Fixed parameters ------------------------------------------------------------
su_fe <- extract_fixed_pars(mod_su_best, su_term_map)
gr_fe <- extract_fixed_pars(mod_gr_best, gr_term_map)
do_fe <- extract_fixed_pars(mod_do_best, do_term_map)
fl_fe <- extract_fixed_pars(mod_fl_best, fl_term_map)
se_fe <- extract_fixed_pars(mod_se_best, se_term_map)
se_n_fe <- extract_fixed_pars(mod_se_n_best, se_n_term_map)
re_fe <- extract_fixed_pars(mod_re_best, re_term_map)

p_dorm_survival <- plogis(
  unname(get_fixef_safe(mod_dorm_su)['(Intercept)']))
p_reactivate <- plogis(
  unname(get_fixef_safe(mod_ra)['(Intercept)']))


# Constants -------------------------------------------------------------------
gr_var_coef <- coef(mod_gr_var)
mesh_limits <- range(
  c(df$logsize_t0, df_gr$logsize_t1, df_ra_size$logsize_reactivate),
  na.rm = TRUE, finite = TRUE)

constants <- tibble::tribble(
  ~coefficient, ~value,
  'recr_sz', recr_sz,
  'recr_sd', recr_sd,
  'react_sz', react_sz,
  'react_sd', react_sd,
  'p_dorm_survival', p_dorm_survival,
  'p_reactivate', p_reactivate,
  'a', unname(gr_var_coef['a']),
  'b', unname(gr_var_coef['b']),
  'L', mesh_limits[1],
  'U', mesh_limits[2],
  'mat_siz', 200,
  'mod_su_index', v_mod_su_index,
  'mod_gr_index', v_mod_gr_index,
  'mod_do_index', v_mod_do_index,
  'mod_fl_index', v_mod_fl_index,
  'mod_se_index', v_mod_se_index,
  'mod_se_n_index', v_mod_se_n_index) %>%
  mutate(coefficient = as.character(coefficient), value = as.numeric(value))

pars_cons <- bind_coef_rows(list(
  su_fe, gr_fe, do_fe, fl_fe, se_fe, se_n_fe, re_fe, constants))
pars_all_mean <- coef_df_to_list(pars_cons)
pars_mean <- pars_all_mean


# Year-varying parameters -----------------------------------------------------
su_out_yr <- extract_group_pars(mod_su_best, 'year', su_term_map)
gr_out_yr <- extract_group_pars(mod_gr_best, 'year', gr_term_map)
do_out_yr <- extract_group_pars(mod_do_best, 'year', do_term_map)
fl_out_yr <- extract_group_pars(mod_fl_best, 'year', fl_term_map)
se_out_yr <- extract_group_pars(mod_se_best, 'year', se_term_map)
se_n_out_yr <- extract_group_pars(mod_se_n_best, 'year', se_n_term_map)

dorm_su_out_yr <- extract_group_pars(
  mod_dorm_su, 'year', list(p_dorm_survival = '(Intercept)')) %>%
  mutate(value = plogis(value))
ra_out_yr <- extract_group_pars(
  mod_ra, 'year', list(p_reactivate = '(Intercept)')) %>%
  mutate(value = plogis(value))

pars_var <- bind_coef_rows(list(
  su_out_yr, gr_out_yr, do_out_yr, fl_out_yr, se_out_yr, se_n_out_yr,
  dorm_su_out_yr, ra_out_yr))
pars_all_year <- coef_df_to_list(pars_var)


# Recruitment-year random-effect deviations -----------------------------------
pars_re_year <- extract_re_ranef_devs(
  mod_re_best, 'recruit_year', 'recruit_year')
pars_all_re_year <- coef_df_to_list(pars_re_year)


# IPM helper functions --------------------------------------------------------
inv_logit <- function(x) {
  plogis(x)}

get_par <- function(pars, par_name, default = 0) {
  if (!is.null(pars[[par_name]])) return(pars[[par_name]])
  default}

make_ipm_pars <- function(
    pars_mean, pars_year = NULL, pars_re_year = NULL,
    year = NULL, recruit_year = NULL) {
  pars <- pars_mean

  if (!is.null(year) && !is.null(pars_year)) {
    year <- as.character(year)
    year_hits <- grep(paste0('_', year, '$'), names(pars_year), value = TRUE)
    for (nm in year_hits) {
      base_nm <- sub(paste0('_', year, '$'), '', nm)
      pars[[base_nm]] <- pars_year[[nm]]}}

  # K_t recruits the seed cohort from t-1 into plants at t+1.
  if (is.null(recruit_year) && !is.null(year)) {
    recruit_year <- as.integer(year) + 1L}

  if (!is.null(recruit_year) && !is.null(pars_re_year)) {
    pars$re_b0 <- get_par(pars, 're_b0') + get_par(
      pars_re_year, paste0('re_u0_recruit_year_', recruit_year))}

  pars}

surv_lp <- function(x, pars) {
  get_par(pars, 'surv_b0') + get_par(pars, 'surv_b1') * x +
    get_par(pars, 'surv_b2') * x^2 + get_par(pars, 'surv_b3') * x^3}

sx <- function(x, pars) {
  inv_logit(surv_lp(x, pars))}

grow_mu <- function(x, pars) {
  get_par(pars, 'grow_b0') + get_par(pars, 'grow_b1') * x +
    get_par(pars, 'grow_b2') * x^2 + get_par(pars, 'grow_b3') * x^3}

# Variance was fitted against fitted next-year size, not starting size.
grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * grow_mu(x, pars)))}

gxy <- function(x, y, pars) {
  dnorm(y, mean = grow_mu(x, pars), sd = grow_sd(x, pars))}

dx <- function(x, pars) {
  eta <- get_par(pars, 'dorm_b0') + get_par(pars, 'dorm_b1') * x +
    get_par(pars, 'dorm_b2') * x^2 + get_par(pars, 'dorm_b3') * x^3
  inv_logit(eta)}

fl_x <- function(x, pars) {
  eta <- get_par(pars, 'fl_b0') + get_par(pars, 'fl_b1') * x +
    get_par(pars, 'fl_b2') * x^2 + get_par(pars, 'fl_b3') * x^3
  inv_logit(eta)}

se_x <- function(x, pars) {
  eta <- get_par(pars, 'se_b0') + get_par(pars, 'se_b1') * x +
    get_par(pars, 'se_b2') * x^2 + get_par(pars, 'se_b3') * x^3
  inv_logit(eta)}

se_n_x <- function(x, pars) {
  eta <- get_par(pars, 'se_n_b0') + get_par(pars, 'se_n_b1') * x +
    get_par(pars, 'se_n_b2') * x^2 + get_par(pars, 'se_n_b3') * x^3
  exp(eta)}

seedx <- function(x, pars) {
  fl_x(x, pars) * se_x(x, pars) * se_n_x(x, pars)}

re_per_seed <- function(pars) {
  exp(get_par(pars, 're_b0'))}

re_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
  if (!is.null(h)) dens <- dens / sum(dens * h)
  dens}

ra_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(y, mean = pars$react_sz, sd = pars$react_sd)
  if (!is.null(h)) dens <- dens / sum(dens * h)
  dens}

pxy <- function(x, y, pars) {
  sx(x, pars) * (1 - dx(x, pars)) * gxy(x, y, pars)}


# Kernel ----------------------------------------------------------------------
# State order: active size-bin counts, dormant count, delayed seed cohort.
# Rows = t+1; columns = t. Plants t -> seed cohort t+1 -> recruits t+2.
# The seed state carries production equivalents; re_per_seed includes losses.
kernel <- function(pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  i_dorm <- n + 1
  i_seed <- n + 2
  Smat <- sx(y, pars)
  Dvec <- dx(y, pars)
  Gmat <- t(outer(y, y, gxy, pars)) * h

  for (i in seq_len(n)) {
    if (i <= n / 2) {
      Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    } else {
      Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])}}

  Tmat <- sweep(Gmat, 2, Smat * (1 - Dvec), '*')
  react_y <- ra_y_dist(y, pars, h = h) * h
  recr_y <- re_y_dist(y, pars, h = h) * h

  Pmat <- matrix(0, n + 2, n + 2)
  Fmat <- matrix(0, n + 2, n + 2)

  Pmat[1:n, 1:n] <- Tmat
  Pmat[i_dorm, 1:n] <- Smat * Dvec
  Pmat[1:n, i_dorm] <-
    pars$p_dorm_survival * pars$p_reactivate * react_y
  Pmat[i_dorm, i_dorm] <-
    pars$p_dorm_survival * (1 - pars$p_reactivate)

  Fmat[i_seed, 1:n] <- seedx(y, pars)
  Pmat[1:n, i_seed] <- re_per_seed(pars) * recr_y

  k_yx <- Pmat + Fmat

  list(
    k_yx = k_yx, Pmat = Pmat, Fmat = Fmat, Tmat = Tmat,
    Gmat = Gmat, Smat = Smat, Dvec = Dvec, meshpts = y,
    dormant_index = i_dorm, seed_index = i_seed, h = h, L = L, U = U)}

lambda_ipm <- function(pars) {
  Re(eigen(kernel(pars)$k_yx)$values[1])}


# Mean IPM --------------------------------------------------------------------
lam_mean <- lambda_ipm(pars_mean)
lam_mean


# Year-specific IPMs ----------------------------------------------------------
lambda_ipm_year <- function(year, recruit_year = NULL) {
  pars_i <- make_ipm_pars(
    pars_mean = pars_all_mean, pars_year = pars_all_year,
    pars_re_year = pars_all_re_year, year = year,
    recruit_year = recruit_year)
  lambda_ipm(pars_i)}

# Retain years represented in all fitted yearly processes.
ipm_years <- sort(Reduce(intersect, list(
  as.integer(levels(df_su$year)),
  as.integer(levels(df_gr$year)),
  as.integer(levels(df_do$year)),
  as.integer(levels(df_dorm_su$year)),
  as.integer(levels(df_ra$year)),
  as.integer(levels(df_fl$year)),
  as.integer(levels(df_se$year)),
  as.integer(levels(df_se_n$year)),
  as.integer(levels(df_re_mod$recruit_year)) - 1L)))

lambda_year <- data.frame(
  year = ipm_years,
  lambda = sapply(ipm_years, lambda_ipm_year))

fig_lambda_year <- ggplot(lambda_year, aes(year, lambda)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(
    title = 'Year-specific asymptotic lambda', subtitle = v_ggp_suffix,
    x = 'Year', y = expression(lambda))

fig_lambda_year


# Observed and projected population growth ------------------------------------
# Match populations at both censuses within the reliable demographic period.
# Plant counts include active plants and dormant plants, but not seeds.
last_reliable_t1 <- max(
  df$year_t1[df$reliable_demography %in% TRUE], na.rm = TRUE)


# Population-year abundance ---------------------------------------------------
# As in POLE, include sized active plants and unsized recruits; add dormancy.
df_counts_population_year <- df %>%
  filter(year <= last_reliable_t1) %>%
  group_by(population, year) %>%
  summarise(
    n_sized = sum(dormant_t0 == 0 & is.finite(logsize_t0), na.rm = TRUE),
    n_unsized_recruits = sum(
      dormant_t0 == 0 & !is.finite(logsize_t0) & recruit == 1,
      na.rm = TRUE),
    n_dormant = sum(dormant_t0 == 1, na.rm = TRUE),
    n_ipm_state = n_sized + n_unsized_recruits + n_dormant,
    .groups = 'drop')


# Matched consecutive population transitions ----------------------------------
df_obs_pgr_population <- df_counts_population_year %>%
  arrange(population, year) %>%
  group_by(population) %>%
  mutate(
    year_t1 = lead(year), n_t1 = lead(n_ipm_state),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1, n_ipm_state > 0, year %in% ipm_years) %>%
  rename(n_t0 = n_ipm_state) %>%
  inner_join(
    df_re_seed %>%
      filter(nr_seed_obs > 0) %>%
      transmute(population, year = year + 1L, seeds_prev = nr_seeds),
    by = c('population', 'year')) %>%
  mutate(obs_pgr = n_t1 / n_t0)


# Initial population size distribution ----------------------------------------
make_initial_n_population <- function(year0, population_i, pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  breaks <- seq(L, U, length.out = n + 1)
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  df0 <- df %>% filter(year == year0, population == population_i)
  sizes <- df0 %>%
    filter(dormant_t0 == 0, is.finite(logsize_t0)) %>%
    pull(logsize_t0)

  size_counts <- hist(
    pmin(pmax(sizes, L), U), breaks = breaks,
    plot = FALSE, include.lowest = TRUE)$counts

  n_unsized_recruits <- sum(
    df0$dormant_t0 == 0 & !is.finite(df0$logsize_t0) & df0$recruit == 1,
    na.rm = TRUE)
  size_counts <- size_counts +
    n_unsized_recruits * re_y_dist(y, pars, h = h) * h

  n_dormant <- sum(df0$dormant_t0 == 1, na.rm = TRUE)
  n_seed <- df_re_seed %>%
    filter(population == population_i, year == year0 - 1L) %>%
    pull(nr_seeds)

  c(size_counts, n_dormant, n_seed)}


# Project one population-year -------------------------------------------------
# The initial seed cohort uses observed production in t-1.
# Seeds predicted during t affect plant abundance at t+2, not t+1.
project_one_population_year <- function(yr, population_i) {
  pars_y <- make_ipm_pars(
    pars_mean = pars_all_mean, pars_year = pars_all_year,
    pars_re_year = pars_all_re_year, year = yr, recruit_year = yr + 1L)

  n_obs <- make_initial_n_population(
    year0 = yr, population_i = population_i, pars = pars_y)
  K <- kernel(pars_y)$k_yx
  n_proj <- K %*% n_obs

  plant_indices <- seq_len(pars_y$mat_siz + 1)
  n_initial <- sum(n_obs[plant_indices])
  n_projected <- sum(n_proj[plant_indices])

  data.frame(
    year = yr, population = population_i,
    n_obs_model = n_initial, n_proj_model = n_projected,
    asym_lambda = Re(eigen(K)$values[1]),
    proj_lambda = as.numeric(n_projected / n_initial))}


# Project all matched population transitions ----------------------------------
df_proj_population <- bind_rows(lapply(
  seq_len(nrow(df_obs_pgr_population)), function(i) {
    project_one_population_year(
      yr = df_obs_pgr_population$year[i],
      population_i = df_obs_pgr_population$population[i])}))


# Population-level observed and modeled growth ---------------------------------
df_compare_population <- df_obs_pgr_population %>%
  left_join(df_proj_population, by = c('year', 'population')) %>%
  mutate(
    error_asymptotic_vs_obs = asym_lambda - obs_pgr,
    error_projected_vs_obs = proj_lambda - obs_pgr)


# Whole-population annual comparison ------------------------------------------
df_compare <- df_compare_population %>%
  group_by(year) %>%
  summarise(
    asym_lambda = weighted.mean(asym_lambda, w = n_obs_model),
    n_t0 = sum(n_t0), n_t1 = sum(n_t1),
    n_obs_model = sum(n_obs_model), n_proj_model = sum(n_proj_model),
    n_populations = n(), .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 / n_t0,
    proj_lambda = n_proj_model / n_obs_model)


# Observed vs modeled plot ----------------------------------------------------
df_plot <- df_compare %>%
  dplyr::select(year, obs_pgr, asym_lambda, proj_lambda) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = 'lambda_type', values_to = 'lambda') %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = 'Year-specific asymptotic lambda',
      proj_lambda = 'One-step lambda from observed population structure'))

fig_mod_vs_obs <- ggplot(df_plot, aes(x = lambda, y = obs_pgr)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = 'free_x') +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = v_ggp_suffix,
    x = expression('Modeled ' * lambda),
    y = 'Observed population growth rate') +
  theme_classic()

fig_mod_vs_obs


# Log-transformed observed vs modeled plot -------------------------------------
df_plot_log <- df_plot %>%
  filter(obs_pgr > 0, lambda > 0) %>%
  mutate(log_obs_pgr = log(obs_pgr), log_lambda = log(lambda))

fig_mod_vs_obs_log <- ggplot(
  df_plot_log, aes(x = log_lambda, y = log_obs_pgr)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = 'free_x') +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = paste(v_ggp_suffix, '- log-transformed lambda'),
    x = expression('log modeled ' * lambda),
    y = 'log observed population growth rate') +
  theme_classic()

fig_mod_vs_obs_log


# Summary statistics ----------------------------------------------------------
df_compare_summary <- df_compare %>%
  summarise(
    n_years = n(),
    arithmetic_mean_obs_pgr = mean(obs_pgr, na.rm = TRUE),
    geometric_mean_obs_pgr = exp(mean(log(obs_pgr), na.rm = TRUE)),
    arithmetic_mean_asym_lambda = mean(asym_lambda, na.rm = TRUE),
    geometric_mean_asym_lambda = exp(mean(log(asym_lambda), na.rm = TRUE)),
    arithmetic_mean_proj_lambda = mean(proj_lambda, na.rm = TRUE),
    geometric_mean_proj_lambda = exp(mean(log(proj_lambda), na.rm = TRUE)),
    mean_error_asymptotic_vs_obs = mean(asym_lambda - obs_pgr, na.rm = TRUE),
    mean_error_projected_vs_obs = mean(proj_lambda - obs_pgr, na.rm = TRUE),
    percent_bias_asymptotic_vs_obs = 100 *
      sum(asym_lambda - obs_pgr, na.rm = TRUE) / sum(obs_pgr, na.rm = TRUE),
    percent_bias_projected_vs_obs = 100 *
      sum(proj_lambda - obs_pgr, na.rm = TRUE) / sum(obs_pgr, na.rm = TRUE),
    rmse_asymptotic_vs_obs = sqrt(mean(
      (asym_lambda - obs_pgr)^2, na.rm = TRUE)),
    rmse_projected_vs_obs = sqrt(mean(
      (proj_lambda - obs_pgr)^2, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(), names_to = 'statistic', values_to = 'value')

df_compare_summary
df_compare


# Population-level summary ----------------------------------------------------
df_compare_population_summary <- df_compare_population %>%
  filter(
    !is.na(obs_pgr), !is.na(asym_lambda), !is.na(proj_lambda),
    obs_pgr > 0, asym_lambda > 0, proj_lambda > 0) %>%
  group_by(population) %>%
  summarise(
    n_year_transitions = n(),
    lambda_obs_geometric = exp(mean(log(obs_pgr), na.rm = TRUE)),
    lambda_obs_arithmetic = mean(obs_pgr, na.rm = TRUE),
    lambda_asymptotic_geometric = exp(mean(log(asym_lambda), na.rm = TRUE)),
    lambda_asymptotic_arithmetic = mean(asym_lambda, na.rm = TRUE),
    lambda_projected_geometric = exp(mean(log(proj_lambda), na.rm = TRUE)),
    lambda_projected_arithmetic = mean(proj_lambda, na.rm = TRUE),
    error_projected_geo_vs_obs_geo =
      lambda_projected_geometric - lambda_obs_geometric,
    rmse_projected_vs_obs = sqrt(mean((proj_lambda - obs_pgr)^2,
                                      na.rm = TRUE)),
    mean_n_initial = mean(n_obs_model, na.rm = TRUE), .groups = 'drop') %>%
  arrange(population)

df_compare_population_summary


# Save key outputs ------------------------------------------------------------
# saveRDS(pars_all_mean,
#         file.path(dir_result, 'mcgraw_paqu_yr_pars_mean.rds'))
# saveRDS(pars_all_year,
#         file.path(dir_result, 'mcgraw_paqu_yr_pars_year.rds'))
# saveRDS(pars_all_re_year,
#         file.path(dir_result, 'mcgraw_paqu_yr_pars_recruit_year.rds'))
# saveRDS(lambda_year,
#         file.path(dir_result, 'mcgraw_paqu_yr_lambda.rds'))
# saveRDS(df_compare,
#         file.path(dir_result, 'mcgraw_paqu_yr_lambda_compare.rds'))
# saveRDS(df_compare_population,
#         file.path(dir_result, 'mcgraw_paqu_yr_lambda_compare_population.rds'))
