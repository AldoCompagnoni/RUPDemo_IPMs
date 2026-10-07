# IPM year-specific; zero = most basic
# McGraw 2017 - Panax quinquefolius

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.10.07

# Study organism: Panax quinquefolius L.
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.9.4
# Meta data link:
# https://portal.edirepository.org/nis/metadataviewer?packageid=edi.9.4
# Citing publication: McGraw et al. 2017. Long Term Research in
# Environmental Biology: Demographic census data for thirty natural
# populations of American Ginseng: 1998-2016 ver 4.
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
  GGally,
  brms,
  loo)


# Specification ---------------------------------------------------------------
v_brm_suffix <- ""
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

# Keep manual model choices empty for AICc selection.
# Values 0:3 restrict survival/growth to the corresponding polynomial degree.
v_mod_set_su <- c()
v_mod_set_gr <- c()


# Directory -------------------------------------------------------------------
dir_wd <- file.path("C:/code/RUPDemo_IPMs")
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
source(file.path(dir_wd, 'helper_functions/line_color_pred_fun.R'))
source(file.path(dir_wd, 'helper_functions/predictor_fun.R'))


# Data ------------------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0(v_head, '_', v_sp_abb, '_df_workdata.csv'))) %>%
  mutate(
    year = as.integer(year)) %>%
  filter(
    state == v_states,
    !is.na(year),
    !(year %in% v_years_re))


# Controls --------------------------------------------------------------------
ctrl_glmer <- glmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

ctrl_lmer <- lmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))


# Survival --------------------------------------------------------------------
df_su <- df %>%
  filter(
    reliable_demography,
    !is.na(survives),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  mutate(
    year = factor(year)) %>%
  dplyr::select(
    state, population, id, year, size_t0, survives,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Survival model --------------------------------------------------------------
mod_su_0 <- glmer(
  survives ~ 1 + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer)

mod_su_1 <- glmer(
  survives ~ logsize_t0 + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer)

mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer)

mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (1 | year),
  data = df_su,
  family = binomial,
  control = ctrl_glmer)

mods_su <- list(
  mod_su_0, mod_su_1, mod_su_2, mod_su_3)

mods_su_dAICc <- bbmle::AICctab(
  mods_su,
  weights = TRUE,
  sort = FALSE)$dAICc

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
  df_i <- df_su %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE),
    length.out = 100)

  pred_i <- tibble(
    logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(
        year_i,
        levels = levels(df_su$year)))

  pred_i$survives <- predict(
    mod_su_best,
    newdata = pred_i,
    type = 'response',
    re.form = NULL)

  pts_i <- df_i %>%
    mutate(
      bin = cut(
        logsize_t0,
        breaks = 8)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      survives = mean(survives, na.rm = TRUE),
      n = n(),
      .groups = 'drop') %>%
    filter(
      !is.na(logsize_t0),
      n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, survives),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, survives),
      linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = expression('Survival probability'[t1])) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      text = element_text(size = 5),
      legend.position = 'none')
}

su_yrs <- lapply(
  levels(df_su$year),
  make_su_year_plot)

fig_su_years <- wrap_plots(su_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Survival - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(
        size = 13,
        face = 'bold'),
      plot.subtitle = element_text(
        size = 9)))

fig_su_years


# Growth data -----------------------------------------------------------------
df_gr <- df %>%
  filter(
    consecutive,
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1)) %>%
  mutate(
    year = factor(year)) %>%
  dplyr::select(
    state, population, id, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

fig_gr_raw <- ggplot(
  df_gr,
  aes(logsize_t0, logsize_t1)) +
  geom_point(
    alpha = 0.5,
    pch = 16,
    size = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1) +
  theme_bw() +
  labs(
    title = 'Growth',
    subtitle = v_ggp_suffix,
    x = expression('log(leaf area)'[t0]),
    y = expression('log(leaf area)'[t1])) +
  theme(
    plot.subtitle = element_text(size = 8))

fig_gr_raw


# Growth model ----------------------------------------------------------------
mod_gr_0 <- lmer(
  logsize_t1 ~ 1 + (logsize_t0 | year),
  data = df_gr,
  REML = FALSE,
  control = ctrl_lmer)

mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + (logsize_t0 | year),
  data = df_gr,
  REML = FALSE,
  control = ctrl_lmer)

mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 +
    (logsize_t0 | year),
  data = df_gr,
  REML = FALSE,
  control = ctrl_lmer)

mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
    (logsize_t0 | year),
  data = df_gr,
  REML = FALSE,
  control = ctrl_lmer)

mods_gr <- list(
  mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)

mods_gr_dAICc <- bbmle::AICctab(
  mods_gr,
  weights = TRUE,
  sort = FALSE)$dAICc

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
  df_i <- df_gr %>%
    filter(year == year_i)

  x <- seq(
    min(df_i$logsize_t0, na.rm = TRUE),
    max(df_i$logsize_t0, na.rm = TRUE),
    length.out = 100)

  pred_i <- tibble(
    logsize_t0 = x) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      year = factor(
        year_i,
        levels = levels(df_gr$year)))

  pred_i$logsize_t1 <- predict(
    mod_gr_best,
    newdata = pred_i,
    re.form = NULL)

  pts_i <- df_i %>%
    mutate(
      bin = cut(
        logsize_t0,
        breaks = 8)) %>%
    group_by(bin) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      logsize_t1 = mean(logsize_t1, na.rm = TRUE),
      n = n(),
      .groups = 'drop') %>%
    filter(
      !is.na(logsize_t0),
      n > 0)

  ggplot() +
    geom_point(
      data = pts_i,
      aes(logsize_t0, logsize_t1),
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(logsize_t0, logsize_t1),
      linewidth = 0.7) +
    geom_abline(
      intercept = 0,
      slope = 1,
      lty = 2) +
    labs(
      title = year_i,
      x = expression('log(leaf area)'[t0]),
      y = expression('log(leaf area)'[t1])) +
    theme_bw() +
    theme(
      text = element_text(size = 5),
      legend.position = 'none')
}

gr_yrs <- lapply(
  levels(df_gr$year),
  make_gr_year_plot)

fig_gr_years <- wrap_plots(gr_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Growth - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(
        size = 13,
        face = 'bold'),
      plot.subtitle = element_text(
        size = 9)))

fig_gr_years


# Growth variance -------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_best)
mod_gr_y <- resid(mod_gr_best)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(
    maxiter = 1000,
    tol = 1e-6,
    warnOnly = TRUE))


# Recruitment -----------------------------------------------------------------
# Constant recruitment per established plant with a recruit-year random effect.

df_re_parent <- df %>%
  filter(
    !is.na(persistence_t0),
    persistence_t0 != 'DEAD',
    size_t0 > 0) %>%
  group_by(
    population,
    year) %>%
  summarise(
    n_parents = n_distinct(id),
    .groups = 'drop') %>%
  mutate(
    year_t1 = year + 1L)

df_re_sampled <- df %>%
  distinct(
    population,
    year) %>%
  rename(
    year_t1 = year)

df_re_count <- df %>%
  filter(
    recruit == 1) %>%
  count(
    population,
    year,
    name = 'nr_recruits') %>%
  rename(
    year_t1 = year)

df_re <- df_re_parent %>%
  inner_join(
    df_re_sampled,
    by = c(
      'population',
      'year_t1')) %>%
  left_join(
    df_re_count,
    by = c(
      'population',
      'year_t1')) %>%
  mutate(
    nr_recruits = replace_na(
      nr_recruits,
      0L),
    year_t1 = factor(year_t1))

mod_re <- glmer.nb(
  nr_recruits ~ 1 + offset(log(n_parents)) +
    (1 | year_t1),
  data = df_re,
  control = ctrl_glmer)

mod_re_ranef <- coef(mod_re)$year_t1

mod_re
summary(mod_re)

fecu_mean <- exp(
  unname(
    fixef(mod_re)[['(Intercept)']]))


# Recruit size distribution ---------------------------------------------------
df_re_size <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0))


# Recruitment plots by year ---------------------------------------------------
make_re_year_plot <- function(year_i) {
  df_i <- df_re %>%
    filter(year_t1 == year_i)

  x <- seq(
    min(df_i$n_parents, na.rm = TRUE),
    max(df_i$n_parents, na.rm = TRUE),
    length.out = 100)

  pred_i <- tibble(
    n_parents = x,
    year_t1 = factor(
      year_i,
      levels = levels(df_re$year_t1)))

  pred_i$nr_recruits <- predict(
    mod_re,
    newdata = pred_i,
    type = 'response',
    re.form = NULL)

  ggplot() +
    geom_point(
      data = df_i,
      aes(n_parents, nr_recruits),
      alpha = 0.6,
      size = 1.1) +
    geom_line(
      data = pred_i,
      aes(n_parents, nr_recruits),
      linewidth = 0.7) +
    labs(
      title = year_i,
      x = expression('Established individuals'[t0]),
      y = expression('New seedlings'[t1])) +
    theme_bw() +
    theme(
      text = element_text(size = 5),
      legend.position = 'none')
}

re_yrs <- lapply(
  levels(df_re$year_t1),
  make_re_year_plot)

fig_re_years <- wrap_plots(re_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = 'Recruitment - year specific',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(
        size = 13,
        face = 'bold'),
      plot.subtitle = element_text(
        size = 9)))

fig_re_years


# Parameter extraction --------------------------------------------------------
get_coef <- function(x, name) {
  if (name %in% names(x)) {
    unname(x[[name]])
  } else {
    0
  }
}


# Mean parameters -------------------------------------------------------------
su_fix <- fixef(mod_su_best)
gr_fix <- fixef(mod_gr_best)
gr_var_coef <- coef(mod_gr_var)

pars_mean <- list(
  prefix = v_script_prefix,
  species = v_species,

  surv_b0 = get_coef(
    su_fix,
    '(Intercept)'),
  surv_b1 = get_coef(
    su_fix,
    'logsize_t0'),
  surv_b2 = get_coef(
    su_fix,
    'logsize_t0_2'),
  surv_b3 = get_coef(
    su_fix,
    'logsize_t0_3'),

  grow_b0 = get_coef(
    gr_fix,
    '(Intercept)'),
  grow_b1 = get_coef(
    gr_fix,
    'logsize_t0'),
  grow_b2 = get_coef(
    gr_fix,
    'logsize_t0_2'),
  grow_b3 = get_coef(
    gr_fix,
    'logsize_t0_3'),

  a = as.numeric(
    gr_var_coef[1]),
  b = as.numeric(
    gr_var_coef[2]),

  fecu_b0 = fecu_mean,

  recr_sz = mean(
    df_re_size$logsize_t0,
    na.rm = TRUE),
  recr_sd = sd(
    df_re_size$logsize_t0,
    na.rm = TRUE),

  L = min(
    df_gr$logsize_t0,
    na.rm = TRUE),
  U = max(
    df_gr$logsize_t0,
    na.rm = TRUE),
  mat_siz = 200,

  mod_gr_index = v_mod_gr_index,
  mod_su_index = v_mod_su_index)


# Year-specific parameters ----------------------------------------------------
make_ipm_pars <- function(
    year,
    re_year_t1 = NULL) {

  pars <- pars_mean
  year <- as.character(year)

  if (year %in% rownames(mod_su_ranef)) {
    pars$surv_b0 <- get_coef(
      mod_su_ranef[year, ],
      '(Intercept)')
    pars$surv_b1 <- get_coef(
      mod_su_ranef[year, ],
      'logsize_t0')
    pars$surv_b2 <- get_coef(
      mod_su_ranef[year, ],
      'logsize_t0_2')
    pars$surv_b3 <- get_coef(
      mod_su_ranef[year, ],
      'logsize_t0_3')
  }

  if (year %in% rownames(mod_gr_ranef)) {
    pars$grow_b0 <- get_coef(
      mod_gr_ranef[year, ],
      '(Intercept)')
    pars$grow_b1 <- get_coef(
      mod_gr_ranef[year, ],
      'logsize_t0')
    pars$grow_b2 <- get_coef(
      mod_gr_ranef[year, ],
      'logsize_t0_2')
    pars$grow_b3 <- get_coef(
      mod_gr_ranef[year, ],
      'logsize_t0_3')
  }

  if (is.null(re_year_t1)) {
    re_year_t1 <- as.integer(year) + 1L
  }

  re_year_t1 <- as.character(re_year_t1)

  if (re_year_t1 %in% rownames(mod_re_ranef)) {
    pars$fecu_b0 <- exp(
      get_coef(
        mod_re_ranef[re_year_t1, ],
        '(Intercept)'))
  }

  pars
}


# IPM functions ---------------------------------------------------------------
grow_sd <- function(x, pars) {
  sqrt(
    pars$a *
      exp(pars$b * x))
}


# Growth from size x to size y.
gxy <- function(
    x,
    y,
    pars,
    num_pars = v_mod_gr_index) {

  mean_value <- 0

  for (i in 0:num_pars) {
    param_name <- paste0(
      'grow_b',
      i)

    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value +
        pars[[param_name]] *
        x^i
    }
  }

  sd_value <- grow_sd(
    x,
    pars)

  dnorm(
    y,
    mean = mean_value,
    sd = sd_value)
}


inv_logit <- function(x) {
  exp(x) /
    (1 + exp(x))
}


# Survival of an x-sized individual to t1.
sx <- function(
    x,
    pars,
    num_pars = v_mod_su_index) {

  survival_value <- pars$surv_b0

  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0(
        'surv_b',
        i)

      if (!is.null(pars[[param_name]])) {
        survival_value <- survival_value +
          pars[[param_name]] *
          x^i
      }
    }
  }

  inv_logit(
    survival_value)
}


# Survival-growth transition.
pxy <- function(
    x,
    y,
    pars) {

  sx(
    x,
    pars) *
    gxy(
      x,
      y,
      pars)
}


# Constant fecundity distributed across recruit sizes.
fy <- function(
    y,
    pars,
    h) {

  n_recr <- pars$fecu_b0

  recr_sd <- max(
    h / 10,
    pars$recr_sd,
    na.rm = TRUE)

  recr_y <- dnorm(
    y,
    pars$recr_sz,
    recr_sd) * h

  recr_y <- recr_y /
    sum(recr_y)

  n_recr *
    recr_y
}


# Kernel ----------------------------------------------------------------------
kernel <- function(pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  y <- 0.5 *
    (b[1:n] + b[2:(n + 1)])

  Fmat <- matrix(
    0,
    n,
    n)

  Fmat[] <- matrix(
    fy(
      y,
      pars,
      h),
    n,
    n)

  Smat <- sx(
    y,
    pars)

  Gmat <- matrix(
    0,
    n,
    n)

  Gmat[] <- t(
    outer(
      y,
      y,
      gxy,
      pars)) * h

  Tmat <- matrix(
    0,
    n,
    n)

  for (i in seq_len(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] +
      1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] *
      Smat[i]
  }

  for (i in ((n / 2) + 1):n) {
    Gmat[n, i] <- Gmat[n, i] +
      1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] *
      Smat[i]
  }

  k_yx <- Fmat +
    Tmat

  list(
    k_yx = k_yx,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    meshpts = y,
    h = h,
    L = L,
    U = U)
}


lambda_ipm <- function(pars) {
  Re(
    eigen(
      kernel(pars)$k_yx)$values[1])
}


# Mean IPM --------------------------------------------------------------------
lam_mean <- lambda_ipm(
  pars_mean)

lam_mean


# Year-specific IPMs ----------------------------------------------------------
lambda_ipm_year <- function(
    year,
    re_year_t1 = NULL) {

  pars_i <- make_ipm_pars(
    year = year,
    re_year_t1 = re_year_t1)

  lambda_ipm(
    pars_i)
}

ipm_years <- sort(
  Reduce(
    intersect,
    list(
      as.integer(levels(df_su$year)),
      as.integer(levels(df_gr$year)),
      as.integer(levels(df_re$year_t1)) - 1L)))

lambda_year <- data.frame(
  year = ipm_years,
  lambda = sapply(
    ipm_years,
    lambda_ipm_year))

fig_lambda_year <- ggplot(
  lambda_year,
  aes(year, lambda)) +
  geom_hline(
    yintercept = 1,
    linetype = 'dashed') +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(
    title = 'Year-specific asymptotic lambda',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = expression(lambda))

fig_lambda_year


# Observed and projected population growth ------------------------------------
# Population-year abundance ---------------------------------------------------
df_counts_population_year <- df %>%
  group_by(
    population,
    year) %>%
  summarise(
    n_sized = sum(
      is.finite(logsize_t0),
      na.rm = TRUE),
    n_unsized_recruits = sum(
      !is.finite(logsize_t0) &
        recruit == 1,
      na.rm = TRUE),
    n_ipm_state = n_sized +
      n_unsized_recruits,
    .groups = 'drop')


# Matched consecutive population transitions ----------------------------------
df_obs_pgr_population <- df_counts_population_year %>%
  arrange(
    population,
    year) %>%
  group_by(
    population) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n_ipm_state),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(
    year_gap == 1,
    n_ipm_state > 0,
    year %in% ipm_years) %>%
  rename(
    n_t0 = n_ipm_state) %>%
  mutate(
    obs_pgr = n_t1 / n_t0)


# Initial population size distribution ----------------------------------------
make_initial_n_population <- function(
    year0,
    population_i,
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
    (b[1:n] +
       b[2:(n + 1)])

  df0 <- df %>%
    filter(
      year == year0,
      population == population_i)

  sizes <- df0 %>%
    filter(
      is.finite(logsize_t0)) %>%
    pull(
      logsize_t0)

  size_counts <- hist(
    pmin(
      pmax(
        sizes,
        L),
      U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE)$counts

  n_density <- size_counts /
    h

  n_unsized_recruits <- df0 %>%
    summarise(
      n = sum(
        !is.finite(logsize_t0) &
          recruit == 1,
        na.rm = TRUE)) %>%
    pull(n)

  if (n_unsized_recruits > 0) {
    recr_y <- dnorm(
      y,
      mean = pars$recr_sz,
      sd = pars$recr_sd)

    recr_y <- recr_y /
      sum(recr_y * h)

    n_density <- n_density +
      n_unsized_recruits *
      recr_y
  }

  n_density
}


# Project one population-year -------------------------------------------------
project_one_population_year <- function(
    yr,
    population_i) {

  pars_y <- make_ipm_pars(
    year = yr,
    re_year_t1 = yr + 1L)

  n_obs <- make_initial_n_population(
    year0 = yr,
    population_i = population_i,
    pars = pars_y)

  K <- kernel(
    pars_y)$k_yx

  h <- (pars_y$U - pars_y$L) /
    pars_y$mat_siz

  n_initial <- sum(n_obs) *
    h

  n_proj <- K %*%
    n_obs

  n_projected <- sum(n_proj) *
    h

  data.frame(
    year = yr,
    population = population_i,
    n_obs_model = n_initial,
    n_proj_model = n_projected,
    asym_lambda = Re(
      eigen(K)$values[1]),
    proj_lambda = as.numeric(
      n_projected /
        n_initial))
}


# Project all matched population transitions ----------------------------------
df_proj_population <- bind_rows(
  lapply(
    seq_len(nrow(df_obs_pgr_population)),
    function(i) {

      project_one_population_year(
        yr = df_obs_pgr_population$year[i],
        population_i =
          df_obs_pgr_population$population[i])
    }))


# Population-level observed and modeled growth --------------------------------
df_compare_population <- df_obs_pgr_population %>%
  left_join(
    df_proj_population,
    by = c(
      'year',
      'population')) %>%
  mutate(
    error_asymptotic_vs_obs =
      asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr)


# Whole-population annual comparison ------------------------------------------
df_compare <- df_compare_population %>%
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
    n_populations = n(),
    .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 /
      n_t0,
    proj_lambda = n_proj_model /
      n_obs_model)


# Observed vs modeled plot ----------------------------------------------------
df_plot <- df_compare %>%
  dplyr::select(
    year,
    obs_pgr,
    asym_lambda,
    proj_lambda) %>%
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
        'Year-specific asymptotic lambda',
      proj_lambda =
        'Projected lambda from observed population size distributions'))

fig_mod_vs_obs <- ggplot(
  df_plot,
  aes(
    x = lambda,
    y = obs_pgr)) +
  geom_point(
    size = 3) +
  geom_abline(
    intercept = 0,
    slope = 1,
    lty = 2) +
  facet_wrap(
    ~ lambda_type,
    scales = 'free_x') +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = v_ggp_suffix,
    x = expression('Modeled ' * lambda),
    y = 'Observed population growth rate') +
  theme_classic()

fig_mod_vs_obs


# Log-transformed observed vs modeled plot ------------------------------------
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
    y = log_obs_pgr)) +
  geom_point(
    size = 3) +
  geom_abline(
    intercept = 0,
    slope = 1,
    lty = 2) +
  facet_wrap(
    ~ lambda_type,
    scales = 'free_x') +
  labs(
    title = 'Observed population growth vs modeled lambda',
    subtitle = paste(
      v_ggp_suffix,
      '- log-transformed lambda'),
    x = expression('log modeled ' * lambda),
    y = 'log observed population growth rate') +
  theme_classic()

fig_mod_vs_obs_log


# Summary statistics ----------------------------------------------------------
df_compare_summary <- df_compare %>%
  summarise(
    n_years = n(),
    arithmetic_mean_obs_pgr = mean(
      obs_pgr,
      na.rm = TRUE),
    geometric_mean_obs_pgr = exp(
      mean(
        log(obs_pgr),
        na.rm = TRUE)),
    arithmetic_mean_asym_lambda = mean(
      asym_lambda,
      na.rm = TRUE),
    geometric_mean_asym_lambda = exp(
      mean(
        log(asym_lambda),
        na.rm = TRUE)),
    arithmetic_mean_proj_lambda = mean(
      proj_lambda,
      na.rm = TRUE),
    geometric_mean_proj_lambda = exp(
      mean(
        log(proj_lambda),
        na.rm = TRUE)),
    mean_error_asymptotic_vs_obs = mean(
      asym_lambda - obs_pgr,
      na.rm = TRUE),
    mean_error_projected_vs_obs = mean(
      proj_lambda - obs_pgr,
      na.rm = TRUE),
    percent_bias_asymptotic_vs_obs = 100 *
      sum(
        asym_lambda - obs_pgr,
        na.rm = TRUE) /
      sum(
        obs_pgr,
        na.rm = TRUE),
    percent_bias_projected_vs_obs = 100 *
      sum(
        proj_lambda - obs_pgr,
        na.rm = TRUE) /
      sum(
        obs_pgr,
        na.rm = TRUE),
    rmse_asymptotic_vs_obs = sqrt(
      mean(
        (asym_lambda - obs_pgr)^2,
        na.rm = TRUE)),
    rmse_projected_vs_obs = sqrt(
      mean(
        (proj_lambda - obs_pgr)^2,
        na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'statistic',
    values_to = 'value')

df_compare_summary

df_compare %>%
  print(n = 100)


# Population-level summary ----------------------------------------------------
df_compare_population_summary <- df_compare_population %>%
  filter(
    !is.na(obs_pgr),
    !is.na(asym_lambda),
    !is.na(proj_lambda),
    obs_pgr > 0,
    asym_lambda > 0,
    proj_lambda > 0) %>%
  group_by(population) %>%
  summarise(
    n_year_transitions = n(),
    lambda_obs_geometric = exp(
      mean(
        log(obs_pgr),
        na.rm = TRUE)),
    lambda_obs_arithmetic = mean(
      obs_pgr,
      na.rm = TRUE),
    lambda_asymptotic_geometric = exp(
      mean(
        log(asym_lambda),
        na.rm = TRUE)),
    lambda_asymptotic_arithmetic = mean(
      asym_lambda,
      na.rm = TRUE),
    lambda_projected_geometric = exp(
      mean(
        log(proj_lambda),
        na.rm = TRUE)),
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
    .groups = 'drop') %>%
  arrange(population)

df_compare_population_summary %>%
  print(
    n = 100,
    width = Inf)


# # Save key outputs -----------------------------------------------------------
# # Uncomment if needed.
# # saveRDS(
# #   pars_mean,
# #   file.path(
# #     dir_result,
# #     'paqu_year_zero_pars_mean.rds'))
# # saveRDS(
# #   lambda_year,
# #   file.path(
# #     dir_result,
# #     'paqu_year_zero_lambda.rds'))
# # saveRDS(
# #   df_compare,
# #   file.path(
# #     dir_result,
# #     'paqu_year_zero_lambda_compare.rds'))
