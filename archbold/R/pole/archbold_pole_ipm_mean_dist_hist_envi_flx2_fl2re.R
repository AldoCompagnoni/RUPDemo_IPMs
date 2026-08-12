# Mean IPM with fire - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.30

# Study organism: Polygala lewtonii
# Citing publication: https://doi.org/10.1071/BT11271
# Time period: 2001-2024, through 2025

# This script fits mean vital-rate models with fire as a fixed effect.
# It contains no year random effects and no site random effects.
# Year and site are retained only for data preparation, recruitment timing,
# disturbance exposure, and observed-growth validation.

# rm(list = ls())


# Setting the stage ------------------------------------------------------------
set.seed(100)
options(stringsAsFactors = FALSE)


# Packages ---------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(
  tidyverse,
  patchwork,
  bbmle,
  scales)


# Specification ---------------------------------------------------------------
v_head <- 'archbold'
v_species <- 'Polygala lewtonii'
v_years_re <- c(2001)

v_sp_abb <- tolower(
  gsub(
    ' ', '',
    paste(
      substr(unlist(strsplit(v_species, ' ')), 1, 2),
      collapse = '')))

v_script_prefix <- v_head
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Empty vectors select the best-supported polynomial by AICc.
v_mod_set_su <- c()
v_mod_set_gr <- c()
v_mod_set_fl <- c()
v_mod_set_fl_n <- c()
v_mod_set_re <- c()


# Directory --------------------------------------------------------------------
dir_pub <- file.path(v_head)
dir_R <- file.path(dir_pub, 'R', v_sp_abb)
dir_data <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

for (dir_i in c(
  file.path(dir_pub, 'R'),
  file.path(dir_pub, 'data'),
  file.path(dir_pub, 'results'),
  dir_R, dir_data, dir_result)) {
  if (!dir.exists(dir_i)) {
    dir.create(dir_i, recursive = TRUE)
  }
}


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


# Vital-rate data --------------------------------------------------------------
df_su <- df %>%
  filter(
    !is.na(survives),
    size_t0 != 0,
    is.finite(logvol_t0),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance)) %>%
  dplyr::select(
    id, site, year, size_t0, survives, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3, disturbance, stage)

df_gr <- df %>%
  filter(
    size_t0 != 0,
    size_t1 != 0,
    is.finite(logvol_t0),
    is.finite(logvol_t1),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance)) %>%
  dplyr::select(
    id, site, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3, disturbance, stage)

# Recruits are retained because some observed recruits flowered.
df_fl <- df %>%
  filter(
    !is.na(flower),
    is.finite(logvol_t0),
    is.finite(logvol_t0_2),
    is.finite(logvol_t0_3),
    !is.na(disturbance_prev)) %>%
  dplyr::select(
    id, site, year, size_t0, flower, fl_nr, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1,
    logvol_t0_2, logvol_t0_3, disturbance_prev, stage, recruits)

df_fl_cond <- df_fl %>%
  filter(
    flower == 1,
    !is.na(fl_nr),
    fl_nr > 0,
    fl_nr == round(fl_nr))


# Flowering stems to recruits -----------------------------------------------
# Total flowering stems in quadrat q during year t predict recruits in the
# same quadrat during year t + 1. Current fire is included as a quadrat-level
# fixed effect.
df_quad_year_sampled <- df %>%
  distinct(site, quad, year)

df_fs2r <- df %>%
  group_by(site, quad, year) %>%
  summarise(
    fs_t0 = sum(fl_nr, na.rm = TRUE),
    disturbance = if_else(
      all(is.na(disturbance)),
      NA_real_,
      max(disturbance, na.rm = TRUE)),
    .groups = "drop") %>%
  mutate(
    year_t1 = year + 1) %>%
  semi_join(
    df_quad_year_sampled %>%
      transmute(
        site,
        quad,
        year_t1 = year),
    by = c(
      "site",
      "quad",
      "year_t1")) %>%
  left_join(
    df %>%
      filter(recruits == 1) %>%
      count(
        site,
        quad,
        year,
        name = "re_t1"),
    by = c(
      "site",
      "quad",
      "year_t1" = "year")) %>%
  mutate(
    re_t1 = replace_na(re_t1, 0L),
    site = factor(site),
    fs_t0_log = log1p(fs_t0),
    re_t1_log = log1p(re_t1)) %>%
  filter(
    !is.na(year),
    !is.na(year_t1),
    !is.na(disturbance),
    !(year %in% v_years_re),
    !(year_t1 %in% v_years_re))


# Survival model ---------------------------------------------------------------
mod_su_00 <- glm(
  survives ~ 1,
  data = df_su, family = binomial)

mod_su_0 <- glm(
  survives ~ disturbance,
  data = df_su, family = binomial)

mod_su_10 <- glm(
  survives ~ logvol_t0,
  data = df_su, family = binomial)

mod_su_1 <- glm(
  survives ~ logvol_t0 + disturbance,
  data = df_su, family = binomial)

mod_su_20 <- glm(
  survives ~ logvol_t0 + logvol_t0_2,
  data = df_su, family = binomial)

mod_su_2 <- glm(
  survives ~ logvol_t0 + logvol_t0_2 + disturbance,
  data = df_su, family = binomial)

mod_su_30 <- glm(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
  data = df_su, family = binomial)

mod_su_3 <- glm(
  survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance,
  data = df_su, family = binomial)

mods_su <- list(
  mod_su_00, mod_su_0, mod_su_10, mod_su_1,
  mod_su_20, mod_su_2, mod_su_30, mod_su_3)
mods_su_dAICc <- bbmle::AICctab(
  mods_su, weights = TRUE, sort = FALSE)$dAICc
mods_su_sorted <- order(mods_su_dAICc)

if (length(v_mod_set_su) == 0) {
  mod_su_index_best <- mods_su_sorted[1]
} else {
  mod_su_index_best <- v_mod_set_su + 1
}

v_mod_su_index <- floor((mod_su_index_best - 1) / 2)
mod_su_best <- mods_su[[mod_su_index_best]]
summary(mod_su_best)
mods_su_dAICc


# Growth model -----------------------------------------------------------------
mod_gr_00 <- lm(
  logvol_t1 ~ 1,
  data = df_gr)

mod_gr_0 <- lm(
  logvol_t1 ~ disturbance,
  data = df_gr)

mod_gr_10 <- lm(
  logvol_t1 ~ logvol_t0,
  data = df_gr)

mod_gr_1 <- lm(
  logvol_t1 ~ logvol_t0 + disturbance,
  data = df_gr)

mod_gr_20 <- lm(
  logvol_t1 ~ logvol_t0 + logvol_t0_2,
  data = df_gr)

mod_gr_2 <- lm(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + disturbance,
  data = df_gr)

mod_gr_30 <- lm(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
  data = df_gr)

mod_gr_3 <- lm(
  logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance,
  data = df_gr)

mods_gr <- list(
  mod_gr_00, mod_gr_0, mod_gr_10, mod_gr_1,
  mod_gr_20, mod_gr_2, mod_gr_30, mod_gr_3)
mods_gr_dAICc <- bbmle::AICctab(
  mods_gr, weights = TRUE, sort = FALSE)$dAICc
mods_gr_sorted <- order(mods_gr_dAICc)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_best <- mods_gr_sorted[1]
} else {
  mod_gr_index_best <- v_mod_set_gr + 1
}

v_mod_gr_index <- floor((mod_gr_index_best - 1) / 2)
mod_gr_best <- mods_gr[[mod_gr_index_best]]
summary(mod_gr_best)
mods_gr_dAICc


# Growth variance --------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_best)
mod_gr_y <- resid(mod_gr_best)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(
    maxiter = 1000,
    tol = 1e-6,
    warnOnly = TRUE))


# Flowering-probability model --------------------------------------------------
mod_fl_00 <- glm(
  flower ~ 1,
  data = df_fl, family = binomial)

mod_fl_0 <- glm(
  flower ~ disturbance_prev,
  data = df_fl, family = binomial)

mod_fl_10 <- glm(
  flower ~ logvol_t0,
  data = df_fl, family = binomial)

mod_fl_1 <- glm(
  flower ~ logvol_t0 + disturbance_prev,
  data = df_fl, family = binomial)

mod_fl_20 <- glm(
  flower ~ logvol_t0 + logvol_t0_2,
  data = df_fl, family = binomial)

mod_fl_2 <- glm(
  flower ~ logvol_t0 + logvol_t0_2 + disturbance_prev,
  data = df_fl, family = binomial)

mod_fl_30 <- glm(
  flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
  data = df_fl, family = binomial)

mod_fl_3 <- glm(
  flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance_prev,
  data = df_fl, family = binomial)

mods_fl <- list(
  mod_fl_00, mod_fl_0, mod_fl_10, mod_fl_1,
  mod_fl_20, mod_fl_2, mod_fl_30, mod_fl_3)
mods_fl_dAICc <- bbmle::AICctab(
  mods_fl, weights = TRUE, sort = FALSE)$dAICc
mods_fl_sorted <- order(mods_fl_dAICc)

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_best <- mods_fl_sorted[1]
} else {
  mod_fl_index_best <- v_mod_set_fl + 1
}

v_mod_fl_index <- floor((mod_fl_index_best - 1) / 2)
mod_fl_best <- mods_fl[[mod_fl_index_best]]
summary(mod_fl_best)
mods_fl_dAICc


# Flower-number model ----------------------------------------------------------
mod_fl_n_00 <- MASS::glm.nb(
  fl_nr ~ 1,
  data = df_fl_cond)

mod_fl_n_0 <- MASS::glm.nb(
  fl_nr ~ disturbance_prev,
  data = df_fl_cond)

mod_fl_n_10 <- MASS::glm.nb(
  fl_nr ~ logvol_t0,
  data = df_fl_cond)

mod_fl_n_1 <- MASS::glm.nb(
  fl_nr ~ logvol_t0 + disturbance_prev,
  data = df_fl_cond)

mod_fl_n_20 <- MASS::glm.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2,
  data = df_fl_cond)

mod_fl_n_2 <- MASS::glm.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2 + disturbance_prev,
  data = df_fl_cond)

mod_fl_n_30 <- MASS::glm.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
  data = df_fl_cond)

mod_fl_n_3 <- MASS::glm.nb(
  fl_nr ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + disturbance_prev,
  data = df_fl_cond)

mods_fl_n <- list(
  mod_fl_n_00, mod_fl_n_0, mod_fl_n_10, mod_fl_n_1,
  mod_fl_n_20, mod_fl_n_2, mod_fl_n_30, mod_fl_n_3)
mods_fl_n_dAICc <- bbmle::AICctab(
  mods_fl_n, weights = TRUE, sort = FALSE)$dAICc
mods_fl_n_sorted <- order(mods_fl_n_dAICc)

if (length(v_mod_set_fl_n) == 0) {
  mod_fl_n_index_best <- mods_fl_n_sorted[1]
} else {
  mod_fl_n_index_best <- v_mod_set_fl_n + 1
}

v_mod_fl_n_index <- floor((mod_fl_n_index_best - 1) / 2)
mod_fl_n_best <- mods_fl_n[[mod_fl_n_index_best]]
summary(mod_fl_n_best)
mods_fl_n_dAICc


# Mean quadrat-level flowering-stem-to-recruit model ---------------------------
# Recruitment can occur from the persistent seedbank when no flowering stems
# were recorded. Candidate models therefore include intercept-only and
# disturbance-only structures.
mod_re_00 <- lm(
  re_t1_log ~ 1,
  data = df_fs2r)

mod_re_0 <- lm(
  re_t1_log ~ disturbance,
  data = df_fs2r)

mod_re_10 <- lm(
  re_t1_log ~ fs_t0_log,
  data = df_fs2r)

mod_re_1 <- lm(
  re_t1_log ~ fs_t0_log + disturbance,
  data = df_fs2r)

mod_re_20 <- lm(
  re_t1_log ~ fs_t0_log * disturbance,
  data = df_fs2r)

mods_re <- list(
  mod_re_00,
  mod_re_0,
  mod_re_10,
  mod_re_1,
  mod_re_20)

mods_re_dAICc <- bbmle::AICctab(
  mods_re,
  weights = TRUE,
  sort = FALSE)$dAICc

mods_re_sorted <- order(mods_re_dAICc)

if (length(v_mod_set_re) == 0) {
  mod_re_index_best <- mods_re_sorted[1]
} else {
  mod_re_index_best <- v_mod_set_re + 1
}

v_mod_re_index <- mod_re_index_best - 1
mod_re_best <- mods_re[[mod_re_index_best]]

summary(mod_re_best)
mods_re_dAICc


# Mean vital-rate plots --------------------------------------------------------
make_mean_prediction <- function(data, model, disturbance_var) {
  x <- seq(
    min(data$logvol_t0, na.rm = TRUE),
    max(data$logvol_t0, na.rm = TRUE),
    length.out = 200)
  
  pred <- expand_grid(
    logvol_t0 = x,
    disturbance_value = c(0, 1)) %>%
    mutate(
      logvol_t0_2 = logvol_t0^2,
      logvol_t0_3 = logvol_t0^3,
      disturbance_plot = factor(
        disturbance_value,
        levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  pred[[disturbance_var]] <- pred$disturbance_value
  pred$prediction <- predict(
    model,
    newdata = pred,
    type = 'response')
  
  pred
}

make_mean_plot <- function(
    data, model, response, disturbance_var, title, y_lab,
    probability = FALSE) {
  pred <- make_mean_prediction(
    data = data,
    model = model,
    disturbance_var = disturbance_var)
  
  data_plot <- data %>%
    mutate(
      disturbance_plot = factor(
        .data[[disturbance_var]],
        levels = c(0, 1),
        labels = c('No fire', 'Fire')))
  
  point_layer <- if (probability) {
    geom_jitter(
      data = data_plot,
      aes(
        x = logvol_t0,
        y = .data[[response]],
        color = disturbance_plot),
      alpha = 0.15,
      size = 0.7,
      width = 0.05,
      height = 0.05)
  } else {
    geom_point(
      data = data_plot,
      aes(
        x = logvol_t0,
        y = .data[[response]],
        color = disturbance_plot),
      alpha = 0.15,
      size = 0.7)
  }
  
  p <- ggplot() +
    point_layer +
    geom_line(
      data = pred,
      aes(
        x = logvol_t0,
        y = prediction,
        color = disturbance_plot),
      linewidth = 0.8) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    theme_bw() +
    labs(
      title = title,
      subtitle = v_ggp_suffix,
      x = expression('log(volume)'[t0]),
      y = y_lab,
      color = NULL)
  
  if (probability) {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }
  
  p
}

fig_su_mean <- make_mean_plot(
  data = df_su,
  model = mod_su_best,
  response = 'survives',
  disturbance_var = 'disturbance',
  title = 'Mean survival',
  y_lab = 'Survival probability',
  probability = TRUE)

fig_gr_mean <- make_mean_plot(
  data = df_gr,
  model = mod_gr_best,
  response = 'logvol_t1',
  disturbance_var = 'disturbance',
  title = 'Mean growth',
  y_lab = expression('log(volume)'[t1])) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

fig_fl_mean <- make_mean_plot(
  data = df_fl,
  model = mod_fl_best,
  response = 'flower',
  disturbance_var = 'disturbance_prev',
  title = 'Mean flowering probability',
  y_lab = 'Flowering probability',
  probability = TRUE)

fig_fl_n_mean <- make_mean_plot(
  data = df_fl_cond,
  model = mod_fl_n_best,
  response = 'fl_nr',
  disturbance_var = 'disturbance_prev',
  title = 'Mean flowering-stem number',
  y_lab = 'Number of flowering stems')

fig_vital_rates_mean <-
  fig_su_mean + fig_gr_mean + fig_fl_mean + fig_fl_n_mean +
  plot_layout(ncol = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

fig_vital_rates_mean


# Flowering-stem-to-recruit plot ----------------------------------------------
df_fs2r_pred <- expand_grid(
  fs_t0_log = seq(
    0,
    max(df_fs2r$fs_t0_log, na.rm = TRUE),
    length.out = 200),
  disturbance = c(0, 1)) %>%
  mutate(
    re_t1_log = predict(
      mod_re_best,
      newdata = .),
    fs_t0 = expm1(fs_t0_log),
    re_t1 = pmax(0, expm1(re_t1_log)),
    disturbance_plot = factor(
      disturbance,
      levels = c(0, 1),
      labels = c('No fire', 'Fire')))

df_fs2r_plot <- df_fs2r %>%
  mutate(
    disturbance_plot = factor(
      disturbance,
      levels = c(0, 1),
      labels = c('No fire', 'Fire')))

fig_fs2r_mean <- ggplot(
  df_fs2r_plot,
  aes(x = fs_t0, y = re_t1, color = disturbance_plot)) +
  geom_point(alpha = 0.4, size = 1.3) +
  geom_line(
    data = df_fs2r_pred,
    aes(x = fs_t0, y = re_t1, color = disturbance_plot),
    linewidth = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1)) +
  scale_color_manual(
    values = c('No fire' = 'black', 'Fire' = 'red')) +
  theme_bw() +
  labs(
    title = 'Mean flowering-stem-to-recruit relationship',
    subtitle = v_ggp_suffix,
    x = expression('Total flowering stems'[t]),
    y = expression('Recruits'[t+1]),
    color = NULL)

fig_fs2r_mean


# Recruit size distribution ----------------------------------------------------
df_recr <- df %>%
  filter(
    recruits == 1,
    volume_t0 > 0,
    is.finite(logvol_t0))

rc_sz <- df_recr %>%
  summarise(
    mean = mean(logvol_t0, na.rm = TRUE),
    sd = sd(logvol_t0, na.rm = TRUE))


# Parameter extraction ---------------------------------------------------------
size_term_map <- c(
  b0 = '(Intercept)',
  b1 = 'logvol_t0',
  b2 = 'logvol_t0_2',
  b3 = 'logvol_t0_3',
  bd = 'disturbance')

lagged_size_term_map <- c(
  b0 = '(Intercept)',
  b1 = 'logvol_t0',
  b2 = 'logvol_t0_2',
  b3 = 'logvol_t0_3',
  bd = 'disturbance_prev')

extract_terms <- function(model, prefix, term_map) {
  coefs <- stats::coef(model)
  
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

extract_coef <- function(model, term, default = 0) {
  coefs <- stats::coef(model)
  
  if (term %in% names(coefs)) {
    unname(coefs[[term]])
  } else {
    default
  }
}

grow_var_coef <- stats::coef(mod_gr_var)

mesh_limits <- range(
  c(
    df_gr$logvol_t0,
    df_gr$logvol_t1,
    df_recr$logvol_t0),
  na.rm = TRUE,
  finite = TRUE)
mesh_limits <- c(mesh_limits[1] - 0.1, mesh_limits[2] + 0.1)

pars_mean <- as.list(c(
  extract_terms(mod_su_best, 'surv_', size_term_map),
  extract_terms(mod_gr_best, 'grow_', size_term_map),
  extract_terms(mod_fl_best, 'fl_', lagged_size_term_map),
  extract_terms(mod_fl_n_best, 'fln_', lagged_size_term_map),
  re_b0 = extract_coef(
    mod_re_best,
    '(Intercept)'),
  re_b1 = extract_coef(
    mod_re_best,
    'fs_t0_log'),
  re_bd = extract_coef(
    mod_re_best,
    'disturbance'),
  re_bf = extract_coef(
    mod_re_best,
    'fs_t0_log:disturbance'),
  re_sigma = stats::sigma(mod_re_best),
  recr_mean = rc_sz$mean,
  recr_sd = rc_sz$sd,
  grow_var_a = unname(grow_var_coef[['a']]),
  grow_var_b = unname(grow_var_coef[['b']]),
  fs_ref = mean(df_fs2r$fs_t0, na.rm = TRUE),
  L = mesh_limits[1],
  U = mesh_limits[2],
  mat_siz = 200,
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index,
  mod_fl_index = v_mod_fl_index,
  mod_fl_n_index = v_mod_fl_n_index,
  mod_re_index = v_mod_re_index))

pars_mean_table <- tibble(
  parameter = names(pars_mean),
  value = as.numeric(unlist(pars_mean)))

pars_mean_table


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


# Flowering probability --------------------------------------------------------
fl_x <- function(x, pars, disturbance_prev = 0) {
  plogis(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fl_',
    disturbance = disturbance_prev))
}


# Flowering-stem number conditional on flowering ------------------------------
fln_x <- function(x, pars, disturbance_prev = 0) {
  exp(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fln_',
    disturbance = disturbance_prev))
}


# Expected flowering stems per individual ------------------------------------
fs_x <- function(x, pars, disturbance_prev = 0) {
  fl_x(
    x = x,
    pars = pars,
    disturbance_prev = disturbance_prev) *
    fln_x(
      x = x,
      pars = pars,
      disturbance_prev = disturbance_prev)
}


# Flowering-stem-to-recruit conversion ----------------------------------------
predict_recruits_from_stems <- function(
    flowering_stems, pars, disturbance = 0) {
  stems_nonnegative <- pmax(flowering_stems, 0)
  
  slope <- get_par(pars, 're_b1') +
    get_par(pars, 're_bf') * disturbance
  
  eta <- get_par(pars, 're_b0') +
    get_par(pars, 're_bd') * disturbance +
    slope * log1p(stems_nonnegative)
  
  recruits <- expm1(
    eta + 0.5 * get_par(pars, 're_sigma')^2)
  
  pmax(recruits, 0)
}

recruits_per_stem <- function(pars, disturbance = 0) {
  fs_ref <- get_par(pars, 'fs_ref', 0)
  
  if (!is.finite(fs_ref) || fs_ref <= 0) {
    return(0)
  }
  
  predict_recruits_from_stems(
    flowering_stems = fs_ref,
    pars = pars,
    disturbance = disturbance) / fs_ref
}

rx <- function(
    x, pars, disturbance = 0, disturbance_prev = 0) {
  fs_x(
    x = x,
    pars = pars,
    disturbance_prev = disturbance_prev) *
    recruits_per_stem(
      pars = pars,
      disturbance = disturbance)
}


# Recruit size distribution ----------------------------------------------------
recr_y <- function(y, pars) {
  dnorm(
    y,
    mean = get_par(pars, 'recr_mean'),
    sd = get_par(pars, 'recr_sd'))
}


# Fertility matrix -------------------------------------------------------------
fy <- function(
    y, x, pars, h, disturbance = 0, disturbance_prev = 0) {
  Pvec <- recr_y(y, pars)
  recruit_mass <- sum(Pvec * h)
  
  if (!is.finite(recruit_mass) || recruit_mass <= 0) {
    stop('Recruit size distribution has no mass inside the mesh')
  }
  
  Pvec <- Pvec / recruit_mass
  
  Rvec <- rx(
    x = x,
    pars = pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  outer(Pvec, Rvec) * h
}


# Kernel -----------------------------------------------------------------------
kernel <- function(
    pars, disturbance = 0, disturbance_prev = 0) {
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  Fmat <- fy(
    y = y,
    x = y,
    pars = pars,
    h = h,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  Smat <- sx(
    x = y,
    pars = pars,
    disturbance = disturbance)
  
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
  
  Tmat <- matrix(0, n, n)
  
  for (i in 1:(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }
  
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

lambda_ipm <- function(
    pars, disturbance = 0, disturbance_prev = 0) {
  K <- kernel(
    pars = pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)$k_yx
  
  Re(eigen(K, only.values = TRUE)$values[1])
}


# Mean observed size distribution ---------------------------------------------
# The vector is the arithmetic mean of annual size-frequency vectors. Unsized
# recruits are distributed over the fitted recruit-size distribution.
make_initial_n_mean <- function(pars, df_init = df) {
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  breaks <- seq(L, U, length.out = n + 1)
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  recruit_density <- recr_y(y, pars)
  recruit_density <- recruit_density / sum(recruit_density * h)
  
  years_i <- df_init %>%
    filter(!is.na(year)) %>%
    distinct(year) %>%
    arrange(year) %>%
    pull(year)
  
  counts_by_year <- purrr::map(years_i, function(year_i) {
    df_i <- df_init %>%
      filter(year == year_i)
    
    x <- df_i %>%
      filter(is.finite(logvol_t0)) %>%
      pull(logvol_t0)
    
    size_counts <- hist(
      pmin(pmax(x, L), U),
      breaks = breaks,
      plot = FALSE,
      include.lowest = TRUE)$counts
    
    n_unsized_recruits <- df_i %>%
      summarise(
        n = sum(
          !is.finite(logvol_t0) & recruits == 1,
          na.rm = TRUE)) %>%
      pull(n)
    
    size_counts + n_unsized_recruits * recruit_density * h
  })
  
  mean_counts <- Reduce(`+`, counts_by_year) / length(counts_by_year)
  mean_counts / h
}


# Quadrat-level flowering-stem reference --------------------------------------
# The flowering-stem-to-recruit model is fitted to quadrat-year totals.
# Its flowering-stem reference must therefore also be on the quadrat scale.
df_fs2r_previous_fire <- df_fs2r %>%
  left_join(
    df %>%
      distinct(
        site,
        quad,
        year,
        disturbance_prev),
    by = c(
      "site",
      "quad",
      "year")) %>%
  filter(!is.na(disturbance_prev))

df_stem_reference <- df_fs2r_previous_fire %>%
  group_by(disturbance_prev) %>%
  summarise(
    fs_ref = mean(fs_t0),
    n_quadrat_transitions = n(),
    .groups = "drop")

get_stem_reference <- function(disturbance_prev_i) {
  df_stem_reference %>%
    filter(
      disturbance_prev ==
        disturbance_prev_i) %>%
    pull(fs_ref)
}

mean_stems_no_previous_fire <- get_stem_reference(0)
mean_stems_previous_fire <- get_stem_reference(1)

df_stem_reference


# Disturbance history ----------------------------------------------------------
# Each sampled quadrat-year receives equal weight. Current disturbance refers
# to the interval t -> t + 1 and previous disturbance to t - 1 -> t.
df_disturbance_regime <- df %>%
  filter(
    !is.na(quad),
    !is.na(year),
    !is.na(disturbance),
    !is.na(disturbance_prev)) %>%
  distinct(
    quad,
    year,
    disturbance,
    disturbance_prev) %>%
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


# Disturbance transition probabilities ----------------------------------------
# Environmental transitions are conditional probabilities:
# P(disturbance_t | disturbance_t-1).
df_disturbance_transition <- df_disturbance_regime %>%
  group_by(disturbance_prev) %>%
  mutate(
    p_transition = n_quadrat_years /
      sum(n_quadrat_years)) %>%
  ungroup()

df_disturbance_transition

get_transition_probability <- function(
    disturbance_i, disturbance_prev_i) {
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


# Previous-disturbance probabilities ------------------------------------------
# These probabilities are used only to divide the initial observed population
# between the two environmental-history states.
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


# Parameters for previous-disturbance states ----------------------------------
# Flowering and flowering-stem production depend on disturbance in the
# preceding interval, so the flowering-stem reference also follows that state.
pars_prev_0 <- pars_mean
pars_prev_0$fs_ref <- mean_stems_no_previous_fire

pars_prev_1 <- pars_mean
pars_prev_1$fs_ref <- mean_stems_previous_fire


# State-specific demographic kernels ------------------------------------------
# First digit  = current disturbance, t -> t + 1.
# Second digit = previous disturbance, t - 1 -> t.
K_00 <- kernel(
  pars = pars_prev_0,
  disturbance = 0,
  disturbance_prev = 0)$k_yx

K_01 <- kernel(
  pars = pars_prev_1,
  disturbance = 0,
  disturbance_prev = 1)$k_yx

K_10 <- kernel(
  pars = pars_prev_0,
  disturbance = 1,
  disturbance_prev = 0)$k_yx

K_11 <- kernel(
  pars = pars_prev_1,
  disturbance = 1,
  disturbance_prev = 1)$k_yx


# Environmental historical megamatrix ----------------------------------------
# Columns represent disturbance during t - 1 -> t.
# Rows represent disturbance during t -> t + 1.
K_env_hist <- rbind(
  cbind(
    p_0_given_0 * K_00,
    p_0_given_1 * K_01),
  cbind(
    p_1_given_0 * K_10,
    p_1_given_1 * K_11))


# Environmental historical asymptotic lambda ----------------------------------
lambda_env_hist <- Re(
  eigen(
    K_env_hist,
    only.values = TRUE)$values[1])

lambda_env_hist



# Asymptotic and projected lambda ---------------------------------------------
# Mean observed size distribution used for the four fixed environmental states.
n_obs <- make_initial_n_mean(pars_mean)

n <- as.integer(pars_mean$mat_siz)
h <- (pars_mean$U - pars_mean$L) / n

# For the environmental historical IPM, the observed population is divided
# between the two previous-disturbance states. The two blocks are mutually
# exclusive parts of the same population and therefore must not be duplicated.
n_env_initial <- c(
  p_prev_0 * n_obs,
  p_prev_1 * n_obs)


# General projection function -------------------------------------------------
project_kernel <- function(K, n_initial, h, label) {
  n_initial_total <- sum(n_initial) * h
  
  n_proj <- as.numeric(K %*% n_initial)
  n_projected <- sum(n_proj) * h
  
  tibble(
    lambda_type = label,
    asym_lambda = Re(eigen(K, only.values = TRUE)$values[1]),
    proj_lambda = n_projected / n_initial_total,
    n_initial = n_initial_total,
    n_projected = n_projected)
}


# Fixed disturbance-state lambdas ---------------------------------------------
make_state_lambda <- function(
    K, label, disturbance, disturbance_prev) {
  project_kernel(
    K = K,
    n_initial = n_obs,
    h = h,
    label = label) %>%
    mutate(
      disturbance = disturbance,
      disturbance_prev = disturbance_prev)
}

df_lambda_states <- bind_rows(
  make_state_lambda(
    K = K_00,
    label = 'Undisturbed',
    disturbance = 0,
    disturbance_prev = 0),
  make_state_lambda(
    K = K_10,
    label = 'Current-year fire',
    disturbance = 1,
    disturbance_prev = 0),
  make_state_lambda(
    K = K_01,
    label = 'Previous-year fire',
    disturbance = 0,
    disturbance_prev = 1),
  make_state_lambda(
    K = K_11,
    label = 'Fire in both years',
    disturbance = 1,
    disturbance_prev = 1))


# Environmental historical regime --------------------------------------------
df_lambda_env_hist <- project_kernel(
  K = K_env_hist,
  n_initial = n_env_initial,
  h = h,
  label = 'Environmental historical regime') %>%
  mutate(
    disturbance = NA_real_,
    disturbance_prev = NA_real_)


# Combined lambda table --------------------------------------------------------
df_lambda_env <- bind_rows(
  df_lambda_states,
  df_lambda_env_hist)

df_lambda_env


# Observed annual abundance change --------------------------------------------
# Counts include sized plants and recruits without a usable size measurement.
# Only quadrats with individual records in both consecutive years are retained.
df_counts_quad_year <- df %>%
  filter(
    !is.na(site),
    !is.na(quad),
    !is.na(year)) %>%
  group_by(site, quad, year) %>%
  summarise(
    n_sized = sum(
      is.finite(logvol_t0),
      na.rm = TRUE),
    n_unsized_recruits = sum(
      !is.finite(logvol_t0) & recruits == 1,
      na.rm = TRUE),
    n_t0 = n_sized + n_unsized_recruits,
    .groups = "drop") %>%
  arrange(site, quad, year) %>%
  group_by(site, quad) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n_t0),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1)

df_obs_growth_year <- df_counts_quad_year %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0),
    n_t1 = sum(n_t1),
    n_quadrats = n(),
    n_sites = n_distinct(site),
    obs_lambda = n_t1 / n_t0,
    .groups = "drop") %>%
  filter(n_t0 > 0)

df_obs_growth_year


df_obs_growth_summary <- df_obs_growth_year %>%
  summarise(
    n_transitions = n(),
    lambda_arithmetic = mean(obs_lambda, na.rm = TRUE),
    lambda_geometric = exp(mean(log(obs_lambda), na.rm = TRUE)),
    lambda_median = median(obs_lambda, na.rm = TRUE),
    lambda_min = min(obs_lambda, na.rm = TRUE),
    lambda_max = max(obs_lambda, na.rm = TRUE))

df_obs_growth_summary


# Observed and modeled population growth through time -------------------------
df_model_lines <- df_lambda_env %>%
  dplyr::select(
    lambda_type,
    asym_lambda,
    proj_lambda) %>%
  pivot_longer(
    cols = c(
      asym_lambda,
      proj_lambda),
    names_to = 'metric',
    values_to = 'lambda') %>%
  mutate(
    metric = recode(
      metric,
      asym_lambda = 'Asymptotic',
      proj_lambda = 'One-step projected'))

x_point <- max(
  df_obs_growth_year$year,
  na.rm = TRUE) + 0.4

fig_lambda_env <- ggplot(
  df_obs_growth_year,
  aes(
    x = year,
    y = obs_lambda)) +
  geom_hline(
    yintercept = 1,
    linetype = 'dotted') +
  geom_line(
    color = 'black') +
  geom_point(
    size = 2,
    color = 'black') +
  geom_hline(
    data = df_model_lines,
    aes(
      yintercept = lambda,
      color = lambda_type),
    linewidth = 0.8) +
  geom_point(
    data = df_model_lines,
    aes(
      x = x_point,
      y = lambda,
      color = lambda_type,
      shape = metric),
    size = 3) +
  expand_limits(
    x = x_point + 1) +
  theme_bw() +
  labs(
    title = paste(
      'Observed annual growth and environmental-history',
      'IPM lambda'),
    subtitle = v_ggp_suffix,
    x = 'Year t0',
    y = expression(lambda),
    color = NULL,
    shape = NULL)

fig_lambda_env


# Kernel diagnostics -----------------------------------------------------------
df_kernel_diagnostics <- tibble(
  kernel = c(
    'Undisturbed',
    'Current-year fire',
    'Previous-year fire',
    'Fire in both years',
    'Environmental historical regime'),
  n_rows = c(
    nrow(K_00),
    nrow(K_10),
    nrow(K_01),
    nrow(K_11),
    nrow(K_env_hist)),
  n_cols = c(
    ncol(K_00),
    ncol(K_10),
    ncol(K_01),
    ncol(K_11),
    ncol(K_env_hist)),
  min_entry = c(
    min(K_00),
    min(K_10),
    min(K_01),
    min(K_11),
    min(K_env_hist)),
  max_entry = c(
    max(K_00),
    max(K_10),
    max(K_01),
    max(K_11),
    max(K_env_hist)),
  any_nonfinite = c(
    any(!is.finite(K_00)),
    any(!is.finite(K_10)),
    any(!is.finite(K_01)),
    any(!is.finite(K_11)),
    any(!is.finite(K_env_hist))))

df_kernel_diagnostics


# Environmental historical structure -----------------------------------------
df_env_hist_structure <- tibble(
  block = c(
    '0 <- 0',
    '0 <- 1',
    '1 <- 0',
    '1 <- 1'),
  correct = c(
    isTRUE(all.equal(
      K_env_hist[1:n, 1:n],
      p_0_given_0 * K_00)),
    isTRUE(all.equal(
      K_env_hist[1:n, (n + 1):(2 * n)],
      p_0_given_1 * K_01)),
    isTRUE(all.equal(
      K_env_hist[(n + 1):(2 * n), 1:n],
      p_1_given_0 * K_10)),
    isTRUE(all.equal(
      K_env_hist[(n + 1):(2 * n), (n + 1):(2 * n)],
      p_1_given_1 * K_11))))

df_env_hist_structure
