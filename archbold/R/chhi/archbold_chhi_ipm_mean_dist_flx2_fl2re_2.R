# Mean IPM with disturbance - Archbold - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.24

# This script fits mean vital-rate models with disturbance as a fixed effect.
# It contains no year random effects and no site random effects.
# The fertility pathway is:
# flowering probability x flowerhead number conditional on flowering x
# flowerhead-to-recruit conversion.

# rm(list = ls())


# Setting the stage ------------------------------------------------------------
set.seed(100)
options(stringsAsFactors = FALSE)


# Packages ---------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(
  MASS,
  tidyverse,
  bbmle,
  patchwork,
  janitor,
  scales)


# Specification ----------------------------------------------------------------
v_head <- 'archbold'
v_species <- 'Chrysopsis highlandsensis'
v_years_re <- c()

v_sp_abb <- tolower(
  gsub(
    ' ', '',
    paste(
      substr(unlist(strsplit(v_species, ' ')), 1, 2),
      collapse = '')))

v_script_prefix <- v_head
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Empty vectors select the best-supported model by AICc.
v_mod_set_su <- c()
v_mod_set_gr <- c()
v_mod_set_fl <- c()
v_mod_set_fl_n <- c()


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
df_og <- read_csv(
  file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>%
  janitor::clean_names()


df_fire <- read_csv(
  file.path(dir_data, 'chrysopsis_highlandsensis_fire.csv')) %>%
  janitor::clean_names() %>%
  rename(
    year = burn_yr,
    fire = treatment) %>%
  select(-any_of(c('year0', 'notes')))


# Recruitment classification --------------------------------------------------
# Recruits include:
# 1. plants classified as seedlings;
# 2. plants first observed after site monitoring began and recorded as a
#    new plant in at least one quarterly census.
#
# New plants from a site's initial census are excluded because their presence
# before monitoring began is unknown.
#
# Diagnostic checks showed that seedling recruits were never counted in more
# than one year and always occurred in the individual's first observed year.
# Quarterly census records identified 1,109 first-observed non-seedlings as
# new plants. Of these, 828 occurred in a site's initial census and were
# excluded. The remaining 280 appeared after monitoring had begun and were
# included as recruits; 264 were vegetative and 16 were bolting at the annual
# census.

df_recruit_flag <- df_og %>%
  rename(
    plant_id = identifier,
    year = year0) %>%
  filter(
    !is.na(site),
    !is.na(plant_id),
    !is.na(year)) %>%
  mutate(
    new_plant_status = if_any(
      c(s_03, s_06, s_09, s_12), ~ .x == 3)) %>%
  group_by(site) %>%
  mutate(site_first_year = min(year, na.rm = TRUE)) %>%
  group_by(site, plant_id) %>%
  mutate(first_plant_year = min(year, na.rm = TRUE)) %>%
  group_by(site, plant_id, year) %>%
  summarise(
    seedling = any(astg == 1, na.rm = TRUE),
    new_plant_status = any(new_plant_status, na.rm = TRUE),
    site_first_year = first(site_first_year),
    first_plant_year = first(first_plant_year),
    .groups = 'drop') %>%
  mutate(
    recruit = as.integer(
      seedling |
        (new_plant_status &
          year == first_plant_year &
          year > site_first_year)),
    recruit_source = case_when(
      seedling ~ 'Seedling stage',
      recruit == 1 ~ 'New plant after site establishment',
      TRUE ~ 'Not recruit'))


# Individual-level data --------------------------------------------------------
df <- df_og %>%
  rename(
    plant_id = identifier,
    year = year0,
    survival = survival_1) %>%
  arrange(site, plant_id, year, survival) %>%
  group_by(site, plant_id, year) %>%
  summarise(
    survives = if (all(is.na(survival))) {
      NA_real_
    } else {
      min(survival, na.rm = TRUE)
    },
    size_t0 = if (all(is.na(dia))) {
      NA_real_
    } else {
      max(dia, na.rm = TRUE)
    },
    size_t1 = if (all(is.na(dia_1))) {
      NA_real_
    } else {
      max(dia_1, na.rm = TRUE)
    },
    fl_nr = if (all(is.na(hd))) {
      NA_real_
    } else {
      max(hd, na.rm = TRUE)
    },
    .groups = 'drop') %>%
  left_join(
    df_recruit_flag %>%
      select(site, plant_id, year, recruit, recruit_source),
    by = c('site', 'plant_id', 'year')) %>%
  mutate(
    plant_id = factor(plant_id),
    recruit = replace_na(recruit, 0L),
    recruit_source = replace_na(recruit_source, 'Not recruit'),
    flower = if_else(!is.na(fl_nr) & fl_nr > 0, 1, 0),
    logsize_t0 = log(size_t0),
    logsize_t1 = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  full_join(df_fire, by = c('site', 'year')) %>%
  mutate(
    fire = case_when(
      is.na(fire) ~ 'No fire',
      fire == 'burn' ~ 'Fire',
      TRUE ~ NA_character_),
    site = factor(site),
    fire = factor(fire, levels = c('No fire', 'Fire')),
    disturbance = as.integer(fire == 'Fire'))


# Vital-rate data --------------------------------------------------------------
df_su <- df %>%
  filter(
    size_t0 > 0,
    !is.na(survives),
    is.finite(logsize_t0),
    !is.na(disturbance))


df_gr <- df %>%
  filter(
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    !is.na(disturbance))


# Flowering probability --------------------------------------------------------
df_fl <- df %>%
  filter(
    !is.na(fl_nr),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(flower = if_else(fl_nr > 0, 1, 0))


# Flowerhead number conditional on flowering ----------------------------------
df_fl_n <- df_fl %>%
  filter(
    flower == 1,
    fl_nr > 0,
    is.finite(fl_nr),
    fl_nr == round(fl_nr))


# Flowerhead-to-recruit data ---------------------------------------------------
# Total flowerheads in year t predict recruits in year t + 1. The mean model is
# forced through zero and contains no year or site random effects.
df_ind <- df %>%
  filter(
    !is.na(site),
    !is.na(plant_id),
    !is.na(year)) %>%
  transmute(
    site,
    plant_id,
    year,
    fh_nr = replace_na(fl_nr, 0),
    recruit = replace_na(recruit, 0L))


df_site_year <- df_ind %>%
  group_by(site, year) %>%
  summarise(
    flowerheads_t0 = sum(fh_nr[recruit == 0], na.rm = TRUE),
    recruits_t0 = sum(recruit == 1, na.rm = TRUE),
    .groups = 'drop')


df_fl2re <- df_site_year %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(
    recruits_t1 = lead(recruits_t0),
    year_t1 = lead(year),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1) %>%
  mutate(
    site = factor(site),
    logflowerheads_t0 = log1p(flowerheads_t0),
    logrecruits_t1 = log1p(recruits_t1))


# Remove selected years --------------------------------------------------------
v_years_og <- sort(unique(df$year))

df <- df %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_su <- df_su %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_gr <- df_gr %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_fl <- df_fl %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_fl_n <- df_fl_n %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_fl2re <- df_fl2re %>%
  filter(
    !is.na(year),
    !is.na(year_t1),
    !(year %in% v_years_re),
    !(year_t1 %in% v_years_re))

v_years <- sort(unique(df$year))


# Survival model ---------------------------------------------------------------
mod_su_00 <- glm(
  survives ~ 1,
  data = df_su, family = binomial)

mod_su_0 <- glm(
  survives ~ disturbance,
  data = df_su, family = binomial)

mod_su_10 <- glm(
  survives ~ logsize_t0,
  data = df_su, family = binomial)

mod_su_1 <- glm(
  survives ~ logsize_t0 + disturbance,
  data = df_su, family = binomial)

mod_su_20 <- glm(
  survives ~ logsize_t0 + logsize_t0_2,
  data = df_su, family = binomial)

mod_su_2 <- glm(
  survives ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_su, family = binomial)

mod_su_30 <- glm(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_su, family = binomial)

mod_su_3 <- glm(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
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
mod_su_coef <- coef(mod_su_best)
summary(mod_su_best)


# Growth model -----------------------------------------------------------------
mod_gr_00 <- lm(
  logsize_t1 ~ 1,
  data = df_gr)

mod_gr_0 <- lm(
  logsize_t1 ~ disturbance,
  data = df_gr)

mod_gr_10 <- lm(
  logsize_t1 ~ logsize_t0,
  data = df_gr)

mod_gr_1 <- lm(
  logsize_t1 ~ logsize_t0 + disturbance,
  data = df_gr)

mod_gr_20 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2,
  data = df_gr)

mod_gr_2 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_gr)

mod_gr_30 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_gr)

mod_gr_3 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 +
    logsize_t0_3 + disturbance,
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
mod_gr_coef <- coef(mod_gr_best)
summary(mod_gr_best)


# Growth variance --------------------------------------------------------------
# Variance is modeled against the predicted growth mean. grow_sd() below uses
# the same mean-scale predictor.
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
  flower ~ disturbance,
  data = df_fl, family = binomial)

mod_fl_10 <- glm(
  flower ~ logsize_t0,
  data = df_fl, family = binomial)

mod_fl_1 <- glm(
  flower ~ logsize_t0 + disturbance,
  data = df_fl, family = binomial)

mod_fl_20 <- glm(
  flower ~ logsize_t0 + logsize_t0_2,
  data = df_fl, family = binomial)

mod_fl_2 <- glm(
  flower ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_fl, family = binomial)

mod_fl_30 <- glm(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fl, family = binomial)

mod_fl_3 <- glm(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
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
mod_fl_coef <- coef(mod_fl_best)
summary(mod_fl_best)


# Flowerhead-number model ------------------------------------------------------
mod_fl_n_00 <- MASS::glm.nb(
  fl_nr ~ 1,
  data = df_fl_n)

mod_fl_n_0 <- MASS::glm.nb(
  fl_nr ~ disturbance,
  data = df_fl_n)

mod_fl_n_10 <- MASS::glm.nb(
  fl_nr ~ logsize_t0,
  data = df_fl_n)

mod_fl_n_1 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + disturbance,
  data = df_fl_n)

mod_fl_n_20 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2,
  data = df_fl_n)

mod_fl_n_2 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_fl_n)

mod_fl_n_30 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fl_n)

mod_fl_n_3 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
  data = df_fl_n)

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
mod_fl_n_coef <- coef(mod_fl_n_best)
summary(mod_fl_n_best)


# Mean flowerhead-to-recruit model --------------------------------------------
# Zero flowerheads predict zero recruits. No site or year random effects are
# fitted.
mod_fl2re_mean <- lm(
  logrecruits_t1 ~ 0 + logflowerheads_t0,
  data = df_fl2re)

summary(mod_fl2re_mean)


# Mean vital-rate plots --------------------------------------------------------
make_mean_prediction <- function(data, model) {
  x <- seq(
    min(data$logsize_t0, na.rm = TRUE),
    max(data$logsize_t0, na.rm = TRUE),
    length.out = 200)

  pred <- expand_grid(
    logsize_t0 = x,
    disturbance = c(0, 1)) %>%
    mutate(
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3,
      disturbance_plot = factor(
        disturbance,
        levels = c(0, 1),
        labels = c('No fire', 'Fire')))

  pred$prediction <- predict(
    model,
    newdata = pred,
    type = 'response')

  pred
}


make_mean_plot <- function(
    data, model, response, title, y_lab,
    probability = FALSE) {
  pred <- make_mean_prediction(data, model)

  data_plot <- data %>%
    mutate(
      disturbance_plot = factor(
        disturbance,
        levels = c(0, 1),
        labels = c('No fire', 'Fire')))

  point_layer <- if (probability) {
    geom_jitter(
      data = data_plot,
      aes(
        x = logsize_t0,
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
        x = logsize_t0,
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
        x = logsize_t0,
        y = prediction,
        color = disturbance_plot),
      linewidth = 0.8) +
    scale_color_manual(
      values = c('No fire' = 'black', 'Fire' = 'red')) +
    theme_bw() +
    labs(
      title = title,
      subtitle = v_ggp_suffix,
      x = expression('log(size)'[t0]),
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
  title = 'Mean survival',
  y_lab = 'Survival probability',
  probability = TRUE)

fig_gr_mean <- make_mean_plot(
  data = df_gr,
  model = mod_gr_best,
  response = 'logsize_t1',
  title = 'Mean growth',
  y_lab = expression('log(size)'[t1]))

fig_fl_mean <- make_mean_plot(
  data = df_fl,
  model = mod_fl_best,
  response = 'flower',
  title = 'Mean flowering probability',
  y_lab = 'Flowering probability',
  probability = TRUE)

fig_fl_n_mean <- make_mean_plot(
  data = df_fl_n,
  model = mod_fl_n_best,
  response = 'fl_nr',
  title = 'Mean flowerhead number given flowering',
  y_lab = 'Number of flowerheads')

fig_vital_rates_mean <-
  fig_su_mean + fig_gr_mean + fig_fl_mean + fig_fl_n_mean +
  plot_layout(ncol = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

fig_vital_rates_mean


# Flowerhead-to-recruit plot ---------------------------------------------------
df_fl2re_pred <- tibble(
  logflowerheads_t0 = seq(
    0,
    max(df_fl2re$logflowerheads_t0, na.rm = TRUE),
    length.out = 200)) %>%
  mutate(
    logrecruits_t1 = predict(
      mod_fl2re_mean,
      newdata = .),
    flowerheads_t0 = expm1(logflowerheads_t0),
    recruits_t1 = pmax(0, expm1(logrecruits_t1)))

fig_fl2re_mean <- ggplot(
  df_fl2re,
  aes(x = flowerheads_t0, y = recruits_t1)) +
  geom_point(alpha = 0.4, size = 1.3) +
  geom_line(
    data = df_fl2re_pred,
    aes(x = flowerheads_t0, y = recruits_t1),
    linewidth = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1)) +
  theme_bw() +
  labs(
    title = 'Mean flowerhead-to-recruit relationship',
    subtitle = v_ggp_suffix,
    x = expression('Total flowerheads'[t]),
    y = expression('Recruits'[t+1]))

fig_fl2re_mean


# Recruit size distribution ----------------------------------------------------
df_recr <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0))

rc_sz <- df_recr %>%
  summarise(
    mean = mean(logsize_t0, na.rm = TRUE),
    sd = sd(logsize_t0, na.rm = TRUE))


# Parameter extraction ---------------------------------------------------------
size_term_map <- c(
  b0 = '(Intercept)',
  b1 = 'logsize_t0',
  b2 = 'logsize_t0_2',
  b3 = 'logsize_t0_3',
  bd = 'disturbance')


extract_terms <- function(model, prefix, term_map = size_term_map) {
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


grow_var_coef <- stats::coef(mod_gr_var)

mesh_limits <- range(
  c(
    df_gr$logsize_t0,
    df_gr$logsize_t1,
    df_recr$logsize_t0),
  na.rm = TRUE,
  finite = TRUE)

pars_mean <- as.list(c(
  extract_terms(mod_su_best, 'surv_'),
  extract_terms(mod_gr_best, 'grow_'),
  extract_terms(mod_fl_best, 'fl_'),
  extract_terms(mod_fl_n_best, 'fln_'),
  fl2re_b0 = 0,
  fl2re_b1 = unname(
    stats::coef(mod_fl2re_mean)[['logflowerheads_t0']]),
  fl2re_sigma = stats::sigma(mod_fl2re_mean),
  recr_mean = rc_sz$mean,
  recr_sd = rc_sz$sd,
  grow_var_a = unname(grow_var_coef[['a']]),
  grow_var_b = unname(grow_var_coef[['b']]),
  flowerheads_ref = mean(df_fl2re$flowerheads_t0, na.rm = TRUE),
  L = mesh_limits[1],
  U = mesh_limits[2],
  mat_siz = 200))

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


# Flowerhead number conditional on flowering ----------------------------------
fln_x <- function(x, pars, disturbance = 0) {
  exp(vital_lp(
    x = x,
    pars = pars,
    prefix = 'fln_',
    disturbance = disturbance))
}


# Expected flowerheads per individual -----------------------------------------
flowerheads_x <- function(x, pars, disturbance = 0) {
  fl_x(x, pars, disturbance) *
    fln_x(x, pars, disturbance)
}


# Flowerhead-to-recruit conversion --------------------------------------------
predict_recruits_from_flowerheads <- function(flowerheads, pars) {
  flowerheads_nonnegative <- pmax(flowerheads, 0)

  eta <- get_par(pars, 'fl2re_b0') +
    get_par(pars, 'fl2re_b1') * log1p(flowerheads_nonnegative)

  recruits <- expm1(
    eta + 0.5 * get_par(pars, 'fl2re_sigma')^2)

  ifelse(
    flowerheads_nonnegative <= 0,
    0,
    pmax(recruits, 0))
}


recruits_per_flowerhead <- function(pars) {
  flowerheads_ref <- get_par(pars, 'flowerheads_ref', 0)

  if (!is.finite(flowerheads_ref) || flowerheads_ref <= 0) {
    return(0)
  }

  predict_recruits_from_flowerheads(flowerheads_ref, pars) /
    flowerheads_ref
}


rx <- function(x, pars, disturbance = 0) {
  flowerheads_x(
    x = x,
    pars = pars,
    disturbance = disturbance) *
    recruits_per_flowerhead(pars)
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

  Fmat <- fy(
    y = y,
    x = y,
    pars = pars,
    h = h,
    disturbance = disturbance)

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


lambda_ipm <- function(pars, disturbance = 0) {
  K <- kernel(
    pars = pars,
    disturbance = disturbance)$k_yx

  Re(eigen(K, only.values = TRUE)$values[1])
}


# Mean observed size distribution ---------------------------------------------
# The vector is the arithmetic mean of annual size-frequency vectors.
make_initial_n_mean <- function(pars, df_init = df) {
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  breaks <- seq(L, U, length.out = n + 1)

  years_i <- df_init %>%
    filter(size_t0 > 0, is.finite(logsize_t0)) %>%
    distinct(year) %>%
    arrange(year) %>%
    pull(year)

  counts_by_year <- purrr::map(years_i, function(year_i) {
    x <- df_init %>%
      filter(
        year == year_i,
        size_t0 > 0,
        is.finite(logsize_t0)) %>%
      pull(logsize_t0)

    hist(
      pmin(pmax(x, L), U),
      breaks = breaks,
      plot = FALSE,
      include.lowest = TRUE)$counts
  })

  mean_counts <- Reduce(`+`, counts_by_year) / length(counts_by_year)
  mean_counts / h
}


# Expected total flowerheads from the mean size distribution ------------------
expected_mean_flowerheads <- function(
    pars, disturbance = 0, df_init = df) {
  n_obs <- make_initial_n_mean(pars, df_init)
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n

  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  sum(
    flowerheads_x(y, pars, disturbance) * n_obs * h,
    na.rm = TRUE)
}


mean_flowerheads_no_fire <- expected_mean_flowerheads(
  pars = pars_mean,
  disturbance = 0)

mean_flowerheads_fire <- expected_mean_flowerheads(
  pars = pars_mean,
  disturbance = 1)


# Mean disturbance exposure ---------------------------------------------------
# Each represented site-year receives equal weight; plant-rich sites do not
# receive greater weight merely because they contain more plant records.
p_disturbance_mean <- df %>%
  filter(
    !is.na(site),
    !is.na(year),
    !is.na(disturbance)) %>%
  distinct(site, year, disturbance) %>%
  summarise(
    p_disturbance = mean(disturbance, na.rm = TRUE)) %>%
  pull(p_disturbance)

mean_flowerheads_regime <-
  (1 - p_disturbance_mean) * mean_flowerheads_no_fire +
  p_disturbance_mean * mean_flowerheads_fire


# State-specific and disturbance-regime kernels -------------------------------
pars_no_fire <- pars_mean
pars_no_fire$flowerheads_ref <- mean_flowerheads_no_fire

pars_fire <- pars_mean
pars_fire$flowerheads_ref <- mean_flowerheads_fire

# Both component kernels in the mixture use the same expected flowerhead
# production under the observed mean disturbance regime.
pars_regime <- pars_mean
pars_regime$flowerheads_ref <- mean_flowerheads_regime

K_no_fire <- kernel(
  pars = pars_no_fire,
  disturbance = 0)$k_yx

K_fire <- kernel(
  pars = pars_fire,
  disturbance = 1)$k_yx

K_regime_no_fire <- kernel(
  pars = pars_regime,
  disturbance = 0)$k_yx

K_regime_fire <- kernel(
  pars = pars_regime,
  disturbance = 1)$k_yx

K_regime <-
  (1 - p_disturbance_mean) * K_regime_no_fire +
  p_disturbance_mean * K_regime_fire


# Asymptotic and projected lambda ---------------------------------------------
project_kernel <- function(K, pars, label, df_init = df) {
  n_obs <- make_initial_n_mean(pars, df_init)
  h <- (pars$U - pars$L) / pars$mat_siz
  n_initial <- sum(n_obs) * h

  n_proj <- K %*% n_obs
  n_projected <- sum(n_proj) * h

  tibble(
    lambda_type = label,
    asym_lambda = Re(eigen(K, only.values = TRUE)$values[1]),
    proj_lambda = n_projected / n_initial,
    n_initial = n_initial,
    n_projected = n_projected)
}


df_lambda_mean <- bind_rows(
  project_kernel(
    K = K_no_fire,
    pars = pars_no_fire,
    label = 'No fire'),
  project_kernel(
    K = K_fire,
    pars = pars_fire,
    label = 'Fire'),
  project_kernel(
    K = K_regime,
    pars = pars_regime,
    label = 'Expected disturbance regime')) %>%
  mutate(
    p_disturbance = c(0, 1, p_disturbance_mean),
    flowerheads_ref = c(
      mean_flowerheads_no_fire,
      mean_flowerheads_fire,
      mean_flowerheads_regime),
    recruits_ref = purrr::map2_dbl(
      flowerheads_ref,
      list(pars_no_fire, pars_fire, pars_regime),
      predict_recruits_from_flowerheads),
    recruits_per_flowerhead = if_else(
      flowerheads_ref > 0,
      recruits_ref / flowerheads_ref,
      0))


df_lambda_mean


# Observed annual population growth -------------------------------------------
# Counts are paired within sites in consecutive years before being pooled to
# the annual population level.
df_counts_site_year <- df %>%
  filter(
    size_t0 > 0,
    is.finite(logsize_t0),
    !is.na(site),
    !is.na(year)) %>%
  group_by(site, year) %>%
  summarise(
    n_t0 = n_distinct(plant_id),
    .groups = 'drop')


df_obs_growth_site <- df_counts_site_year %>%
  arrange(site, year) %>%
  group_by(site) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n_t0),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1)


df_obs_growth_year <- df_obs_growth_site %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    n_sites = n(),
    obs_lambda = n_t1 / n_t0,
    .groups = 'drop') %>%
  filter(n_t0 > 0, n_t1 > 0)


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
df_model_lines <- df_lambda_mean %>%
  select(lambda_type, asym_lambda, proj_lambda) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = 'metric',
    values_to = 'lambda') %>%
  mutate(
    metric = recode(
      metric,
      asym_lambda = 'Asymptotic',
      proj_lambda = 'Projected'),
    line_type = paste(lambda_type, metric, sep = ' - '))


fig_lambda_mean <- ggplot(
  df_obs_growth_year,
  aes(x = year, y = obs_lambda)) +
  geom_hline(
    yintercept = 1,
    linetype = 'dotted') +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(
    data = df_model_lines,
    aes(
      yintercept = lambda,
      linetype = line_type),
    linewidth = 0.6) +
  theme_bw() +
  labs(
    title = 'Observed annual growth and mean-IPM lambda',
    subtitle = v_ggp_suffix,
    x = 'Year t0',
    y = expression(lambda),
    linetype = NULL)


fig_lambda_mean


# Kernel diagnostics -----------------------------------------------------------
df_kernel_diagnostics <- tibble(
  kernel = c('No fire', 'Fire', 'Expected regime'),
  min_entry = c(
    min(K_no_fire),
    min(K_fire),
    min(K_regime)),
  max_entry = c(
    max(K_no_fire),
    max(K_fire),
    max(K_regime)),
  any_nonfinite = c(
    any(!is.finite(K_no_fire)),
    any(!is.finite(K_fire)),
    any(!is.finite(K_regime))))


df_kernel_diagnostics
