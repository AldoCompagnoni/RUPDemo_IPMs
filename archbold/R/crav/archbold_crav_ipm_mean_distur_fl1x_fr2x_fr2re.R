# Mean IPM with disturbance - Archbold - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diana Spurite, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.07.23

# This script fits mean vital-rate models with disturbance as a fixed effect.
# It contains no year random effects and no site random effects.
# Year, site, and quadrat identifiers are retained only for data preparation,
# recruitment timing, disturbance exposure, and observed-growth validation.

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
  binom,
  janitor,
  scales)


# Specification ----------------------------------------------------------------
v_head <- 'archbold'
v_species <- 'Crotalaria avonensis'
v_years_re <- c()

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
v_mod_set_fr <- c()
v_mod_set_fr_n <- c()


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


# Functions --------------------------------------------------------------------
source('helper_functions/plot_binned_prop.R')


# Data -------------------------------------------------------------------------
df_og <- read_delim(
  file.path(dir_data, 'crotalaria_avonensis_data_v260515.csv'),
  delim = ';', escape_double = FALSE, trim_ws = TRUE)

df_og_quad <- read_delim(
  file.path(dir_data, 'crotalaria_avonensis_data_quad_v260515.csv'),
  delim = ';', escape_double = FALSE, trim_ws = TRUE)


# Disturbance data -------------------------------------------------------------
# Disturbance is assigned to year t0 and affects the t0-to-t1 transition.
df_disturbance <- df_og_quad %>%
  select(
    site, mp, quad,
    burn221, burn1222, burn0123, burn05,
    burn2014_2015, burn2016, burn0217,
    burnoct2017, burnjun2018) %>%
  pivot_longer(
    cols = starts_with('burn'),
    names_to = 'burn_event',
    values_to = 'value') %>%
  mutate(
    disturbance = case_when(
      is.na(value) ~ 0,
      value == FALSE ~ 0,
      value == TRUE ~ 1,
      suppressWarnings(as.numeric(value)) > 0 ~ 1,
      TRUE ~ 0)) %>%
  filter(disturbance == 1) %>%
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
    site = factor(site),
    mp = factor(mp),
    quad = factor(quad)) %>%
  select(site, mp, quad, year, disturbance) %>%
  distinct() %>%
  arrange(site, mp, quad, year)


# Individual-level data --------------------------------------------------------
df <- df_og %>%
  janitor::clean_names() %>%
  rename(
    plant_id = id,
    size_t0 = maxbr,
    flower = maxfl,
    fruit = maxfr) %>%
  mutate(
    plant_id = factor(plant_id),
    site = factor(site),
    mp = factor(mp),
    quad = factor(quad),
    quad_id = factor(paste(site, mp, quad, sep = '_'))) %>%
  arrange(site, quad_id, plant_id, year) %>%
  group_by(plant_id) %>%
  mutate(
    survives = lead(alive),
    size_t1 = case_when(
      survives == 1 ~ lead(size_t0),
      TRUE ~ NA_real_),
    latest_alive_date = if_else(
      any(s > 0 & s != 2 & s < 6, na.rm = TRUE),
      max(year[s > 0 & s != 2 & s < 6], na.rm = TRUE),
      NA_real_),
    earliest_recorded_date = if_else(
      any(s > 0 & s != 2 & s != 6, na.rm = TRUE),
      min(year[s > 0 & s != 2 & s != 6], na.rm = TRUE),
      NA_real_),
    recruit = case_when(
      stage == 1 ~ 1,
      year == earliest_recorded_date &
        earliest_recorded_date > 2003 ~ 1,
      TRUE ~ NA_real_)) %>%
  ungroup() %>%
  group_by(plant_id) %>%
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0) ~ 1,
      size_t0 > 0 & !is.na(survives) ~ 0,
      survives == 0 & year <= max(df_og$year) - 4 ~ 0,
      TRUE ~ NA_real_),
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 &
        lag(dormancy, 2) == 1 & lag(dormancy, 3) == 1 &
        lag(dormancy, 4) == 1 & lag(dormancy, 5) == 1 ~ 6,
      dormancy == 1 & lag(dormancy, 1) == 1 &
        lag(dormancy, 2) == 1 & lag(dormancy, 3) == 1 &
        lag(dormancy, 4) == 1 ~ 5,
      dormancy == 1 & lag(dormancy, 1) == 1 &
        lag(dormancy, 2) == 1 & lag(dormancy, 3) == 1 ~ 4,
      dormancy == 1 & lag(dormancy, 1) == 1 &
        lag(dormancy, 2) == 1 ~ 3,
      dormancy == 1 & lag(dormancy, 1) == 1 &
        lag(dormancy, 2) == 0 ~ 2,
      dormancy == 1 ~ 1,
      TRUE ~ dormancy)) %>%
  ungroup() %>%
  mutate(
    survives = if_else(
      survives == 0 & year > max(df_og$year) - 4,
      NA_real_, survives),
    logsize_t0 = log(size_t0),
    logsize_t1 = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  left_join(
    df_disturbance,
    by = c('site', 'mp', 'quad', 'year')) %>%
  mutate(disturbance = replace_na(disturbance, 0)) %>%
  select(
    site, mp, quad, quad_id, plant_id, year,
    s, survives, size_t0, size_t1,
    logsize_t0, logsize_t1,
    logsize_t0_2, logsize_t0_3,
    flower, fruit, recruit, dormancy, disturbance)


# Recruit validity ------------------------------------------------------------
# A recruit is retained only if its quadrat was represented in the previous
# year. This avoids treating first observations in newly sampled quadrats as
# known recruits.
df_plot_year_sampled <- df %>%
  mutate(year = as.integer(year)) %>%
  group_by(site, mp, quad, quad_id, year) %>%
  summarise(
    n_records_plot_year = n(),
    n_plants_plot_year = n_distinct(plant_id),
    .groups = 'drop')

df_recruit_valid_ids <- df %>%
  filter(recruit == 1) %>%
  mutate(
    year = as.integer(year),
    year_prev = year - 1) %>%
  left_join(
    df_plot_year_sampled %>%
      transmute(
        quad_id,
        year_prev = year,
        sampled_prev_year = TRUE),
    by = c('quad_id', 'year_prev')) %>%
  mutate(sampled_prev_year = replace_na(sampled_prev_year, FALSE)) %>%
  filter(sampled_prev_year) %>%
  select(site, mp, quad, quad_id, plant_id, year) %>%
  mutate(recruit_plot_valid = 1)

df <- df %>%
  mutate(year = as.integer(year)) %>%
  left_join(
    df_recruit_valid_ids,
    by = c('site', 'mp', 'quad', 'quad_id', 'plant_id', 'year')) %>%
  mutate(
    recruit_plot_valid = case_when(
      recruit == 1 & recruit_plot_valid == 1 ~ 1,
      TRUE ~ NA_real_))


# Vital-rate data --------------------------------------------------------------
df_su <- df %>%
  filter(
    size_t0 > 0,
    !is.na(survives),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(disturbance = as.numeric(disturbance))

df_gr <- df %>%
  filter(
    size_t0 > 0,
    survives == 1,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    !is.na(disturbance)) %>%
  mutate(disturbance = as.numeric(disturbance))

# Flowering is one binary process. Flower number is not modeled separately.
df_fl <- df %>%
  filter(
    !is.na(flower),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    disturbance = as.numeric(disturbance),
    flower = if_else(flower > 0, 1, 0))

# Fruiting probability is conditional on flowering.
df_fr <- df %>%
  filter(
    flower > 0,
    is.finite(fruit),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(
    disturbance = as.numeric(disturbance),
    fruiting = if_else(fruit > 0, 1, 0))

# Fruit number is conditional on fruiting.
df_fr_n <- df %>%
  filter(
    fruit > 0,
    is.finite(fruit),
    is.finite(logsize_t0),
    !is.na(disturbance)) %>%
  mutate(disturbance = as.numeric(disturbance))


# Fruit-to-recruit data --------------------------------------------------------
# Total fruits in year t predict valid recruits in year t + 1. The mean model
# is forced through zero and has no year or site random effects.
df_fr2re_quad <- df %>%
  group_by(site, quad_id, year) %>%
  summarise(
    fruits_t0 = sum(fruit, na.rm = TRUE),
    recruits_t0 = sum(recruit_plot_valid == 1, na.rm = TRUE),
    .groups = 'drop') %>%
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
    .groups = 'drop') %>%
  mutate(
    site = factor(site),
    logfruits_t0 = log1p(fruits_t0),
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

df_fr <- df_fr %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_fr_n <- df_fr_n %>%
  filter(!is.na(year), !(year %in% v_years_re))

df_fr2re <- df_fr2re %>%
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

mods_su <- list(mod_su_00, mod_su_0, mod_su_10, mod_su_1, 
                mod_su_20, mod_su_2, mod_su_30, mod_su_3)
mods_su_dAIC <- bbmle::AICctab(mods_su, weights = TRUE, sort = FALSE)$dAIC
mods_su_sorted <- order(mods_su_dAIC)

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
mod_gr_00 <- lm(logsize_t1 ~ 1,
               data = df_gr)

mod_gr_0 <- lm(logsize_t1 ~ disturbance,
               data = df_gr)

mod_gr_10 <- lm(logsize_t1 ~ logsize_t0,
               data = df_gr)

mod_gr_1 <- lm(logsize_t1 ~ logsize_t0 + disturbance,
               data = df_gr)

mod_gr_20 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2,
               data = df_gr)

mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + disturbance,
               data = df_gr)

mod_gr_30 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
               data = df_gr)

mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                 logsize_t0_3 + disturbance,
               data = df_gr)

mods_gr <- list(mod_gr_00, mod_gr_0, mod_gr_10, mod_gr_1, 
                mod_gr_20, mod_gr_2, mod_gr_30, mod_gr_3)
mods_gr_dAICc <- AICctab(mods_gr, weights = T, sort = F)$dAICc
mods_gr_sorted <- order(mods_gr_dAICc)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_best <- mods_gr_sorted[1]
} else {
  mod_gr_index_best <- v_mod_set_gr +1
}

v_mod_gr_index <- floor((mod_gr_index_best - 1) / 2)
mod_gr_best  <- mods_gr[[mod_gr_index_best]]
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
mod_fl_ranef <- coef(mod_fl_best)

summary(mod_fl_best)
mods_fl_dAICc


# Fruiting-probability model ---------------------------------------------------
mod_fr_00 <- glm(
  fruiting ~ 1,
  data = df_fr, family = binomial)

mod_fr_0 <- glm(
  fruiting ~ disturbance,
  data = df_fr, family = binomial)

mod_fr_10 <- glm(
  fruiting ~ logsize_t0,
  data = df_fr, family = binomial)

mod_fr_1 <- glm(
  fruiting ~ logsize_t0 + disturbance,
  data = df_fr, family = binomial)

mod_fr_20 <- glm(
  fruiting ~ logsize_t0 + logsize_t0_2,
  data = df_fr, family = binomial)

mod_fr_2 <- glm(
  fruiting ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_fr, family = binomial)

mod_fr_30 <- glm(
  fruiting ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fr, family = binomial)

mod_fr_3 <- glm(
  fruiting ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
  data = df_fr, family = binomial)

mods_fr <- list(
  mod_fr_00, mod_fr_0, mod_fr_10, mod_fr_1,
  mod_fr_20, mod_fr_2, mod_fr_30, mod_fr_3)

mods_fr_dAICc <- bbmle::AICctab(
  mods_fr, weights = TRUE, sort = FALSE)$dAICc
mods_fr_sorted <- order(mods_fr_dAICc)

if (length(v_mod_set_fr) == 0) {
  mod_fr_index_best <- mods_fr_sorted[1]
} else {
  mod_fr_index_best <- v_mod_set_fr + 1
}

v_mod_fr_index <- floor((mod_fr_index_best - 1) / 2)
mod_fr_best <- mods_fr[[mod_fr_index_best]]
mod_fr_ranef <- coef(mod_fr_best)

summary(mod_fr_best)
mods_fr_dAICc


# Fruit-number model -----------------------------------------------------------
mod_fr_n_00 <- MASS::glm.nb(
  fruit ~ 1,
  data = df_fr_n)

mod_fr_n_0 <- MASS::glm.nb(
  fruit ~ disturbance,
  data = df_fr_n)

mod_fr_n_10 <- MASS::glm.nb(
  fruit ~ logsize_t0,
  data = df_fr_n)

mod_fr_n_1 <- MASS::glm.nb(
  fruit ~ logsize_t0 + disturbance,
  data = df_fr_n)

mod_fr_n_20 <- MASS::glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2,
  data = df_fr_n)

mod_fr_n_2 <- MASS::glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + disturbance,
  data = df_fr_n)

mod_fr_n_30 <- MASS::glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fr_n)

mod_fr_n_3 <- MASS::glm.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + disturbance,
  data = df_fr_n)

mods_fr_n <- list(
  mod_fr_n_00, mod_fr_n_0, mod_fr_n_10, mod_fr_n_1,
  mod_fr_n_20, mod_fr_n_2, mod_fr_n_30, mod_fr_n_3)

mods_fr_n_dAICc <- bbmle::AICctab(
  mods_fr_n, weights = TRUE, sort = FALSE)$dAICc
mods_fr_n_sorted <- order(mods_fr_n_dAICc)

if (length(v_mod_set_fr_n) == 0) {
  mod_fr_n_index_best <- mods_fr_n_sorted[1]
} else {
  mod_fr_n_index_best <- v_mod_set_fr_n + 1
}

v_mod_fr_n_index <- floor((mod_fr_n_index_best - 1) / 2)
mod_fr_n_best <- mods_fr_n[[mod_fr_n_index_best]]
mod_fr_n_ranef <- coef(mod_fr_n_best)

summary(mod_fr_n_best)
mods_fr_n_dAICc


# Mean fruit-to-recruit model --------------------------------------------------
# Zero fruits predict zero recruits. No site or year random effects are fitted.
mod_fr2re_mean <- lm(
  logrecruits_t1 ~ 0 + logfruits_t0,
  data = df_fr2re)

summary(mod_fr2re_mean)


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

fig_fr_mean <- make_mean_plot(
  data = df_fr,
  model = mod_fr_best,
  response = 'fruiting',
  title = 'Mean fruiting probability given flowering',
  y_lab = 'Fruiting probability',
  probability = TRUE)

fig_fr_n_mean <- make_mean_plot(
  data = df_fr_n,
  model = mod_fr_n_best,
  response = 'fruit',
  title = 'Mean fruit number given fruiting',
  y_lab = 'Number of fruits')

fig_vital_rates_mean <-
  fig_su_mean + fig_gr_mean + fig_fl_mean +
  fig_fr_mean + fig_fr_n_mean +
  plot_layout(ncol = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

fig_vital_rates_mean


# Fruit-to-recruit plot --------------------------------------------------------
df_fr2re_pred <- tibble(
  logfruits_t0 = seq(
    0,
    max(df_fr2re$logfruits_t0, na.rm = TRUE),
    length.out = 200)) %>%
  mutate(
    logrecruits_t1 = predict(
      mod_fr2re_mean,
      newdata = .),
    fruits_t0 = expm1(logfruits_t0),
    recruits_t1 = pmax(0, expm1(logrecruits_t1)))

fig_fr2re_mean <- ggplot(
  df_fr2re,
  aes(x = fruits_t0, y = recruits_t1)) +
  geom_point(alpha = 0.4, size = 1.3) +
  geom_line(
    data = df_fr2re_pred,
    aes(x = fruits_t0, y = recruits_t1),
    linewidth = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1)) +
  theme_bw() +
  labs(
    title = 'Mean fruit-to-recruit relationship',
    subtitle = v_ggp_suffix,
    x = expression('Total fruits'[t]),
    y = expression('Valid recruits'[t+1]))

fig_fr2re_mean


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
  extract_terms(mod_fr_best, 'fr_'),
  extract_terms(mod_fr_n_best, 'frn_'),
  fr2re_b0 = 0,
  fr2re_b1 = unname(
    stats::coef(mod_fr2re_mean)[['logfruits_t0']]),
  fr2re_sigma = stats::sigma(mod_fr2re_mean),
  recr_mean = rc_sz$mean,
  recr_sd = rc_sz$sd,
  grow_var_a = unname(grow_var_coef[['a']]),
  grow_var_b = unname(grow_var_coef[['b']]),
  fruits_ref = mean(df_fr2re$fruits_t0, na.rm = TRUE),
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
# Flower count is not multiplied into this expression. The fruit-number model
# predicts total fruits for fruiting individuals.
fruits_x <- function(x, pars, disturbance = 0) {
  fl_x(x, pars, disturbance) *
    fr_x(x, pars, disturbance) *
    frn_x(x, pars, disturbance)
}


# Fruit-to-recruit conversion --------------------------------------------------
predict_recruits_from_fruits <- function(fruits, pars) {
  fruits_nonnegative <- pmax(fruits, 0)

  eta <- get_par(pars, 'fr2re_b0') +
    get_par(pars, 'fr2re_b1') * log1p(fruits_nonnegative)

  recruits <- expm1(
    eta + 0.5 * get_par(pars, 'fr2re_sigma')^2)

  ifelse(fruits_nonnegative <= 0, 0, pmax(recruits, 0))
}

recruits_per_fruit <- function(pars) {
  fruits_ref <- get_par(pars, 'fruits_ref', 0)

  if (!is.finite(fruits_ref) || fruits_ref <= 0) {
    return(0)
  }

  predict_recruits_from_fruits(fruits_ref, pars) / fruits_ref
}

rx <- function(x, pars, disturbance = 0) {
  fruits_x(
    x = x,
    pars = pars,
    disturbance = disturbance) *
    recruits_per_fruit(pars)
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


# Expected total fruits from the mean size distribution -----------------------
expected_mean_fruits <- function(
    pars, disturbance = 0, df_init = df) {
  n_obs <- make_initial_n_mean(pars, df_init)
  n <- as.integer(pars$mat_siz)
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n

  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  sum(
    fruits_x(y, pars, disturbance) * n_obs * h,
    na.rm = TRUE)
}

mean_fruits_no_fire <- expected_mean_fruits(
  pars = pars_mean,
  disturbance = 0)

mean_fruits_fire <- expected_mean_fruits(
  pars = pars_mean,
  disturbance = 1)


# Mean disturbance exposure ---------------------------------------------------
# Each sampled quadrat-year receives equal weight; plant-rich quadrats do not
# receive greater weight merely because they contain more plant records.
p_disturbance_mean <- df %>%
  filter(
    !is.na(year),
    !is.na(quad_id),
    !is.na(disturbance)) %>%
  distinct(year, quad_id, disturbance) %>%
  summarise(
    p_disturbance = mean(disturbance, na.rm = TRUE)) %>%
  pull(p_disturbance)

mean_fruits_regime <-
  (1 - p_disturbance_mean) * mean_fruits_no_fire +
  p_disturbance_mean * mean_fruits_fire


# State-specific and disturbance-regime kernels -------------------------------
pars_no_fire <- pars_mean
pars_no_fire$fruits_ref <- mean_fruits_no_fire

pars_fire <- pars_mean
pars_fire$fruits_ref <- mean_fruits_fire

# For the mixture kernel, both component kernels use the same expected total
# fruit production under the observed mean disturbance regime.
pars_regime <- pars_mean
pars_regime$fruits_ref <- mean_fruits_regime

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
    fruits_ref = c(
      mean_fruits_no_fire,
      mean_fruits_fire,
      mean_fruits_regime),
    recruits_ref = purrr::map2_dbl(
      fruits_ref,
      list(pars_no_fire, pars_fire, pars_regime),
      predict_recruits_from_fruits),
    recruits_per_fruit = if_else(
      fruits_ref > 0,
      recruits_ref / fruits_ref,
      0))

df_lambda_mean


# Observed annual population growth -------------------------------------------
# Counts are restricted to quadrats represented in consecutive years. Because
# quadrats containing no plant records cannot be identified from the individual
# table, transitions to a true zero count may remain unavailable.
df_counts_quad_year <- df %>%
  filter(size_t0 > 0, is.finite(logsize_t0)) %>%
  group_by(quad_id, year) %>%
  summarise(
    n_t0 = n_distinct(plant_id),
    .groups = 'drop') %>%
  arrange(quad_id, year) %>%
  group_by(quad_id) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n_t0),
    year_gap = year_t1 - year) %>%
  ungroup() %>%
  filter(year_gap == 1)

df_obs_growth_year <- df_counts_quad_year %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    n_quads = n(),
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
