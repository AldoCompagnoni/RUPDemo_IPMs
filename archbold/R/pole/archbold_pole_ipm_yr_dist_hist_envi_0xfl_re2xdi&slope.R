# IPM year-specific, environmental and historical, with fire 
# Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.08.14

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


# Recruitment data ------------------------------------------------------------
df_quad_year_sampled <- df %>%
  distinct(site, quad, year)

df_re <- df %>%
  group_by(site, quad, year) %>%
  summarise(
    disturbance = first(
      disturbance[!is.na(disturbance)], default = NA_real_),
    disturbance_prev = first(
      disturbance_prev[!is.na(disturbance_prev)], default = NA_real_),
    .groups = 'drop') %>%
  mutate(year_t1 = year + 1) %>%
  semi_join(
    df_quad_year_sampled %>%
      transmute(site, quad, year_t1 = year),
    by = c('site', 'quad', 'year_t1')) %>%
  left_join(
    df %>%
      filter(recruits == 1) %>%
      count(site, quad, year, name = 're_t1'),
    by = c('site', 'quad', 'year_t1' = 'year')) %>%
  mutate(re_t1 = replace_na(re_t1, 0L)) %>%
  filter(
    !is.na(year),
    !is.na(year_t1),
    !is.na(disturbance),
    !is.na(disturbance_prev),
    !(year %in% v_years_re),
    !(year_t1 %in% v_years_re)) %>%
  mutate(
    site = factor(site),
    year_t0 = factor(year),
    year_t1 = factor(year_t1))


# Recruit model -------------------------------------------------------
ctrl_re <- glmerControl(
  optimizer = 'bobyqa',
  optCtrl = list(maxfun = 2e5))

# Intercept only
mod_re_00 <- glmer.nb(
  re_t1 ~ 1 +
    (1 | year_t1),
  data = df_re,
  control = ctrl_re)


# Current disturbance only
mod_re_0 <- glmer.nb(
  re_t1 ~ disturbance +
    (1 | year_t1),
  data = df_re,
  control = ctrl_re)


# Current + previous disturbance
mod_re_01 <- glmer.nb(
  re_t1 ~ disturbance + disturbance_prev +
    (1 | year_t1),
  data = df_re,
  control = ctrl_re)


# Current fire effect varies among recruitment years
mod_re_yrfire <- glmer.nb(
  re_t1 ~ disturbance + disturbance_prev +
    (1 + disturbance || year_t1),
  data = df_re,
  control = ctrl_re)

mods_re <- list(
  mod_re_00,
  mod_re_0,
  mod_re_01,
  mod_re_yrfire)

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
fire_levels <- c(
  'No fire',
  'Current-year fire',
  'Previous-year fire',
  'Fire in both years')

df_re_plot <- df_re %>%
  mutate(
    fire_group = case_when(
      disturbance == 0 & disturbance_prev == 0 ~ 'No fire',
      disturbance == 1 & disturbance_prev == 0 ~ 'Current-year fire',
      disturbance == 0 & disturbance_prev == 1 ~ 'Previous-year fire',
      disturbance == 1 & disturbance_prev == 1 ~ 'Fire in both years'),
    fire_group = factor(fire_group, levels = fire_levels),
    re_t1_log = log1p(re_t1))

df_re_pred <- expand_grid(
  year_t1 = levels(df_re$year_t1),
  disturbance = c(0, 1),
  disturbance_prev = c(0, 1)) %>%
  mutate(
    year_t1 = factor(year_t1, levels = levels(df_re$year_t1)),
    fire_group = case_when(
      disturbance == 0 & disturbance_prev == 0 ~ 'No fire',
      disturbance == 1 & disturbance_prev == 0 ~ 'Current-year fire',
      disturbance == 0 & disturbance_prev == 1 ~ 'Previous-year fire',
      disturbance == 1 & disturbance_prev == 1 ~ 'Fire in both years'),
    fire_group = factor(fire_group, levels = fire_levels),
    pred_re = predict(
      mod_re_best, newdata = ., type = 'response',
      re.form = NULL, allow.new.levels = TRUE),
    pred_re_log = log1p(pred_re))

fig_re_years <- ggplot(
  df_re_plot,
  aes(x = fire_group, y = re_t1_log, color = fire_group)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1) +
  geom_point(
    data = df_re_pred,
    aes(y = pred_re_log),
    shape = 4, size = 3, stroke = 1) +
  scale_color_manual(
    values = c(
      'No fire' = 'black',
      'Current-year fire' = 'red',
      'Previous-year fire' = 'purple',
      'Fire in both years' = 'magenta')) +
  facet_wrap(~ year_t1, ncol = 4) +
  labs(
    title = 'Recruitment - year specific',
    subtitle = paste(
      v_ggp_suffix,
      '- crosses are model predictions'),
    x = 'Fire history',
    y = 'log(1 + recruits)',
    color = 'Fire history') +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 6),
    legend.position = 'bottom')

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
      coefficient = paste0(
        're_u0_', group_label, '_', rownames(ranef_matrix)),
      value = ranef_matrix[, '(Intercept)'])
  }
  
  if ('disturbance' %in% colnames(ranef_matrix)) {
    out[['disturbance']] <- data.frame(
      coefficient = paste0(
        're_ud_', group_label, '_', rownames(ranef_matrix)),
      value = ranef_matrix[, 'disturbance'])
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

make_term_map <- function(prefix, map) {
  out <- map
  names(out) <- paste0(prefix, names(out))
  out
}

su_term_map <- make_term_map('surv_', size_term_map)
gr_term_map <- make_term_map('grow_', size_term_map)

re_term_map <- list(
  re_b0 = c('(Intercept)'),
  re_bd = c('disturbance'),
  re_bdp = c('disturbance_prev'))


# Fixed parameters -------------------------------------------------------------
su_fe <- extract_fixed_pars(mod_su_best, su_term_map, fill_missing = TRUE)
gr_fe <- extract_fixed_pars(mod_gr_best, gr_term_map, fill_missing = TRUE)
re_fe <- extract_fixed_pars(mod_re_best, re_term_map, fill_missing = TRUE)


# Constants --------------------------------------------------------------------
gr_var_coef <- coef(mod_gr_var)

constants <- tibble::tribble(
  ~coefficient, ~value,
  'recr_sz', rc_sz$mean,
  'recr_sd', rc_sz$sd,
  'a', as.numeric(gr_var_coef[1]),
  'b', as.numeric(gr_var_coef[2]),
  'L', min(c(
    df_su$logvol_t0,
    df_gr$logvol_t0,
    df_recr$logvol_t0), na.rm = TRUE) - 0.1,
  'U', max(c(
    df_su$logvol_t0,
    df_gr$logvol_t0,
    df_recr$logvol_t0), na.rm = TRUE) + 0.1,
  'mat_siz', 200,
  'mod_su_index', v_mod_su_index,
  'mod_gr_index', v_mod_gr_index,
  'mod_re_index', v_mod_re_index) %>%
  mutate(coefficient = as.character(coefficient), value = as.numeric(value))

pars_cons <- bind_coef_rows(list(
  su_fe,
  gr_fe,
  re_fe,
  constants))

check_duplicate_coefficients(pars_cons, object_name = 'pars_cons')
pars_all_mean <- coef_df_to_list(pars_cons)
pars_mean <- pars_all_mean


# Year-varying parameters ------------------------------------------------------
su_out_yr <- extract_group_pars(
  mod_su_best, 'year', su_term_map, fill_missing = TRUE)

gr_out_yr <- extract_group_pars(
  mod_gr_best, 'year', gr_term_map, fill_missing = TRUE)

pars_var <- bind_coef_rows(list(
  su_out_yr,
  gr_out_yr))

check_duplicate_coefficients(pars_var, object_name = 'pars_var')
pars_all_year <- coef_df_to_list(pars_var)


# Recruitment-year random-effect deviations ---------------------------
pars_re_year_t1 <- extract_re_ranef_devs(
  mod_re_best, 'year_t1', 'year_t1') %>%
  mutate(coefficient = as.character(coefficient))
check_duplicate_coefficients(
  pars_re_year_t1, object_name = 'pars_re_year_t1')
pars_all_re_year_t1 <- coef_df_to_list(pars_re_year_t1)


# Disturbance transition probabilities ----------------------------------------
# Current disturbance refers to t -> t + 1 and previous disturbance to
# t - 1 -> t.

df_disturbance_regime <- df %>%
  filter(
    !is.na(site),
    !is.na(quad),
    !is.na(year),
    !is.na(disturbance),
    !is.na(disturbance_prev)) %>%
  distinct(
    site,
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


df_disturbance_transition <- df_disturbance_regime %>%
  group_by(disturbance_prev) %>%
  mutate(
    p_transition = n_quadrat_years /
      sum(n_quadrat_years)) %>%
  ungroup()

df_disturbance_transition


get_transition_probability <- function(
    disturbance_i,
    disturbance_prev_i) {
  
  df_disturbance_transition %>%
    filter(
      disturbance == disturbance_i,
      disturbance_prev == disturbance_prev_i) %>%
    pull(p_transition)
}


p_0_given_0 <- get_transition_probability(
  disturbance_i = 0,
  disturbance_prev_i = 0)

p_1_given_0 <- get_transition_probability(
  disturbance_i = 1,
  disturbance_prev_i = 0)

p_0_given_1 <- get_transition_probability(
  disturbance_i = 0,
  disturbance_prev_i = 1)

p_1_given_1 <- get_transition_probability(
  disturbance_i = 1,
  disturbance_prev_i = 1)


df_disturbance_transition


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
    
    pars$re_bd <- get_par(pars, 're_bd', 0) +
      get_par(
        pars_re_year_t1,
        paste0('re_ud_year_t1_', re_year_t1), 0)
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


re_total_ref <- function(
    pars,
    disturbance = 0,
    disturbance_prev = 0) {
  
  re_eta <- get_par(pars, 're_b0', 0) +
    get_par(pars, 're_bd', 0) * disturbance +
    get_par(pars, 're_bdp', 0) * disturbance_prev
  
  pmax(exp(re_eta), 0)
}

re_y_dist <- function(y, pars, h = NULL) {
  dens <- dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
  
  if (!is.null(h)) {
    dens <- dens / sum(dens * h)
  }
  
  dens
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
  
  re_total <- re_total_ref(
    pars = pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)
  
  R <- re_total *
    re_y_dist(
      y,
      pars,
      h = h)
  
  A_aug <- rbind(
    cbind(T, R),
    c(rep(0, n), 1))
  
  list(
    A_aug = A_aug,
    T = T,
    R = R,
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
  
  A <- kernel(
    pars = pars,
    disturbance = disturbance,
    disturbance_prev = disturbance_prev)$A_aug
  
  max(Mod(
    eigen(
      A,
      only.values = TRUE)$values))
}


# Environmental-history kernel ------------------------------------------------
kernel_env_hist <- function(pars) {
  
  A_00 <- kernel(
    pars = pars,
    disturbance = 0,
    disturbance_prev = 0)$A_aug
  
  A_10 <- kernel(
    pars = pars,
    disturbance = 1,
    disturbance_prev = 0)$A_aug
  
  A_01 <- kernel(
    pars = pars,
    disturbance = 0,
    disturbance_prev = 1)$A_aug
  
  A_11 <- kernel(
    pars = pars,
    disturbance = 1,
    disturbance_prev = 1)$A_aug
  
  A_env <- rbind(
    cbind(
      p_0_given_0 * A_00,
      p_0_given_1 * A_01),
    cbind(
      p_1_given_0 * A_10,
      p_1_given_1 * A_11))
  
  list(
    A_env = A_env,
    A_00 = A_00,
    A_10 = A_10,
    A_01 = A_01,
    A_11 = A_11)
}


# Environmental-history asymptotic lambda -------------------------------------
lambda_env_hist <- function(pars) {
  
  A <- kernel_env_hist(
    pars = pars)$A_env
  
  max(Mod(
    eigen(
      A,
      only.values = TRUE)$values))
}


# Mean IPMs --------------------------------------------------------------------
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 0)
lambda_ipm(pars_mean, disturbance = 0, disturbance_prev = 1)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 1)
lambda_ipm(pars_mean, disturbance = 1, disturbance_prev = 0)

lambda_env_hist(pars_mean)


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


# Year-specific environmental-history lambda -----------------------------------
lambda_env_hist_year <- function(
    year,
    re_year_t1 = NULL) {
  
  pars_i <- make_ipm_pars(
    pars_mean = pars_all_mean,
    pars_year = pars_all_year,
    pars_re_year_t1 = pars_all_re_year_t1,
    year = year,
    re_year_t1 = re_year_t1)
  
  lambda_env_hist(
    pars = pars_i)
}


# Year-specific disturbance scenarios ------------------------------------------
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
      disturbance_prev = 1)}),
  lambda_environmental_history = sapply(ipm_years, function(i) {
    lambda_env_hist_year(
      year = i)})
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
      lambda_consecutive_fire = 'Fire in both years',
      lambda_environmental_history = 'Environmental-history regime')) %>%
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
  
  K_i <- kernel(
    pars = pars_y,
    disturbance = disturbance_y,
    disturbance_prev = disturbance_prev_y)
  
  h <- K_i$h
  
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
        n_surv_growth_model = NA_real_,
        n_recr_model = NA_real_,
        n_proj_model = NA_real_,
        state_asym_lambda = NA_real_,
        proj_lambda = NA_real_))
  }
  
  # Existing individuals
  n_surv_growth <- K_i$T %*% n_obs
  
  # Direct quadrat-level recruitment
  n_recruits <- K_i$R
  
  # Historical one-year projection
  n_proj <- n_surv_growth + n_recruits
  
  n_surv_growth_total <- sum(n_surv_growth) * h
  n_recruits_total <- sum(n_recruits) * h
  n_projected <- sum(n_proj) * h
  
  state_asym_lambda <- max(Mod(
    eigen(
      K_i$A_aug,
      only.values = TRUE)$values))
  
  data.frame(
    year = yr,
    site = site_i,
    quad = quad_i,
    disturbance_num = disturbance_y,
    disturbance_prev_num = disturbance_prev_y,
    n_obs_model = n_initial,
    n_surv_growth_model = n_surv_growth_total,
    n_recr_model = n_recruits_total,
    n_proj_model = n_projected,
    state_asym_lambda = state_asym_lambda,
    proj_lambda = n_projected / n_initial)
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
    error_state_asymptotic_vs_obs =
      state_asym_lambda - obs_pgr,
    error_projected_vs_obs =
      proj_lambda - obs_pgr)


# Disturbance-state kernel used for each quadrat-year --------------------------
df_compare_quad <- df_compare_quad %>%
  mutate(
    kernel_state = case_when(
      disturbance_num == 0 & disturbance_prev_num == 0 ~ 'K_00',
      disturbance_num == 1 & disturbance_prev_num == 0 ~ 'K_10',
      disturbance_num == 0 & disturbance_prev_num == 1 ~ 'K_01',
      disturbance_num == 1 & disturbance_prev_num == 1 ~ 'K_11',
      TRUE ~ NA_character_))


df_kernel_state_use <- df_compare_quad %>%
  count(
    year,
    kernel_state,
    name = 'n_quadrat_transitions') %>%
  arrange(
    year,
    kernel_state)

df_kernel_state_use


# Site-year comparison from matched quadrats ----------------------------------
df_compare_site <- df_compare_quad %>%
  group_by(
    site,
    year) %>%
  summarise(
    n_t0 = sum(
      n_t0,
      na.rm = TRUE),
    n_t1 = sum(
      n_t1,
      na.rm = TRUE),
    n_obs_model = sum(
      n_obs_model,
      na.rm = TRUE),
    n_surv_growth_model = sum(
      n_surv_growth_model,
      na.rm = TRUE),
    n_recr_model = sum(
      n_recr_model,
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
    n_t0 = sum(
      n_t0,
      na.rm = TRUE),
    n_t1 = sum(
      n_t1,
      na.rm = TRUE),
    n_obs_model = sum(
      n_obs_model,
      na.rm = TRUE),
    n_surv_growth_model = sum(
      n_surv_growth_model,
      na.rm = TRUE),
    n_recr_model = sum(
      n_recr_model,
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
      asym_lambda = 'Environmental-history asymptotic lambda',
      proj_lambda =
        'Projected growth from observed size distributions + direct recruitment'))

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
      sqrt(
        mean(
          (asym_lambda - obs_pgr)^2,
          na.rm = TRUE)),
    rmse_projected_vs_obs =
      sqrt(
        mean(
          (proj_lambda - obs_pgr)^2,
          na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'statistic',
    values_to = 'value')

df_compare_summary

df_compare %>% print(n = 100)

# Historical projection components -------------------------------------------
df_projection_components <- df_compare %>%
  transmute(
    year,
    observed_t0 = n_t0,
    observed_t1 = n_t1,
    predicted_surv_growth = n_surv_growth_model,
    predicted_recruits = n_recr_model,
    predicted_t1 = n_proj_model,
    observed_growth = obs_pgr,
    projected_growth = proj_lambda)

df_projection_components %>%
  print(n = 100)


# Recruitment and survival-growth diagnostic ---------------------------------
df_recruits_quad_year <- df %>%
  group_by(site, quad, year) %>%
  summarise(
    obs_recruits = sum(recruits == 1, na.rm = TRUE),
    .groups = 'drop')

df_components_quad <- df_obs_pgr_quad %>%
  left_join(
    df_recruits_quad_year %>%
      rename(year_t1 = year),
    by = c('site', 'quad', 'year_t1')) %>%
  mutate(
    obs_recruits = replace_na(obs_recruits, 0),
    obs_nonrecruits = pmax(n_t1 - obs_recruits, 0)) %>%
  left_join(
    df_proj_quad %>%
      select(
        year,
        site,
        quad,
        disturbance_num,
        disturbance_prev_num,
        n_surv_growth_model,
        n_recr_model),
    by = c(
      'year',
      'site',
      'quad',
      'disturbance_num',
      'disturbance_prev_num'))

df_components_year <- df_components_quad %>%
  group_by(year) %>%
  summarise(
    obs_t1 = sum(n_t1),
    obs_recruits = sum(obs_recruits),
    pred_recruits = sum(n_recr_model),
    obs_nonrecruits = sum(obs_nonrecruits),
    pred_surv_growth = sum(n_surv_growth_model),
    recruit_error = pred_recruits - obs_recruits,
    surv_growth_error = pred_surv_growth - obs_nonrecruits,
    .groups = 'drop')

df_components_year %>%
  print(n = 100)

df_components_year %>%
  summarise(
    RMSE_recruitment =
      sqrt(mean(recruit_error^2)),
    RMSE_surv_growth =
      sqrt(mean(surv_growth_error^2)),
    MAE_recruitment =
      mean(abs(recruit_error)),
    MAE_surv_growth =
      mean(abs(surv_growth_error)))


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
  print(
    n = 100,
    width = Inf)


# Yearly lambda comparison diagnostic -----------------------------------------
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
  print(
    n = 100,
    width = Inf)


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
