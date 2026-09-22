# IPM mean
# McGraw 2017 - Panax quinquefolius

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.09.22

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
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data ------------------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0(v_head, '_', v_sp_abb, '_df_workdata.csv'))) %>% 
  filter(state == v_states)

required_cols <- c(
  'state', 'population', 'id', 'year', 'year_t1', 'size_t0', 'size_t1',
  'logsize_t0', 'logsize_t1', 'logsize_t0_2', 'logsize_t0_3', 'survives',
  'recruit', 'consecutive', 'reliable_survival', 'persistence_t0')

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste(
    'Missing columns in working data:',
    paste(missing_cols, collapse = ', ')))}


# Survival --------------------------------------------------------------------
df_su <- df %>%
  filter(
    reliable_survival,
    !is.na(survives),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  dplyr::select(
    state, population, id, year, size_t0, survives,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Survival model --------------------------------------------------------------
mod_su_0 <- glm(
  survives ~ 1,
  data = df_su, family = 'binomial')

mod_su_1 <- glm(
  survives ~ logsize_t0,
  data = df_su, family = 'binomial')

mod_su_2 <- glm(
  survives ~ logsize_t0 + logsize_t0_2,
  data = df_su, family = 'binomial')

mod_su_3 <- glm(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_su, family = 'binomial')

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

mod_su_bestfit <- mods_su[[mod_su_index_bestfit]]

mod_su_bestfit
mods_su_dAICc


# Survival plot ---------------------------------------------------------------
df_su_newdata <- tibble(
  logsize_t0 = seq(
    min(df_su$logsize_t0),
    max(df_su$logsize_t0),
    length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_su_newdata$predicted <- predict(
  mod_su_bestfit,
  newdata = df_su_newdata,
  type = 'response')

fig_su_line <- ggplot(
  df_su,
  aes(x = logsize_t0, y = survives)) +
  geom_jitter(
    height = 0.05, width = 0,
    alpha = 0.3) +
  geom_line(
    data = df_su_newdata,
    aes(y = predicted),
    linewidth = 1) +
  labs(
    title = 'Survival probability by size',
    subtitle = v_ggp_suffix,
    x = expression('log(leaf area)'[t0]),
    y = 'Probability of survival') +
  theme_bw()


df_su_bindata <- plot_binned_prop(
  df_su, 10, logsize_t0, survives)

fig_su_bin <- ggplot() +
  geom_point(
    data = df_su_bindata,
    aes(x = logsize_t0, y = survives)) +
  geom_errorbar(
    data = df_su_bindata,
    aes(
      x = logsize_t0,
      ymin = lwr,
      ymax = upr),
    width = 0.1) +
  geom_line(
    data = df_su_newdata,
    aes(x = logsize_t0, y = predicted),
    linewidth = 1.2) +
  labs(
    x = expression('log(leaf area)'[t0]),
    y = '') +
  theme_bw() +
  ylim(0, 1)

fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su


# Growth data -----------------------------------------------------------------
df_gr <- df %>%
  filter(
    consecutive,
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1)) %>%
  dplyr::select(
    state, population, id, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)


# Growth model ----------------------------------------------------------------
mod_gr_0 <- lm(
  logsize_t1 ~ 1,
  data = df_gr)

mod_gr_1 <- lm(
  logsize_t1 ~ logsize_t0,
  data = df_gr)

mod_gr_2 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2,
  data = df_gr)

mod_gr_3 <- lm(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_gr)

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

mod_gr_bestfit <- mods_gr[[mod_gr_index_bestfit]]

mod_gr_bestfit
mods_gr_dAICc


# Growth plot -----------------------------------------------------------------
df_gr_newdata <- tibble(
  logsize_t0 = seq(
    min(df_gr$logsize_t0),
    max(df_gr$logsize_t0),
    length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_gr_newdata <- cbind(
  df_gr_newdata,
  as.data.frame(
    predict(
      mod_gr_bestfit,
      newdata = df_gr_newdata,
      interval = 'confidence')))

fig_gr <- ggplot(
  df_gr,
  aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(
    data = df_gr_newdata,
    aes(y = fit),
    linewidth = 1,
    color = "lightgreen") +
  geom_ribbon(
    data = df_gr_newdata,
    aes(y = fit, ymin = lwr, ymax = upr),
    alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    title = 'Growth prediction',
    subtitle = v_ggp_suffix,
    x = expression('log(leaf area)'[t0]),
    y = expression('log(leaf area)'[t1])) +
  theme_bw()

fig_gr


# Growth variance -------------------------------------------------------------
mod_gr_x <- fitted(mod_gr_bestfit)
mod_gr_y <- resid(mod_gr_bestfit)^2

mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x),
  start = list(a = 1, b = 0),
  control = nls.control(
    maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Recruitment -----------------------------------------------------------------
# Constant mean recruitment is estimated as recruits per established plant.
# Parent abundance is measured in year t0 and recruits enter in year t1.

df_re_parent <- df %>%
  filter(
    !is.na(persistence_t0),
    persistence_t0 != 'DEAD',
    size_t0 > 0) %>%
  group_by(population, year) %>%
  summarise(
    n_parents = n_distinct(id),
    .groups = 'drop') %>%
  mutate(year_t1 = year + 1L)

# Years in which a population was actually sampled.
df_re_sampled <- df %>%
  distinct(population, year) %>%
  rename(year_t1 = year)

# New seedlings observed in year t1.
df_re_count <- df %>%
  filter(recruit == 1) %>%
  count(population, year, name = 'nr_recruits') %>%
  rename(year_t1 = year)

# Keep only population-year transitions with sampling at both t0 and t1.
df_re <- df_re_parent %>%
  inner_join(
    df_re_sampled,
    by = c('population', 'year_t1')) %>%
  left_join(
    df_re_count,
    by = c('population', 'year_t1')) %>%
  mutate(
    nr_recruits = replace_na(nr_recruits, 0L))

# Negative-binomial model with exposure offset. The intercept is the constant
# mean number of recruits produced per established plant per year.
mod_re <- MASS::glm.nb(
  nr_recruits ~ 1 + offset(log(n_parents)),
  data = df_re)

mod_re

df_re$pred_recruits <- predict(
  mod_re,
  newdata = df_re,
  type = 'response')

fecu_mean <- exp(unname(coef(mod_re)[['(Intercept)']]))

# Recruit size distribution at first observation.
df_re_size <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0))


# Recruitment plot ------------------------------------------------------------
fig_re <- ggplot(
  df_re,
  aes(x = n_parents, y = nr_recruits)) +
  geom_point(alpha = 0.6) +
  geom_point(
    aes(y = pred_recruits),
    shape = 1) +
  labs(
    title = 'Mean recruitment - Predicted recruits',
    subtitle = v_ggp_suffix,
    x = expression('Established individuals'[t0]),
    y = expression('New seedlings'[t1])) +
  theme_bw()

fig_re


# Exporting parameter estimates -----------------------------------------------
# Growth
pars_gr <- tibble(
  coefficient = names(coef(mod_gr_bestfit)),
  value = unname(coef(mod_gr_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))

pars_gr_var <- tibble(
  coefficient = names(coef(mod_gr_var)),
  value = unname(coef(mod_gr_var)))

pars_gr <- bind_rows(pars_gr, pars_gr_var)

# Survival
pars_su <- tibble(
  coefficient = names(coef(mod_su_bestfit)),
  value = unname(coef(mod_su_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))

# Recruitment and mesh limits
pars_other <- tibble(
  coefficient = c(
    'recr_sz', 'recr_sd', 'max_siz', 'min_siz', 'fecu_b0'),
  value = c(
    mean(df_re_size$logsize_t0, na.rm = TRUE),
    sd(df_re_size$logsize_t0, na.rm = TRUE),
    max(df_gr$logsize_t0, na.rm = TRUE),
    min(df_gr$logsize_t0, na.rm = TRUE),
    fecu_mean))

write.csv(
  pars_gr,
  file.path(
    dir_data,
    paste0(v_script_prefix, '_', v_sp_abb, '_grow_pars_mean.csv')),
  row.names = FALSE)

write.csv(
  pars_su,
  file.path(
    dir_data,
    paste0(v_script_prefix, '_', v_sp_abb, '_surv_pars_mean.csv')),
  row.names = FALSE)

write.csv(
  pars_other,
  file.path(
    dir_data,
    paste0(v_script_prefix, '_', v_sp_abb, '_other_pars_mean.csv')),
  row.names = FALSE)


# Building the IPM from scratch -----------------------------------------------
extr_value <- function(x, field) {
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix = v_script_prefix,
  species = v_species,
  surv_b0 = extr_value(pars_su, 'b0'),
  surv_b1 = extr_value(pars_su, 'logsize_t0'),
  surv_b2 = extr_value(pars_su, 'logsize_t0_2'),
  surv_b3 = extr_value(pars_su, 'logsize_t0_3'),
  grow_b0 = extr_value(pars_gr, 'b0'),
  grow_b1 = extr_value(pars_gr, 'logsize_t0'),
  grow_b2 = extr_value(pars_gr, 'logsize_t0_2'),
  grow_b3 = extr_value(pars_gr, 'logsize_t0_3'),
  a = extr_value(pars_gr, 'a'),
  b = extr_value(pars_gr, 'b'),
  fecu_b0 = extr_value(pars_other, 'fecu_b0'),
  recr_sz = extr_value(pars_other, 'recr_sz'),
  recr_sd = extr_value(pars_other, 'recr_sd'),
  L = extr_value(pars_other, 'min_siz'),
  U = extr_value(pars_other, 'max_siz'),
  mat_siz = 200,
  mod_gr_index = v_mod_gr_index,
  mod_su_index = v_mod_su_index))

write.csv(
  as.data.frame(pars),
  file.path(
    dir_data,
    paste0(v_script_prefix, '_', v_sp_abb, '_pars.csv')),
  row.names = FALSE)


# IPM functions ---------------------------------------------------------------
grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * x))
}

# Growth from size x to size y.
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- 0
  for (i in 0:num_pars) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i}
  }
  sd_value <- grow_sd(x, pars)
  dnorm(y, mean = mean_value, sd = sd_value)
}

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Survival of an x-sized individual to t1.
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('surv_b', i)
      if (!is.null(pars[[param_name]])) {
        survival_value <- survival_value + pars[[param_name]] * x^i}
    }
  }
  inv_logit(survival_value)
}

# Survival-growth transition.
pxy <- function(x, y, pars) {
  sx(x, pars) * gxy(x, y, pars)
}

# Constant mean fecundity distributed across recruit sizes.
fy <- function(y, pars, h) {
  n_recr <- pars$fecu_b0
  recr_sd <- max(h / 10, pars$recr_sd, na.rm = TRUE)
  recr_y <- dnorm(y, pars$recr_sz, recr_sd) * h
  recr_y <- recr_y / sum(recr_y)
  n_recr * recr_y
}


# Kernel ----------------------------------------------------------------------
kernel <- function(pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  Fmat <- matrix(0, n, n)
  Fmat[] <- matrix(fy(y, pars, h), n, n)

  Smat <- sx(y, pars)

  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h

  Tmat <- matrix(0, n, n)

  for (i in seq_len(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }

  for (i in ((n / 2) + 1):n) {
    Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    Tmat[, i] <- Gmat[, i] * Smat[i]
  }

  k_yx <- Fmat + Tmat

  list(
    k_yx = k_yx,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    meshpts = y)
}

lambda_ipm <- function(i) {
  Re(eigen(kernel(i)$k_yx)$values[1])
}


# Mean IPM population growth rate ---------------------------------------------
lam_mean <- lambda_ipm(pars)
lam_mean


# Observed population growth rate ---------------------------------------------
# Count living individuals in each sampled population-year.
df_pop_n <- df %>%
  filter(
    !is.na(persistence_t0),
    persistence_t0 != 'DEAD') %>%
  group_by(population, year) %>%
  summarise(
    n = n_distinct(id),
    .groups = 'drop')

pop_counts_t0 <- df_pop_n %>%
  transmute(
    population,
    year = year + 1L,
    n_t0 = n)

pop_counts_t1 <- df_pop_n %>%
  transmute(
    population,
    year,
    n_t1 = n)

pop_counts <- inner_join(
  pop_counts_t0,
  pop_counts_t1,
  by = c('population', 'year')) %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0),
    n_t1 = sum(n_t1),
    .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 / n_t0)

lam_mean_count <- exp(mean(log(pop_counts$obs_pgr), na.rm = TRUE))
lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)

lam_mean_count
lam_mean_overall


# Building the IPM with ipmr --------------------------------------------------
proto_ipm_p <- init_ipm(
  sim_gen = 'simple',
  di_dd = 'di',
  det_stoch = 'det') %>%
  define_kernel(
    name = 'P',
    family = 'CC',
    formula = s * g,
    s = plogis(
      surv_b0 +
        (if (mod_su_index >= 1) surv_b1 * size_1 else 0) +
        (if (mod_su_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3 * size_1^3 else 0)),
    mu_g = grow_b0 +
      (if (mod_gr_index >= 1) grow_b1 * size_1 else 0) +
      (if (mod_gr_index >= 2) grow_b2 * size_1^2 else 0) +
      (if (mod_gr_index >= 3) grow_b3 * size_1^3 else 0),
    g = dnorm(size_2, mu_g, grow_sig),
    grow_sig = sqrt(a * exp(b * size_1)),
    data_list = pars,
    states = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions(
      fun = 'norm', target = 'g')) %>%
  define_kernel(
    name = 'F',
    family = 'CC',
    formula = fecu_b0 * r_d,
    r_d = dnorm(size_2, recr_sz, recr_sd),
    data_list = pars,
    states = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'r_d')) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c('P', 'F'),
      int_rule = rep('midpoint', 2),
      state_start = rep('size', 2),
      state_end = rep('size', 2))) %>%
  define_domains(
    size = c(
      pars$L,
      pars$U,
      pars$mat_siz)) %>%
  define_pop_state(
    n_size = rep(
      1 / pars$mat_siz,
      pars$mat_siz))

ipmr_p <- make_ipm(
  proto_ipm = proto_ipm_p,
  iterations = 200)

lam_mean_ipmr <- lambda(ipmr_p)

lam_out <- data.frame(
  coefficient = names(lam_mean_ipmr),
  value = lam_mean_ipmr)


lam_out_wide <- as.list(
  pivot_wider(
    lam_out,
    names_from = 'coefficient',
    values_from = 'value'))
