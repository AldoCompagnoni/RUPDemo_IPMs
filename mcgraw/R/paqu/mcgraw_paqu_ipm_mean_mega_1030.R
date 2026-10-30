# IPM mean; entry-size memory and at most two dormant censuses
# McGraw 2017 - Panax quinquefolius

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.10.29

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
v_output_suffix <- "_mega_1030"
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
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
v_mod_set_do   <- c()
v_mod_set_fl   <- c()
v_mod_set_se   <- c()
v_mod_set_se_n <- c()


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
  filter(state == v_states)


# Survival --------------------------------------------------------------------
df_su <- df %>%
  filter(
    reliable_demography,
    dormant_t0 == 0,
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
    reliable_demography,
    survives == 1,
    dormant_t0 == 0,
    dormant_t1 == 0,
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


# Recruit size distribution ---------------------------------------------------
df_re_size <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0))


# Dormancy entry data ----------------------------------------------------------
df_do <- df %>%
  filter(
    reliable_demography,
    dormant_t0 == 0,
    survives == 1,
    !is.na(dormant_t1),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  mutate(enter_dormancy = dormant_t1) %>%
  dplyr::select(
    id, year, size_t0, enter_dormancy,
    logsize_t0, logsize_t0_2, logsize_t0_3)

df_do %>%
  count(enter_dormancy)


# Dormancy entry model ---------------------------------------------------------
mod_do_0 <- glm(
  enter_dormancy ~ 1,
  data = df_do, family = "binomial")

mod_do_1 <- glm(
  enter_dormancy ~ logsize_t0,
  data = df_do, family = "binomial")

mod_do_2 <- glm(
  enter_dormancy ~ logsize_t0 + logsize_t0_2,
  data = df_do, family = "binomial")

mod_do_3 <- glm(
  enter_dormancy ~
    logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_do, family = "binomial")

mods_do <- list(
  mod_do_0, mod_do_1, mod_do_2, mod_do_3)

mods_do_dAICc <- AICctab(
  mods_do, weights = TRUE, sort = FALSE)$dAICc

mods_do_sorted <- order(mods_do_dAICc)

if (length(v_mod_set_do) == 0) {
  mod_do_index_bestfit <- mods_do_sorted[1]
  v_mod_do_index <- mod_do_index_bestfit - 1
} else {
  mod_do_index_bestfit <- v_mod_set_do + 1
  v_mod_do_index <- v_mod_set_do
}

mod_do_bestfit <- mods_do[[mod_do_index_bestfit]]

mod_do_bestfit
mods_do_dAICc

# Dormancy entry plot ----------------------------------------------------------
df_do_newdata <- df_do %>%
  summarise(
    x = list(seq(
      min(logsize_t0),
      max(logsize_t0),
      length.out = 100)),
    .groups = "drop") %>%
  unnest(x) %>%
  mutate(
    logsize_t0 = x,
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  mutate(
    predicted = predict(
      mod_do_bestfit,
      newdata = .,
      type = "response"))

fig_do_line <- ggplot(
  df_do,
  aes(x = logsize_t0, y = enter_dormancy)) +
  geom_jitter(
    height = 0.05, width = 0,
    alpha = 0.3) +
  geom_line(
    data = df_do_newdata,
    aes(y = predicted),
    linewidth = 1) +
  labs(
    title = "Dormancy probability by size",
    subtitle = v_ggp_suffix,
    x = expression("log(diameter)"[t0]),
    y = "Probability of entering dormancy") +
  theme_bw()


df_do_bindata <- plot_binned_prop(
  df_do, 10, logsize_t0, enter_dormancy)

df_do_pred <- df_do %>%
  reframe(
    logsize_t0 = seq(
      min(logsize_t0),
      max(logsize_t0),
      length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  mutate(
    enter_dormancy = predict(
      mod_do_bestfit,
      newdata = .,
      type = "response"))

fig_do_bin <- ggplot() +
  geom_point(
    data = df_do_bindata,
    aes(x = logsize_t0, y = enter_dormancy)) +
  geom_errorbar(
    data = df_do_bindata,
    aes(
      x = logsize_t0,
      ymin = lwr,
      ymax = upr),
    width = 0.1) +
  geom_line(
    data = df_do_pred,
    aes(x = logsize_t0, y = enter_dormancy),
    linewidth = 1.2) +
  labs(
    x = expression("log(diameter)"[t0]),
    y = "") +
  theme_bw() +
  ylim(0, 1)

fig_do <- fig_do_line + fig_do_bin + plot_layout()
fig_do


# Dormant survival and reactivation  ------------------------------------------
# An active census, a gap, or an unconfirmed transition starts a new episode.
df_do_history <- df %>%
  group_by(state, population, id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    dormancy_spell = cumsum(
      coalesce(dormant_t0 != 1, TRUE) |
        coalesce(year != lag(year_t1), TRUE) |
        !coalesce(lag(
          reliable_demography & survives == 1 & dormant_t1 == 1 &
            year_t1 == year + 1L), FALSE))) %>%
  ungroup()

# Identify entry even when the last active size was not measured.
# year_entry is the last active census; the first dormant census is +1.
df_do_entry <- df_do_history %>%
  filter(
    reliable_demography, dormant_t0 == 0, dormant_t1 == 1,
    survives == 1, year_t1 == year + 1L) %>%
  transmute(state, population, id, dormancy_spell, year_entry = year)

df_do_history <- df_do_history %>%
  left_join(
    df_do_entry, by = c("state", "population", "id", "dormancy_spell")) %>%
  mutate(dormancy_year = if_else(
    dormant_t0 == 1, as.integer(year - year_entry), NA_integer_))

df_ra <- df_do_history %>%
  filter(
    reliable_demography,
    dormant_t0 == 1,
    year_t1 == year + 1L,
    !is.na(survives))

# Share one annual dormant-survival estimate across both dormant stages.
mod_dorm_su <- glm(
  survives ~ 1,
  data = df_ra,
  family = "binomial")

p_dorm_survival <- predict(
  mod_dorm_su,
  newdata = data.frame(x = 1),
  type = "response")[1]


df_ra <- df_ra %>%
  filter(
    survives == 1,
    !is.na(dormant_t1)) %>%
  mutate(
    reactivate = as.integer(dormant_t1 == 0))

# Estimate first-year return from all surviving D1 transitions.
# Surviving D2 plants must return under the imposed two-census limit.
df_ra_1 <- df_ra %>%
  filter(dormancy_year == 1L)

mod_ra <- glm(
  reactivate ~ 1,
  data = df_ra_1,
  family = "binomial")

p_reactivate_1 <- predict(
  mod_ra,
  newdata = data.frame(x = 1),
  type = "response")[1]

# Diagnose observed spells extending beyond the imposed two-census limit.
# Include a confirmed D2 -> D3 even if no later source census is available.
df_dorm_over2 <- df_do_history %>%
  filter(
    dormant_t0 == 1,
    dormancy_year > 2L |
      (dormancy_year == 2L & reliable_demography & survives == 1 &
         dormant_t1 == 1 & year_t1 == year + 1L)) %>%
  distinct(state, population, id, dormancy_spell)

df_ra %>%
  count(dormancy_year, reactivate, name = "n_transitions")

nrow(df_dorm_over2)


df_ra_size <- df_do_history %>%
  filter(
    reliable_demography,
    dormant_t0 == 1,
    survives == 1,
    dormant_t1 == 0,
    year_t1 == year + 1L,
    dormancy_year %in% c(1L, 2L),
    size_t1 > 0) %>%
  mutate(
    logsize_reactivate = log(size_t1))


# Dormancy entry vs reactivation size -----------------------------------------
# Fit return size to completed one- or two-census dormancy episodes.
# Longer observed spells are counted above, not relabelled as two years.

df_do_entry_size <- df_do_history %>%
  filter(
    reliable_demography,
    dormant_t0 == 0,
    dormant_t1 == 1,
    survives == 1,
    year_t1 == year + 1L,
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  transmute(
    state, population, id, dormancy_spell,
    year_entry = year, logsize_entry = logsize_t0)

df_ra_entry_size <- df_ra_size %>%
  filter(year_t1 == year + 1L, is.finite(logsize_reactivate)) %>%
  transmute(
    state, population, id, dormancy_spell, year_ra = year,
    dormancy_year, logsize_reactivate) %>%
  inner_join(
    df_do_entry_size,
    by = c("state", "population", "id", "dormancy_spell"))

mod_ra_size <- lm(
  logsize_reactivate ~ logsize_entry,
  data = df_ra_entry_size)

summary(mod_ra_size)

df_ra_entry_size %>%
  summarise(
    n = n(),
    correlation = cor(logsize_entry, logsize_reactivate))

fig_ra_size <- ggplot(
  df_ra_entry_size,
  aes(x = logsize_entry, y = logsize_reactivate)) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, linewidth = 1) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed") +
  labs(
    title = "Dormancy entry vs reactivation size",
    subtitle = v_ggp_suffix,
    x = "log(leaf area) before dormancy",
    y = "log(leaf area) at reactivation") +
  theme_bw()

fig_ra_size


# Conditional reactivation size: mean depends on entry size.
# Use residual SD, not the marginal SD of all reactivation sizes.
react_b0 <- unname(coef(mod_ra_size)[["(Intercept)"]])
react_b1 <- unname(coef(mod_ra_size)[["logsize_entry"]])
react_sd <- sigma(mod_ra_size)


# Flower data -----------------------------------------------------------------
df_fl <- df %>%
  filter(
    !is.na(flower),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  dplyr::select(
    state, population, id, year, size_t0, flower,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Flower model ----------------------------------------------------------------
mod_fl_0 <- glm(
  flower ~ 1,
  data = df_fl,
  family = "binomial")

mod_fl_1 <- glm(
  flower ~ logsize_t0,
  data = df_fl,
  family = "binomial")

mod_fl_2 <- glm(
  flower ~ logsize_t0 + logsize_t0_2,
  data = df_fl,
  family = "binomial")

mod_fl_3 <- glm(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fl,
  family = "binomial")

mods_fl <- list(
  mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)

mods_fl_dAICc <- AICctab(
  mods_fl,
  weights = TRUE,
  sort = FALSE)$dAICc

mods_fl_sorted <- order(mods_fl_dAICc)

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index <- mod_fl_index_bestfit - 1
} else {
  mod_fl_index_bestfit <- v_mod_set_fl + 1
  v_mod_fl_index <- v_mod_set_fl
}

mod_fl_bestfit <- mods_fl[[mod_fl_index_bestfit]]

mod_fl_bestfit
mods_fl_dAICc


# Flowering probability plot --------------------------------------------------
df_fl_newdata <- tibble(
  logsize_t0 = seq(
    min(df_fl$logsize_t0),
    max(df_fl$logsize_t0),
    length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_fl_newdata$flower <- predict(
  mod_fl_bestfit,
  newdata = df_fl_newdata,
  type = "response")

df_fl_bindata <- plot_binned_prop(
  df_fl, 10, logsize_t0, flower)

fig_fl_line <- ggplot(
  df_fl,
  aes(x = logsize_t0, y = flower)) +
  geom_jitter(
    height = 0.05,
    width = 0,
    alpha = 0.3) +
  geom_line(
    data = df_fl_newdata,
    aes(y = flower),
    linewidth = 1) +
  labs(
    title = "Flowering probability by size",
    subtitle = v_ggp_suffix,
    x = expression("log(leaf area)"[t0]),
    y = "Probability of flowering") +
  theme_bw()

fig_fl_bin <- ggplot() +
  geom_point(
    data = df_fl_bindata,
    aes(x = logsize_t0, y = flower)) +
  geom_errorbar(
    data = df_fl_bindata,
    aes(
      x = logsize_t0,
      ymin = lwr,
      ymax = upr),
    width = 0.1) +
  geom_line(
    data = df_fl_newdata,
    aes(x = logsize_t0, y = flower),
    linewidth = 1.2) +
  labs(
    x = expression("log(leaf area)"[t0]),
    y = "") +
  theme_bw() +
  ylim(0, 1)

fig_fl <- fig_fl_line + fig_fl_bin + plot_layout()
fig_fl


# Seed production data --------------------------------------------------------
# Conditional on flowering.
# seed_prod = 1 indicates production of at least one seed.

df_se <- df %>%
  filter(
    flower == 1,
    !is.na(seed_nr),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  mutate(
    seed_prod = as.integer(seed_nr > 0)) %>%
  dplyr::select(
    state, population, id, year, size_t0,
    flower, seed_prod, seed_nr,
    logsize_t0, logsize_t0_2, logsize_t0_3)


# Seed production probability model ------------------------------------------
mod_se_0 <- glm(
  seed_prod ~ 1,
  data = df_se,
  family = "binomial")

mod_se_1 <- glm(
  seed_prod ~ logsize_t0,
  data = df_se,
  family = "binomial")

mod_se_2 <- glm(
  seed_prod ~ logsize_t0 + logsize_t0_2,
  data = df_se,
  family = "binomial")

mod_se_3 <- glm(
  seed_prod ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_se,
  family = "binomial")

mods_se <- list(
  mod_se_0, mod_se_1, mod_se_2, mod_se_3)

mods_se_dAICc <- AICctab(
  mods_se,
  weights = TRUE,
  sort = FALSE)$dAICc

mods_se_sorted <- order(mods_se_dAICc)

if (length(v_mod_set_se) == 0) {
  mod_se_index_bestfit <- mods_se_sorted[1]
  v_mod_se_index <- mod_se_index_bestfit - 1
} else {
  mod_se_index_bestfit <- v_mod_set_se + 1
  v_mod_se_index <- v_mod_set_se
}

mod_se_bestfit <- mods_se[[mod_se_index_bestfit]]

mod_se_bestfit
mods_se_dAICc


# Seed production probability plot -------------------------------------------
df_se_newdata <- tibble(
  logsize_t0 = seq(
    min(df_se$logsize_t0),
    max(df_se$logsize_t0),
    length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_se_newdata$seed_prod <- predict(
  mod_se_bestfit,
  newdata = df_se_newdata,
  type = "response")

df_se_bindata <- plot_binned_prop(
  df_se, 10, logsize_t0, seed_prod)

fig_se_line <- ggplot(
  df_se,
  aes(x = logsize_t0, y = seed_prod)) +
  geom_jitter(
    height = 0.05,
    width = 0,
    alpha = 0.3) +
  geom_line(
    data = df_se_newdata,
    aes(y = seed_prod),
    linewidth = 1) +
  labs(
    title = "Seed production probability by size",
    subtitle = v_ggp_suffix,
    x = expression("log(leaf area)"[t0]),
    y = "Probability of producing seeds") +
  theme_bw()

fig_se_bin <- ggplot() +
  geom_point(
    data = df_se_bindata,
    aes(x = logsize_t0, y = seed_prod)) +
  geom_errorbar(
    data = df_se_bindata,
    aes(
      x = logsize_t0,
      ymin = lwr,
      ymax = upr),
    width = 0.1) +
  geom_line(
    data = df_se_newdata,
    aes(x = logsize_t0, y = seed_prod),
    linewidth = 1.2) +
  labs(
    x = expression("log(leaf area)"[t0]),
    y = "") +
  theme_bw() +
  ylim(0, 1)

fig_se <- fig_se_line + fig_se_bin + plot_layout()
fig_se


# Seed number data ------------------------------------------------------------
# Conditional on flowering and positive seed production.

df_se_n <- df_se %>%
  filter(
    seed_prod == 1,
    seed_nr > 0)


# Seed number model -----------------------------------------------------------
mod_se_n_0 <- MASS::glm.nb(
  seed_nr ~ 1,
  data = df_se_n)

mod_se_n_1 <- MASS::glm.nb(
  seed_nr ~ logsize_t0,
  data = df_se_n)

mod_se_n_2 <- MASS::glm.nb(
  seed_nr ~ logsize_t0 + logsize_t0_2,
  data = df_se_n)

mod_se_n_3 <- MASS::glm.nb(
  seed_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_se_n)

mods_se_n <- list(
  mod_se_n_0,
  mod_se_n_1,
  mod_se_n_2,
  mod_se_n_3)

mods_se_n_dAICc <- AICctab(
  mods_se_n,
  weights = TRUE,
  sort = FALSE)$dAICc

mods_se_n_sorted <- order(mods_se_n_dAICc)

if (length(v_mod_set_se_n) == 0) {
  mod_se_n_index_bestfit <- mods_se_n_sorted[1]
  v_mod_se_n_index <- mod_se_n_index_bestfit - 1
} else {
  mod_se_n_index_bestfit <- v_mod_set_se_n + 1
  v_mod_se_n_index <- v_mod_set_se_n
}

mod_se_n_bestfit <- mods_se_n[[mod_se_n_index_bestfit]]

mod_se_n_bestfit
mods_se_n_dAICc


# Seed number plot ------------------------------------------------------------
df_se_n_pred <- tibble(
  logsize_t0 = seq(
    min(df_se_n$logsize_t0),
    max(df_se_n$logsize_t0),
    length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_se_n_pred$seed_nr <- predict(
  mod_se_n_bestfit,
  newdata = df_se_n_pred,
  type = "response")

df_se_n_binned <- df_se_n %>%
  mutate(
    bin = cut(
      logsize_t0,
      breaks = 10)) %>%
  group_by(bin) %>%
  summarise(
    logsize_t0 = mean(
      logsize_t0,
      na.rm = TRUE),
    seed_nr = mean(
      seed_nr,
      na.rm = TRUE),
    se = sd(
      seed_nr,
      na.rm = TRUE) / sqrt(n()),
    .groups = "drop") %>%
  mutate(
    lwr = seed_nr - 1.96 * se,
    upr = seed_nr + 1.96 * se)

fig_se_n_line <- ggplot() +
  geom_jitter(
    data = df_se_n,
    aes(
      x = logsize_t0,
      y = seed_nr),
    alpha = 0.25,
    width = 0.08,
    height = 0.3) +
  geom_line(
    data = df_se_n_pred,
    aes(
      x = logsize_t0,
      y = seed_nr),
    linewidth = 0.9) +
  labs(
    x = expression("log(leaf area)"[t0]),
    y = "Number of seeds") +
  theme_bw()

fig_se_n_bin <- ggplot() +
  geom_point(
    data = df_se_n_binned,
    aes(
      x = logsize_t0,
      y = seed_nr)) +
  geom_errorbar(
    data = df_se_n_binned,
    aes(
      x = logsize_t0,
      ymin = lwr,
      ymax = upr),
    width = 0.2) +
  geom_line(
    data = df_se_n_pred,
    aes(
      x = logsize_t0,
      y = seed_nr),
    linewidth = 0.9) +
  labs(
    x = expression("log(leaf area)"[t0]),
    y = "Number of seeds") +
  theme_bw()

fig_se_n <- fig_se_n_line + fig_se_n_bin +
  plot_annotation(
    title = "Seed number conditional on seed production",
    subtitle = v_ggp_suffix)

fig_se_n


# Seed to recruit transition --------------------------------------------------
# Seeds produced in year t are related to recruits observed in year t + 1.

# Total observed seed production per population and year.
df_re_seed <- df %>%
  group_by(
    population, year) %>%
  summarise(
    nr_seeds = sum(
      seed_nr,
      na.rm = TRUE),
    nr_seed_obs = sum(
      !is.na(seed_nr)),
    .groups = "drop") %>%
  filter(
    nr_seed_obs > 0) %>%
  mutate(
    year_recruit = year + 1L)


# Number of recruits per sampled population and year.
df_re_sampled <- df %>%
  distinct(
    population, year)

df_re_count <- df %>%
  filter(
    recruit == 1) %>%
  count(
    population, year,
    name = "nr_recruits")

df_re_count <- df_re_sampled %>%
  left_join(
    df_re_count,
    by = c(
      "population",
      "year")) %>%
  mutate(
    nr_recruits = replace_na(
      nr_recruits,
      0L))


# Match seed production at t to recruits at t + 1.
df_re <- df_re_seed %>%
  inner_join(
    df_re_count,
    by = c(
      "population",
      "year_recruit" = "year"))


# Seed to recruit model -------------------------------------------------------
df_re_mod <- df_re %>%
  filter(nr_seeds > 0)

mod_re <- lm(
  nr_recruits ~ nr_seeds,
  data = df_re_mod)

summary(mod_re)

recr_per_seed <- unname(coef(mod_re)[["nr_seeds"]])

recr_per_seed


# Seed to recruit plot --------------------------------------------------------
df_re_newdata <- tibble(
  nr_seeds = seq(
    min(df_re_mod$nr_seeds),
    max(df_re_mod$nr_seeds), length.out = 200))

df_re_newdata <- bind_cols(
  df_re_newdata,
  as_tibble(predict(
    mod_re, newdata = df_re_newdata,
    interval = "confidence", level = 0.95)))

fig_re <- ggplot(
  df_re_mod, aes(x = nr_seeds, y = nr_recruits)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(
    data = df_re_newdata,
    aes(x = nr_seeds, ymin = lwr, ymax = upr),
    inherit.aes = FALSE, alpha = 0.2) +
  geom_line(
    data = df_re_newdata,
    aes(x = nr_seeds, y = fit),
    inherit.aes = FALSE, linewidth = 1) +
  labs(
    title = "Seed to recruit transition",
    subtitle = v_ggp_suffix,
    x = expression("Number of seeds"[t]),
    y = expression("Number of recruits"[t+1])) +
  theme_bw()

fig_re


# Recruit size distribution ---------------------------------------------------
df_re_size <- df %>%
  filter(
    recruit == 1,
    size_t0 > 0,
    is.finite(logsize_t0))

recr_sz <- mean(
  df_re_size$logsize_t0,
  na.rm = TRUE)

recr_sd <- sd(
  df_re_size$logsize_t0,
  na.rm = TRUE)


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

pars_gr <- bind_rows(
  pars_gr,
  pars_gr_var)


# Survival
pars_su <- tibble(
  coefficient = names(coef(mod_su_bestfit)),
  value = unname(coef(mod_su_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))


# Dormancy entry
pars_do <- tibble(
  coefficient = names(coef(mod_do_bestfit)),
  value = unname(coef(mod_do_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))


# Flowering
pars_fl <- tibble(
  coefficient = names(coef(mod_fl_bestfit)),
  value = unname(coef(mod_fl_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))


# Seed production probability
pars_se <- tibble(
  coefficient = names(coef(mod_se_bestfit)),
  value = unname(coef(mod_se_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))


# Seed number conditional on seed production
pars_se_n <- tibble(
  coefficient = names(coef(mod_se_n_bestfit)),
  value = unname(coef(mod_se_n_bestfit))) %>%
  mutate(
    coefficient = if_else(
      coefficient == '(Intercept)', 'b0', coefficient))


# Mesh limits
mesh_limits <- range(
  c(
    df_gr$logsize_t0,
    df_gr$logsize_t1,
    df_re_size$logsize_t0,
    df_ra_size$logsize_reactivate,
    df_do_entry_size$logsize_entry),
  na.rm = TRUE,
  finite = TRUE)


# Other parameters
pars_other <- tibble(
  coefficient = c(
    'recr_sz',
    'recr_sd',
    'react_b0',
    'react_b1',
    'react_sd',
    'p_dorm_survival',
    'p_reactivate_1',
    'recr_per_seed',
    'max_siz',
    'min_siz'),
  value = c(
    recr_sz,
    recr_sd,
    react_b0,
    react_b1,
    react_sd,
    p_dorm_survival,
    p_reactivate_1,
    recr_per_seed,
    mesh_limits[2],
    mesh_limits[1]))


# Save parameter estimates ----------------------------------------------------
write.csv(
  pars_gr,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_grow_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_su,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_surv_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_do,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_dorm_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_fl,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_flower_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_se,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_seed_prod_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_se_n,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_seed_num_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)

write.csv(
  pars_other,
  file.path(
    dir_data,
    paste0(
      v_script_prefix, '_', v_sp_abb,
      '_other_pars_mean', v_output_suffix, '.csv')),
  row.names = FALSE)


# Building the IPM from scratch -----------------------------------------------
extr_value <- function(x, field) {
  subset(x, coefficient == field)$value}

# Reference active abundance in the recruitment-fitting population-years.
df_re_n <- df %>%
  filter(
    !is.na(persistence_t0), persistence_t0 != 'DEAD', dormant_t0 == 0) %>%
  group_by(population, year) %>%
  summarise(n_active = n_distinct(id), .groups = 'drop') %>%
  semi_join(df_re_mod, by = c('population', 'year'))

pars <- Filter(
  function(x) length(x) > 0,
  list(
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
    
    dorm_b0 = extr_value(pars_do, 'b0'),
    dorm_b1 = extr_value(pars_do, 'logsize_t0'),
    dorm_b2 = extr_value(pars_do, 'logsize_t0_2'),
    dorm_b3 = extr_value(pars_do, 'logsize_t0_3'),
    
    flower_b0 = extr_value(pars_fl, 'b0'),
    flower_b1 = extr_value(pars_fl, 'logsize_t0'),
    flower_b2 = extr_value(pars_fl, 'logsize_t0_2'),
    flower_b3 = extr_value(pars_fl, 'logsize_t0_3'),
    
    seed_b0 = extr_value(pars_se, 'b0'),
    seed_b1 = extr_value(pars_se, 'logsize_t0'),
    seed_b2 = extr_value(pars_se, 'logsize_t0_2'),
    seed_b3 = extr_value(pars_se, 'logsize_t0_3'),
    
    seed_n_b0 = extr_value(pars_se_n, 'b0'),
    seed_n_b1 = extr_value(pars_se_n, 'logsize_t0'),
    seed_n_b2 = extr_value(pars_se_n, 'logsize_t0_2'),
    seed_n_b3 = extr_value(pars_se_n, 'logsize_t0_3'),
    
    recr_sz = extr_value(pars_other, 'recr_sz'),
    recr_sd = extr_value(pars_other, 'recr_sd'),
    recr_b0 = unname(coef(mod_re)[['(Intercept)']]),
    recr_per_seed = unname(coef(mod_re)[['nr_seeds']]),
    recr_n_ref = mean(df_re_n$n_active),
    
    react_b0 = extr_value(pars_other, 'react_b0'),
    react_b1 = extr_value(pars_other, 'react_b1'),
    react_sd = extr_value(pars_other, 'react_sd'),
    
    p_dorm_survival = extr_value(pars_other, 'p_dorm_survival'),
    p_reactivate_1 = extr_value(pars_other, 'p_reactivate_1'),
    
    L = extr_value(pars_other, 'min_siz'),
    U = extr_value(pars_other, 'max_siz'),
    mat_siz = 200,
    
    mod_gr_index = v_mod_gr_index,
    mod_su_index = v_mod_su_index,
    mod_do_index = v_mod_do_index,
    mod_fl_index = v_mod_fl_index,
    mod_se_index = v_mod_se_index,
    mod_se_n_index = v_mod_se_n_index))

write.csv(
  as.data.frame(pars),
  file.path(
    dir_data,
    paste0(v_script_prefix, '_', v_sp_abb,
           '_pars', v_output_suffix, '.csv')),
  row.names = FALSE)


# IPM functions ---------------------------------------------------------------
grow_sd <- function(x, pars) {
  sqrt(pars$a * exp(pars$b * x))}


# Growth from size x to size y.
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- pars$grow_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('grow_b', i)
      if (!is.null(pars[[param_name]])) {
        mean_value <- mean_value + pars[[param_name]] * x^i}}}
  
  sd_value <- grow_sd(x, pars)
  dnorm(y, mean = mean_value, sd = sd_value)}


inv_logit <- function(x) {
  plogis(x)}


# Survival of an x-sized active individual.
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('surv_b', i)
      if (!is.null(pars[[param_name]])) {
        survival_value <- survival_value + pars[[param_name]] * x^i}}}
  
  inv_logit(survival_value)}


# Dormancy entry conditional on survival.
dx <- function(x, pars, num_pars = v_mod_do_index) {
  dorm_value <- pars$dorm_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('dorm_b', i)
      if (!is.null(pars[[param_name]])) {
        dorm_value <- dorm_value + pars[[param_name]] * x^i}}}
  
  inv_logit(dorm_value)}


# Flowering probability.
flx <- function(x, pars, num_pars = v_mod_fl_index) {
  flower_value <- pars$flower_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('flower_b', i)
      if (!is.null(pars[[param_name]])) {
        flower_value <- flower_value + pars[[param_name]] * x^i}}}
  
  inv_logit(flower_value)}


# Positive seed production conditional on flowering.
sex <- function(x, pars, num_pars = v_mod_se_index) {
  seed_value <- pars$seed_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('seed_b', i)
      if (!is.null(pars[[param_name]])) {
        seed_value <- seed_value + pars[[param_name]] * x^i}}}
  
  inv_logit(seed_value)}


# Seed number conditional on positive seed production.
se_nx <- function(x, pars, num_pars = v_mod_se_n_index) {
  seed_n_value <- pars$seed_n_b0
  
  if (num_pars >= 1) {
    for (i in seq_len(num_pars)) {
      param_name <- paste0('seed_n_b', i)
      if (!is.null(pars[[param_name]])) {
        seed_n_value <- seed_n_value + pars[[param_name]] * x^i}}}
  
  exp(seed_n_value)}


# Expected seed production of an active individual.
seedx <- function(x, pars) {
  flx(x, pars) * sex(x, pars) * se_nx(x, pars)}


# Active survival-growth transition.
pxy <- function(x, y, pars) {
  sx(x, pars) * (1 - dx(x, pars)) * gxy(x, y, pars)}


# Kernel ----------------------------------------------------------------------
# Rows = destinations; columns = sources.
# 1:n = active; (n+1):(2*n) = dormant.
# Dormant plants retain size before dormancy.
# Non-return at the next census is treated as death.
kernel <- function(pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  i_dorm <- n + seq_len(n)
  
  # Growth distribution
  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Eviction correction
  for (i in seq_len(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])}
  
  for (i in ((n / 2) + 1):n) {
    Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])}
  
  # Active survival and dormancy entry
  Smat <- sx(y, pars)
  Dvec <- dx(y, pars)
  Tmat <- matrix(0, n, n)
  
  for (i in seq_len(n)) {
    Tmat[, i] <- Gmat[, i] * Smat[i] * (1 - Dvec[i])}
  
  # Reactivation size conditional on entry size
  react_sd_use <- max(h / 10, pars$react_sd, na.rm = TRUE)
  
  Rmat <- vapply(y, function(x) {
    react_y <- dnorm(
      y, mean = pars$react_b0 + pars$react_b1 * x,
      sd = react_sd_use, log = TRUE)
    react_y <- exp(react_y - max(react_y))
    react_y / sum(react_y)}, numeric(n))
  
  # Recruit size distribution
  recr_sd_use <- max(h / 10, pars$recr_sd, na.rm = TRUE)
  recr_y <- dnorm(y, mean = pars$recr_sz, sd = recr_sd_use) * h
  recr_y <- recr_y / sum(recr_y)
  
  Pmat <- matrix(0, 2 * n, 2 * n)
  Fmat <- matrix(0, 2 * n, 2 * n)
  
  # Active -> active
  Pmat[1:n, 1:n] <- Tmat
  
  # Active -> dormant: diagonal entry-size memory.
  Pmat[i_dorm, 1:n] <- diag(Smat * Dvec, nrow = n, ncol = n)
  
  # Dormant -> active: return probability and reactivation size.
  Pmat[1:n, i_dorm] <-
    pars$p_dorm_survival * pars$p_reactivate_1 * Rmat
  
  # Dormant -> dormant: no persistence; non-return is treated as death.
  Pmat[i_dorm, i_dorm] <- 0
  
  # Population intercept allocated per active plant at reference abundance.
  Fmat[1:n, 1:n] <- outer(
    recr_y,
    pars$recr_b0 / pars$recr_n_ref + pars$recr_per_seed * seedx(y, pars))
  
  k_yx <- Pmat + Fmat
  
  list(
    k_yx = k_yx, Pmat = Pmat, Fmat = Fmat,
    Tmat = Tmat, Gmat = Gmat, Rmat = Rmat,
    meshpts = y,
    dormant_index = i_dorm)}


lambda_ipm <- function(i) {
  Re(eigen(kernel(i)$k_yx)$values[1])}


# Mean IPM population growth rate ---------------------------------------------
lam_mean <- lambda_ipm(pars)
lam_mean


# Observed population growth rate ---------------------------------------------
# Count living individuals in each sampled population-year.
df_pop_n <- df %>%
  filter(!is.na(persistence_t0), persistence_t0 != 'DEAD') %>%
  group_by(population, year) %>%
  summarise(n = n_distinct(id), .groups = 'drop')

pop_counts_t0 <- df_pop_n %>%
  transmute(population, year = year + 1L, n_t0 = n)

pop_counts_t1 <- df_pop_n %>%
  transmute(population, year, n_t1 = n)

pop_counts <- inner_join(
  pop_counts_t0, pop_counts_t1, by = c('population', 'year')) %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0), n_t1 = sum(n_t1), .groups = 'drop') %>%
  mutate(obs_pgr = n_t1 / n_t0)

lam_mean_count <- exp(mean(log(pop_counts$obs_pgr), na.rm = TRUE))
lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)

lam_mean_count
lam_mean_overall


# Observed population growth rate ---------------------------------------------
# Restrict counts to the reliable demographic period.
last_reliable_t1 <- max(
  df$year_t1[df$reliable_demography], na.rm = TRUE)

df_pop_n <- df %>%
  filter(!is.na(persistence_t0), persistence_t0 != 'DEAD') %>%
  group_by(population, year) %>%
  summarise(n = n_distinct(id), .groups = 'drop')

pop_counts_t0 <- df_pop_n %>%
  transmute(population, year_t1 = year + 1L, n_t0 = n)

pop_counts_t1 <- df_pop_n %>%
  transmute(population, year_t1 = year, n_t1 = n)

pop_counts <- inner_join(
  pop_counts_t0, pop_counts_t1, by = c('population', 'year_t1')) %>%
  filter(year_t1 <= last_reliable_t1) %>%
  group_by(year_t1) %>%
  summarise(
    n_t0 = sum(n_t0), n_t1 = sum(n_t1), .groups = 'drop') %>%
  mutate(obs_pgr = n_t1 / n_t0)

lam_mean_count <- exp(mean(log(pop_counts$obs_pgr), na.rm = TRUE))
lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)

lam_mean_count
lam_mean_overall


# Building the IPM with ipmr ---------------------------------------------------
proto_ipm_p <- init_ipm(
  sim_gen = 'general', di_dd = 'di', det_stoch = 'det') %>%
  
  # Active -> active
  define_kernel(
    name = 'P',
    family = 'CC',
    formula = s * (1 - p_do) * g * d_size,
    
    s = plogis(
      surv_b0 +
        (if (mod_su_index >= 1) surv_b1 * size_1 else 0) +
        (if (mod_su_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3 * size_1^3 else 0)),
    
    p_do = plogis(
      dorm_b0 +
        (if (mod_do_index >= 1) dorm_b1 * size_1 else 0) +
        (if (mod_do_index >= 2) dorm_b2 * size_1^2 else 0) +
        (if (mod_do_index >= 3) dorm_b3 * size_1^3 else 0)),
    
    mu_g = grow_b0 +
      (if (mod_gr_index >= 1) grow_b1 * size_1 else 0) +
      (if (mod_gr_index >= 2) grow_b2 * size_1^2 else 0) +
      (if (mod_gr_index >= 3) grow_b3 * size_1^3 else 0),
    
    g = dnorm(size_2, mu_g, grow_sig),
    grow_sig = sqrt(a * exp(b * size_1)),
    
    data_list = pars,
    states = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun = 'norm', target = 'g')) %>%
  
  # Active -> dormant
  define_kernel(
    name = 'D_enter',
    family = 'CC',
    formula = s * p_do * entry_memory * d_size,
    
    entry_memory = as.numeric(size_1 == dorm_2) / d_size,
    
    s = plogis(
      surv_b0 +
        (if (mod_su_index >= 1) surv_b1 * size_1 else 0) +
        (if (mod_su_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3 * size_1^3 else 0)),
    
    p_do = plogis(
      dorm_b0 +
        (if (mod_do_index >= 1) dorm_b1 * size_1 else 0) +
        (if (mod_do_index >= 2) dorm_b2 * size_1^2 else 0) +
        (if (mod_do_index >= 3) dorm_b3 * size_1^3 else 0)),
    
    data_list = pars,
    states = list(c('size', 'dorm')),
    uses_par_sets = FALSE,
    evict_cor = FALSE) %>%
  
  # Dormant -> active
  define_kernel(
    name = 'D_leave',
    family = 'CC',
    formula = p_dorm_survival * p_reactivate_1 * r_d * d_dorm,
    
    react_mu = react_b0 + react_b1 * dorm_1,
    react_sig = max(d_size / 10, react_sd, na.rm = TRUE),
    r_d = dnorm(size_2, react_mu, react_sig),
    
    data_list = pars,
    states = list(c('dorm', 'size')),
    uses_par_sets = FALSE,
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun = 'norm', target = 'r_d')) %>%
  
  # Dormant -> dormant: the all-zero block.
  define_kernel(
    name = 'D_stay',
    family = 'CC',
    formula = 0 * dorm_1 * dorm_2 * d_dorm,
    
    data_list = pars,
    states = list(c('dorm')),
    uses_par_sets = FALSE,
    evict_cor = FALSE) %>%
  
  # Active -> recruits: intercept and seed production both contribute.
  define_kernel(
    name = 'F',
    family = 'CC',
    formula =
      (recr_b0 / recr_n_ref + p_fl * p_se * n_se * recr_per_seed) *
      r_d * d_size,
    
    p_fl = plogis(
      flower_b0 +
        (if (mod_fl_index >= 1) flower_b1 * size_1 else 0) +
        (if (mod_fl_index >= 2) flower_b2 * size_1^2 else 0) +
        (if (mod_fl_index >= 3) flower_b3 * size_1^3 else 0)),
    
    p_se = plogis(
      seed_b0 +
        (if (mod_se_index >= 1) seed_b1 * size_1 else 0) +
        (if (mod_se_index >= 2) seed_b2 * size_1^2 else 0) +
        (if (mod_se_index >= 3) seed_b3 * size_1^3 else 0)),
    
    n_se = exp(
      seed_n_b0 +
        (if (mod_se_n_index >= 1) seed_n_b1 * size_1 else 0) +
        (if (mod_se_n_index >= 2) seed_n_b2 * size_1^2 else 0) +
        (if (mod_se_n_index >= 3) seed_n_b3 * size_1^3 else 0)),
    
    recr_sig = max(d_size / 10, recr_sd, na.rm = TRUE),
    r_d = dnorm(size_2, recr_sz, recr_sig),
    
    data_list = pars,
    states = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun = 'norm', target = 'r_d')) %>%
  
  # Kernel implementation
  define_impl(
    make_impl_args_list(
      kernel_names = c('P', 'D_enter', 'D_leave', 'D_stay', 'F'),
      int_rule = rep('midpoint', 5),
      state_start = c('size', 'size', 'dorm', 'dorm', 'size'),
      state_end = c('size', 'dorm', 'size', 'dorm', 'size'))) %>%
  
  # Domains
  define_domains(
    size = c(pars$L, pars$U, pars$mat_siz),
    dorm = c(pars$L, pars$U, pars$mat_siz)) %>%
  
  # Initial population
  define_pop_state(
    pop_vectors = list(
      n_size = rep(1 / pars$mat_siz, pars$mat_siz),
      n_dorm = rep(0, pars$mat_siz)))


ipmr_p <- make_ipm(proto_ipm = proto_ipm_p, iterations = 200)


# Mean IPM population growth rate ---------------------------------------------
lam_mean_ipmr <- lambda(ipmr_p)
lam_mean_ipmr


# Compare hand-built and ipmr implementations ----------------------------------
lam_mean
lam_mean_ipmr


# Export lambda ---------------------------------------------------------------
lam_out <- data.frame(
  coefficient = names(lam_mean_ipmr), value = lam_mean_ipmr)

lam_out_wide <- as.list(
  pivot_wider(lam_out, names_from = 'coefficient', values_from = 'value'))