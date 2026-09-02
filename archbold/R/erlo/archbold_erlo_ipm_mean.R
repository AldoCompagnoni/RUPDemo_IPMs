# IPM mean - Archbold - Eriogonum longifolium

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.08.19

# Study organism: Eriogonum longifolium var. gnaphalifolium
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.226.1
# Meta data link:
# https://portal.edirepository.org/nis/metadataviewer?packageid=edi.226.1
# Citing publication: Satterthwaite et al. 2002, Ecological Applications
# Time period: 1990-2013


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(MASS, tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, GGally)


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Eriogonum longifolium')
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c()

# Create a unique species abbreviation for file naming
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Define script prefix
v_script_prefix <- str_c(v_head)

# Plot subtitle
v_ggp_suffix    <- paste(
  tools::toTitleCase(v_head), '-', v_species)

# Models
v_mod_set_su   <- c(2)
v_mod_set_gr   <- c()
v_mod_set_do   <- c()
v_mod_set_fl   <- c()
v_mod_set_fl_n <- c()


# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_head))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

if (!dir.exists(dir_R     )) {dir.create(dir_R     )}
if (!dir.exists(dir_data  )) {dir.create(dir_data  )}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Functions --------------------------------------------------------------------
# function to plot your survival data 'binned' (instead of 'jittered')
source('helper_functions/plot_binned_prop.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0("ab_", v_sp_abb, "_df_workdata_260820.csv"))) %>%
  mutate(
    year = as.numeric(year),
    row_type = as.character(row_type))

df_ind <- df %>%
  filter(row_type == "individual")

df_re <- df %>%
  filter(
    row_type == "recruitment",
    recruitment_complete,
    !is.na(recruits_simple))


# Survival ---------------------------------------------------------------------
df_su <- df %>%
  filter(
    state == "active",
    !is.na(survives),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  dplyr::select(
    id, year, size_t0, survives,
    logsize_t0, logsize_t0_2, logsize_t0_3)

fig_su_raw <- ggplot(
  data = plot_binned_prop(df_su, 10, logsize_t0, survives)) +
  geom_jitter(
    data = df_su,
    aes(x = logsize_t0, y = survives),
    position = position_jitter(width = 0.1, height = 0.3),
    alpha = 0.1) +
  geom_point(
    aes(x = logsize_t0, y = survives),
    pch = 16, color = "red") +
  geom_errorbar(
    aes(x = logsize_t0, ymin = lwr, ymax = upr),
    linewidth = 0.5, width = 0.1) +
  scale_y_continuous(limits = c(0, 1.01)) +
  theme_bw() +
  labs(
    title = "Survival",
    subtitle = v_ggp_suffix,
    x = expression("log(diameter)"[t0]),
    y = expression("Survival to time t1"))

fig_su_raw


# Survival model ---------------------------------------------------------------
mod_su_0 <- glm(
  survives ~ 1,
  data = df_su, family = "binomial")

mod_su_1 <- glm(
  survives ~ logsize_t0,
  data = df_su, family = "binomial")

mod_su_2 <- glm(
  survives ~ logsize_t0 + logsize_t0_2,
  data = df_su, family = "binomial")

mod_su_3 <- glm(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_su, family = "binomial")

mods_su <- list(
  mod_su_0, mod_su_1, mod_su_2, mod_su_3)

mods_su_dAICc <- AICctab(
  mods_su, weights = TRUE, sort = FALSE)$dAICc

mods_su_sorted <- order(mods_su_dAICc)

if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index <- mod_su_index_bestfit - 1
} else {
  mod_su_index_bestfit <- v_mod_set_su + 1
  v_mod_su_index <- v_mod_set_su
}

mod_su_bestfit <- mods_su[[mod_su_index_bestfit]]

mod_su_bestfit
mods_su_dAICc


# Survival plot ----------------------------------------------------------------
df_su_newdata <- df_su %>%
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
      mod_su_bestfit,
      newdata = .,
      type = "response"))

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
    title = "Survival probability by size",
    subtitle = v_ggp_suffix,
    x = expression("log(diameter)"[t0]),
    y = "Probability of survival") +
  theme_bw()


df_su_bindata <- plot_binned_prop(
  df_su, 10, logsize_t0, survives)

df_su_pred <- df_su %>%
  reframe(
    logsize_t0 = seq(
      min(logsize_t0),
      max(logsize_t0),
      length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  mutate(
    survives = predict(
      mod_su_bestfit,
      newdata = .,
      type = "response"))

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
    data = df_su_pred,
    aes(x = logsize_t0, y = survives),
    linewidth = 1.2) +
  labs(
    x = expression("log(diameter)"[t0]),
    y = "") +
  theme_bw() +
  ylim(0, 1)

fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(
    state == "active",
    state_t1 == "active",
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1)) %>%
  dplyr::select(
    id, year, size_t0, size_t1,
    logsize_t0, logsize_t1,
    logsize_t0_2, logsize_t0_3)


# Growth model -----------------------------------------------------------------
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

mods_gr_dAICc <- AICctab(
  mods_gr, weights = TRUE, sort = FALSE)$dAICc

mods_gr_sorted <- order(mods_gr_dAICc)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index <- mod_gr_index_bestfit - 1
} else {
  mod_gr_index_bestfit <- v_mod_set_gr + 1
  v_mod_gr_index <- v_mod_set_gr
}

mod_gr_bestfit <- mods_gr[[mod_gr_index_bestfit]]

mod_gr_bestfit
mods_gr_dAICc


# Growth plot ------------------------------------------------------------------
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
      interval = "confidence")))

fig_gr <- ggplot(
  df_gr,
  aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.4) +
  geom_line(
    data = df_gr_newdata,
    aes(y = fit),
    linewidth = 1) +
  geom_ribbon(
    data = df_gr_newdata,
    aes(y = fit, ymin = lwr, ymax = upr),
    alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    title = "Growth prediction",
    subtitle = v_ggp_suffix,
    x = expression("log(diameter)"[t0]),
    y = expression("log(diameter)"[t1])) +
  theme_bw()

fig_gr


# Growth variance --------------------------------------------------------------
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)  
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2  
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE)) 



# Dormancy entry data ----------------------------------------------------------
df_do <- df %>%
  filter(
    !is.na(enter_dormancy),
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  dplyr::select(
    id, year, size_t0, enter_dormancy,
    logsize_t0, logsize_t0_2, logsize_t0_3)


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


# Reactivation -----------------------------------------------------------------
df_ra <- df %>%
  filter(!is.na(reactivate))

mod_ra <- glm(
  reactivate ~ 1,
  data = df_ra, family = "binomial")

p_reactivate <- predict(
  mod_ra,
  newdata = data.frame(x = 1),
  type = "response")[1]

p_reactivate

df_ra_size <- df %>%
  filter(
    reactivate == 1,
    size_reactivate_t1 > 0) %>%
  mutate(
    logsize_reactivate =
      log(size_reactivate_t1))

react_sz <- mean(
  df_ra_size$logsize_reactivate,
  na.rm = TRUE)

react_sd <- sd(
  df_ra_size$logsize_reactivate,
  na.rm = TRUE)


# Reactivation plot ------------------------------------------------------------
df_ra_plot <- df_ra %>%
  summarise(
    n = n(),
    nr_reactivate = sum(reactivate == 1),
    reactivation = mean(reactivate)) %>%
  mutate(
    lwr = binom.confint(
      nr_reactivate, n, methods = "wilson")$lower,
    upr = binom.confint(
      nr_reactivate, n, methods = "wilson")$upper)

fig_ra <- ggplot(
  df_ra_plot,
  aes(x = 1, y = reactivation)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    width = 0.1) +
  geom_hline(
    yintercept = p_reactivate,
    linetype = "dashed") +
  scale_x_continuous(
    breaks = 1,
    labels = "Dormant plants") +
  scale_y_continuous(
    limits = c(0, 1)) +
  labs(
    title = "Reactivation probability",
    subtitle = v_ggp_suffix,
    x = "",
    y = "Probability of reactivation") +
  theme_bw()

fig_ra


# Flower data ------------------------------------------------------------------
df_fl <- df %>%
  filter(
    state == "active",
    !is.na(flower),
    size_t0 > 0,
    is.finite(logsize_t0))


# Flower model -----------------------------------------------------------------
mod_fl_0 <- glm(
  flower ~ 1,
  data = df_fl, family = "binomial")

mod_fl_1 <- glm(
  flower ~ logsize_t0,
  data = df_fl, family = "binomial")

mod_fl_2 <- glm(
  flower ~ logsize_t0 + logsize_t0_2,
  data = df_fl, family = "binomial")

mod_fl_3 <- glm(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fl, family = "binomial")

mods_fl <- list(
  mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)

mods_fl_dAICc <- AICctab(
  mods_fl, weights = TRUE, sort = FALSE)$dAICc

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


# Flowering probability plot ---------------------------------------------------
df_fl_newdata <- df_fl %>%
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
      mod_fl_bestfit,
      newdata = .,
      type = "response"))

fig_fl_line <- ggplot(
  df_fl,
  aes(x = logsize_t0, y = flower)) +
  geom_jitter(
    height = 0.05, width = 0,
    alpha = 0.3) +
  geom_line(
    data = df_fl_newdata,
    aes(y = predicted),
    linewidth = 1) +
  labs(
    title = "Flowering probability by size",
    subtitle = v_ggp_suffix,
    x = expression("log(diameter)"[t0]),
    y = "Probability of flowering") +
  theme_bw()


df_fl_bindata <- plot_binned_prop(
  df_fl, 10, logsize_t0, flower)

df_fl_pred <- df_fl %>%
  reframe(
    logsize_t0 = seq(
      min(logsize_t0),
      max(logsize_t0),
      length.out = 100)) %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  mutate(
    flower = predict(
      mod_fl_bestfit,
      newdata = .,
      type = "response"))

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
    data = df_fl_pred,
    aes(x = logsize_t0, y = flower),
    linewidth = 1.2) +
  labs(
    x = expression("log(diameter)"[t0]),
    y = "") +
  theme_bw() +
  ylim(0, 1)

fig_fl <- fig_fl_line + fig_fl_bin + plot_layout()
fig_fl


# Number of scapes conditional on flowering data -------------------------------
df_fl_cond <- df_fl %>%
  filter(
    flower == 1,
    !is.na(fl_nr),
    fl_nr > 0,
    fl_nr %% 1 == 0)


# Number of scapes models ------------------------------------------------------
mod_fl_n_0 <- glm.nb(
  fl_nr ~ 1,
  data = df_fl_cond)

mod_fl_n_1 <- glm.nb(
  fl_nr ~ logsize_t0,
  data = df_fl_cond)

mod_fl_n_2 <- glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2,
  data = df_fl_cond)

mod_fl_n_3 <- glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
  data = df_fl_cond)

mods_fl_n <- list(
  mod_fl_n_0, mod_fl_n_1, mod_fl_n_2, mod_fl_n_3)

mods_fl_n_dAICc <- AICctab(
  mods_fl_n, weights = TRUE, sort = FALSE)$dAICc

mods_fl_n_sorted <- order(mods_fl_n_dAICc)

if (length(v_mod_set_fl_n) == 0) {
  mod_fl_n_index_bestfit <- mods_fl_n_sorted[1]
  v_mod_fl_n_index <- mod_fl_n_index_bestfit - 1
} else {
  mod_fl_n_index_bestfit <- v_mod_set_fl_n + 1
  v_mod_fl_n_index <- v_mod_set_fl_n
}

mod_fl_n_bestfit <- mods_fl_n[[mod_fl_n_index_bestfit]]

mod_fl_n_bestfit
mods_fl_n_dAICc


# Predictions for flower number -----------------------------------------------
# Create prediction grid
df_fl_n_pred <- expand.grid(
  logsize_t0 = seq(
    min(df_fl_cond$logsize_t0),
    max(df_fl_cond$logsize_t0),
    length.out = 100))

df_fl_n_pred <- df_fl_n_pred %>%
  mutate(
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

# Predict
df_fl_n_pred$fl_nr <- predict(
  mod_fl_n_bestfit,
  newdata = df_fl_n_pred,
  type = "response")


# Binned observed data
df_fl_n_binned <- df_fl_cond %>%
  mutate(bin = cut(logsize_t0, breaks = 10)) %>%
  group_by(bin) %>%
  summarise(
    logsize_t0 = mean(logsize_t0, na.rm = TRUE),
    fl_nr = mean(fl_nr, na.rm = TRUE),
    se = sd(fl_nr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop") %>%
  mutate(
    lwr = fl_nr - 1.96 * se,
    upr = fl_nr + 1.96 * se)


# Flower number plots
# Plot 1: Raw jitter + prediction
fig_fl_n_line <- ggplot() +
  geom_jitter(
    data = df_fl_cond,
    aes(x = logsize_t0, y = fl_nr),
    alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(
    data = df_fl_n_pred,
    aes(x = logsize_t0, y = fl_nr),
    linewidth = 0.9) +
  theme_bw() +
  labs(
    title = NULL,
    x = expression("log(diameter)"[t0]),
    y = "Number of flowering scapes")


# Plot 2: Binned + prediction
fig_fl_n_bin <- ggplot() +
  geom_point(
    data = df_fl_n_binned,
    aes(x = logsize_t0, y = fl_nr)) +
  geom_errorbar(
    data = df_fl_n_binned,
    aes(x = logsize_t0, ymin = lwr, ymax = upr),
    width = 0.2) +
  geom_line(
    data = df_fl_n_pred,
    aes(x = logsize_t0, y = fl_nr),
    linewidth = 0.9) +
  theme_bw() +
  labs(
    title = NULL,
    x = expression("log(diameter)"[t0]),
    y = "Number of flowering scapes")


# Combine
fig_fl_n <- fig_fl_n_line + fig_fl_n_bin +
  plot_annotation(
    title = "Flower number",
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic")))

fig_fl_n


# Recruitment data -------------------------------------------------------------
df_sc2re <- df_re %>%
  filter(
    !is.na(nr_scapes),
    !is.na(recruits_simple)) %>%
  group_by(site, pop, year) %>%
  summarise(
    total_scapes = sum(nr_scapes),
    recruit_count = sum(recruits_simple),
    .groups = "drop")

df_sc2re_mod <- df_sc2re %>%
  filter(
    !is.na(total_scapes),
    !is.na(recruit_count))


# Recruitment model ------------------------------------------------------------
# Recruitment can occur without current flowering, consistent with recruitment
# from previously accumulated seed.

mod_re_0 <- glm.nb(
  recruit_count ~ 1,
  data = df_sc2re_mod)

mod_re_sc <- glm.nb(
  recruit_count ~ total_scapes,
  data = df_sc2re_mod)

mod_re_logsc <- glm.nb(
  recruit_count ~ log1p(total_scapes),
  data = df_sc2re_mod)

mods_re <- list(
  mod_re_0,
  mod_re_sc,
  mod_re_logsc)

mods_re_dAICc <- AICctab(
  mods_re,
  weights = TRUE,
  sort = FALSE)$dAICc

mods_re_sorted <- order(mods_re_dAICc)
mod_re_bestfit <- mods_re[[mods_re_sorted[1]]]

mods_re_dAICc
summary(mod_re_bestfit)


# Mean recruitment conversion for the IPM -------------------------------------
# The recruitment model includes baseline recruitment, which can represent
# recruitment from the seedbank. For the linear mean IPM, collapse the fitted
# recruitment relationship to one mean recruits-per-scape conversion.

df_sc2re_mod <- df_sc2re_mod %>%
  mutate(
    recruits_pred = predict(
      mod_re_bestfit,
      type = "response"))

recr_per_scape <- sum(df_sc2re_mod$recruits_pred) /
  sum(df_sc2re_mod$total_scapes)

recr_per_scape


# Flowering scapes to recruits plot --------------------------------------------
df_sc2re_pred <- data.frame(
  total_scapes = seq(
    min(df_sc2re_mod$total_scapes),
    max(df_sc2re_mod$total_scapes),
    length.out = 100))

mod_sc2re_preds <- predict(
  mod_re_bestfit,
  newdata = df_sc2re_pred,
  type = "link",
  se.fit = TRUE)

df_sc2re_pred <- df_sc2re_pred %>%
  mutate(
    fit = exp(mod_sc2re_preds$fit),
    lower = exp(
      mod_sc2re_preds$fit -
        1.96 * mod_sc2re_preds$se.fit),
    upper = exp(
      mod_sc2re_preds$fit +
        1.96 * mod_sc2re_preds$se.fit))

fig_sc2re <- ggplot(
  df_sc2re_mod,
  aes(x = total_scapes, y = recruit_count)) +
  geom_jitter(
    height = 0.2, width = 0.5,
    alpha = 0.4) +
  geom_ribbon(
    data = df_sc2re_pred,
    aes(
      x = total_scapes,
      ymin = lower,
      ymax = upper),
    inherit.aes = FALSE,
    alpha = 0.2) +
  geom_line(
    data = df_sc2re_pred,
    aes(x = total_scapes, y = fit),
    inherit.aes = FALSE,
    linewidth = 1.2) +
  labs(
    title = "Recruits t1 by flowering scapes t0",
    subtitle = v_ggp_suffix,
    x = expression("Total flowering scapes "[t0]),
    y = expression("Number of recruits "[t1])) +
  theme_bw()

fig_sc2re


# Recruit size -----------------------------------------------------------------
df_re_size <- df_ind %>%
  filter(
    recruit_type == "seedling",
    size_t0 > 0,
    is.finite(logsize_t0))

recr_sz <- mean(
  df_re_size$logsize_t0,
  na.rm = TRUE)

recr_sd <- sd(
  df_re_size$logsize_t0,
  na.rm = TRUE)



# Extracting parameter estimates -----------------------------------------------

# Survival
coef_su <- data.frame(
  coefficient = names(coef(mod_su_bestfit)),
  value = coef(mod_su_bestfit)) %>%
  mutate(
    coefficient = as.character(coefficient),
    coefficient = replace(
      coefficient, grepl("Intercept", coefficient), "b0"))

# Growth
coef_gr_fe <- data.frame(
  coefficient = names(coef(mod_gr_bestfit)),
  value = coef(mod_gr_bestfit))

coef_gr_var <- data.frame(
  coefficient = names(coef(mod_gr_var)),
  value = coef(mod_gr_var))

coef_gr <- bind_rows(coef_gr_fe, coef_gr_var) %>%
  mutate(
    coefficient = as.character(coefficient),
    coefficient = replace(
      coefficient, grepl("Intercept", coefficient), "b0"))

# Dormancy entry
coef_do <- data.frame(
  coefficient = names(coef(mod_do_bestfit)),
  value = coef(mod_do_bestfit)) %>%
  mutate(
    coefficient = as.character(coefficient),
    coefficient = replace(
      coefficient, grepl("Intercept", coefficient), "b0"))

# Flowering probability
coef_fl <- data.frame(
  coefficient = names(coef(mod_fl_bestfit)),
  value = coef(mod_fl_bestfit)) %>%
  mutate(
    coefficient = as.character(coefficient),
    coefficient = replace(
      coefficient, grepl("Intercept", coefficient), "b0"))

# Scape number conditional on flowering
coef_fln <- data.frame(
  coefficient = names(coef(mod_fl_n_bestfit)),
  value = coef(mod_fl_n_bestfit)) %>%
  mutate(
    coefficient = as.character(coefficient),
    coefficient = replace(
      coefficient, grepl("Intercept", coefficient), "b0"))

extr_value <- function(x, field) {
  subset(x, coefficient == field)$value
}

# IPM constants ----------------------------------------------------------------
coef_misc <- data.frame(
  coefficient = c(
    "recr_sz",
    "recr_sd",
    "react_sz",
    "react_sd",
    "p_reactivate",
    "fecu_b0",
    "max_siz",
    "min_siz"),
  value = c(
    recr_sz,
    recr_sd,
    react_sz,
    react_sd,
    p_reactivate,
    recr_per_scape,
    max(df_gr$logsize_t0),
    min(df_gr$logsize_t0)))

# Parameter list ---------------------------------------------------------------
pars <- Filter(function(x) length(x) > 0, list(
  prefix = v_script_prefix,
  species = v_species,
  
  surv_b0 = extr_value(coef_su, "b0"),
  surv_b1 = extr_value(coef_su, "logsize_t0"),
  surv_b2 = extr_value(coef_su, "logsize_t0_2"),
  surv_b3 = extr_value(coef_su, "logsize_t0_3"),
  
  grow_b0 = extr_value(coef_gr, "b0"),
  grow_b1 = extr_value(coef_gr, "logsize_t0"),
  grow_b2 = extr_value(coef_gr, "logsize_t0_2"),
  grow_b3 = extr_value(coef_gr, "logsize_t0_3"),
  a = extr_value(coef_gr, "a"),
  b = extr_value(coef_gr, "b"),
  
  dorm_b0 = extr_value(coef_do, "b0"),
  dorm_b1 = extr_value(coef_do, "logsize_t0"),
  dorm_b2 = extr_value(coef_do, "logsize_t0_2"),
  dorm_b3 = extr_value(coef_do, "logsize_t0_3"),
  
  fl_b0 = extr_value(coef_fl, "b0"),
  fl_b1 = extr_value(coef_fl, "logsize_t0"),
  fl_b2 = extr_value(coef_fl, "logsize_t0_2"),
  fl_b3 = extr_value(coef_fl, "logsize_t0_3"),
  
  fln_b0 = extr_value(coef_fln, "b0"),
  fln_b1 = extr_value(coef_fln, "logsize_t0"),
  fln_b2 = extr_value(coef_fln, "logsize_t0_2"),
  fln_b3 = extr_value(coef_fln, "logsize_t0_3"),
  
  fecu_b0 = extr_value(coef_misc, "fecu_b0"),
  recr_sz = extr_value(coef_misc, "recr_sz"),
  recr_sd = extr_value(coef_misc, "recr_sd"),
  react_sz = extr_value(coef_misc, "react_sz"),
  react_sd = extr_value(coef_misc, "react_sd"),
  p_reactivate = extr_value(coef_misc, "p_reactivate"),
  
  L = extr_value(coef_misc, "min_siz"),
  U = extr_value(coef_misc, "max_siz"),
  mat_siz = 200,
  
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index,
  mod_do_index = v_mod_do_index,
  mod_fl_index = v_mod_fl_index,
  mod_fl_n_index = v_mod_fl_n_index))


# Building the IPM -------------------------------------------------------------
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}


# Survival ---------------------------------------------------------------------
sx <- function(x, pars, num_pars = pars$mod_su_index) {
  val <- pars$surv_b0
  
  for (i in seq_len(num_pars)) {
    param <- paste0("surv_b", i)
    
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  
  inv_logit(val)
}


# Growth -----------------------------------------------------------------------
grow_mu <- function(x, pars, num_pars = pars$mod_gr_index) {
  val <- pars$grow_b0
  
  for (i in seq_len(num_pars)) {
    param <- paste0("grow_b", i)
    
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  
  val
}


# Growth variation -------------------------------------------------------------
grow_sd <- function(mu, pars) {
  sqrt(pars$a * exp(pars$b * mu))
}


# Growth transition ------------------------------------------------------------
gxy <- function(x, y, pars) {
  mu <- grow_mu(x, pars)
  
  dnorm(
    y,
    mean = mu,
    sd = grow_sd(mu, pars))
}

# Dormancy entry ---------------------------------------------------------------
dorm_x <- function(x, pars, num_pars = pars$mod_do_index) {
  val <- pars$dorm_b0
  
  for (i in seq_len(num_pars)) {
    param <- paste0("dorm_b", i)
    
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  
  inv_logit(val)
}


# Flowering probability --------------------------------------------------------
fl_x <- function(x, pars, num_pars = pars$mod_fl_index) {
  val <- pars$fl_b0
  
  for (i in seq_len(num_pars)) {
    param <- paste0("fl_b", i)
    
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  
  inv_logit(val)
}


# Scape number conditional on flowering ---------------------------------------
fl_n_x <- function(x, pars, num_pars = pars$mod_fl_n_index) {
  val <- pars$fln_b0
  
  for (i in seq_len(num_pars)) {
    param <- paste0("fln_b", i)
    
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  
  exp(val)
}

# Recruitment size distribution ------------------------------------------------
re_y_dist <- function(y, pars) {
  norm <- pnorm(
    pars$U, mean = pars$recr_sz, sd = pars$recr_sd) -
    pnorm(
      pars$L, mean = pars$recr_sz, sd = pars$recr_sd)
  
  dnorm(
    y,
    mean = pars$recr_sz,
    sd = pars$recr_sd) / norm
}


# Reactivation size distribution -----------------------------------------------
react_y_dist <- function(y, pars) {
  norm <- pnorm(
    pars$U, mean = pars$react_sz, sd = pars$react_sd) -
    pnorm(
      pars$L, mean = pars$react_sz, sd = pars$react_sd)
  
  dnorm(
    y,
    mean = pars$react_sz,
    sd = pars$react_sd) / norm
}

# F-kernel ---------------------------------------------------------------------
fyx <- function(y, x, pars) {
  fl_x(x, pars) *
    fl_n_x(x, pars) *
    pars$fecu_b0 *
    re_y_dist(y, pars)
}


# Kernel -----------------------------------------------------------------------
kernel <- function(pars) {
  
  # Continuous active-size classes
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  # Active survival and dormancy
  Smat <- sx(y, pars)
  Dmat <- dorm_x(y, pars)
  
  # Growth matrix
  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
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
    "*")
  
  # Active -> dormant
  Dorm_row <- matrix(
    Smat * Dmat,
    nrow = 1)
  
  # Dormant -> active
  React_vec <- react_y_dist(y, pars) * h
  React_vec <- React_vec / sum(React_vec)
  
  React_col <- matrix(
    pars$p_reactivate * React_vec,
    ncol = 1)
  
  # Dormant -> dormant
  Dorm_stasis <- 1 - pars$p_reactivate
  
  # Fertility from active individuals
  Fmat <- outer(
    y, y,
    Vectorize(function(y, x) {
      fyx(y, x, pars)
    })) * h
  
  # Full 201 x 201 kernel
  K <- rbind(
    cbind(Tmat + Fmat, React_col),
    cbind(Dorm_row, Dorm_stasis))
  
  return(list(
    k_yx = K,
    Fmat = Fmat,
    Tmat = Tmat,
    Gmat = Gmat,
    Dorm_row = Dorm_row,
    React_col = React_col,
    Smat = Smat,
    Dmat = Dmat,
    meshpts = y,
    h = h))
}

# Mean population growth rate --------------------------------------------------
lambda_ipm <- function(pars) {
  Re(eigen(
    kernel(pars)$k_yx,
    only.values = TRUE)$values[1])
}

lam_mean <- lambda_ipm(pars)
lam_mean


# Kernel checks ----------------------------------------------------------------
K_mean <- kernel(pars)

dim(K_mean$k_yx)

range(K_mean$Smat)
range(K_mean$Dmat)

range(colSums(K_mean$Gmat))
sum(K_mean$React_col)
1 - pars$p_reactivate

range(K_mean$k_yx)
any(!is.finite(K_mean$k_yx))

lam_mean


# Stable population structure -------------------------------------------------
eig_mean <- eigen(K_mean$k_yx)

stable_dist <- Re(eig_mean$vectors[, 1])

if (sum(stable_dist) < 0) {
  stable_dist <- -stable_dist
}

stable_dist <- stable_dist / sum(stable_dist)

stable_active <- stable_dist[seq_len(pars$mat_siz)]
stable_dormant <- stable_dist[pars$mat_siz + 1]

df_stable <- tibble(
  logsize = K_mean$meshpts,
  size = exp(logsize),
  stable = stable_active,
  stable_density = stable_active / K_mean$h)

df_stable_summary <- tibble(
  state = c("Active", "Dormant"),
  proportion = c(
    sum(stable_active),
    stable_dormant))

df_stable_summary


# Stable size distribution -----------------------------------------------------
fig_stable <- ggplot(
  df_stable,
  aes(x = size, y = stable_density)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Stable active size distribution",
    subtitle = v_ggp_suffix,
    x = "Rosette diameter (cm)",
    y = "Stable density")

fig_stable


# Observed population growth ---------------------------------------------------
# Count known living individuals. New adults are additionally counted one year
# earlier, matching the simplified recruitment assumption used in this IPM.

df_counts_alive <- df_ind %>%
  filter(state %in% c("active", "dormant", "alive_or_dormant")) %>%
  count(site, pop, qu, year, name = "n_alive")

df_newadult_backfill <- df_ind %>%
  filter(recruit_type == "new_adult") %>%
  transmute(site, pop, qu, year = year - 1) %>%
  count(site, pop, qu, year, name = "n_newadult_backfill")

# Use recruitment rows as the monitored quadrat-year inventory
df_counts_quad <- df %>%
  filter(row_type == "recruitment") %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_counts_alive,
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_newadult_backfill,
    by = c("site", "pop", "qu", "year")) %>%
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

# The final three transitions cannot fully resolve terminal dormancy
dormancy_cutoff <- 3
last_complete_transition <- max(df_ind$year, na.rm = TRUE) -
  dormancy_cutoff - 1

df_counts_year <- df_counts_quad %>%
  filter(year <= last_complete_transition) %>%
  group_by(year) %>%
  summarise(
    n_quads = n(),
    n_t0 = sum(n),
    n_t1 = sum(n_t1), .groups = "drop") %>%
  mutate(lambda_obs = n_t1 / n_t0)

df_counts_year

# Mean observed population growth
lam_obs_y <- df_counts_year$lambda_obs

lam_obs_mean <- mean(lam_obs_y, na.rm = TRUE)
lam_obs_geo <- exp(mean(log(lam_obs_y), na.rm = TRUE))

lam_obs_mean
lam_obs_geo

# Modeled versus observed mean population growth
c(
  lambda_ipm = lam_mean,
  lambda_obs_arithmetic = lam_obs_mean,
  lambda_obs_geometric = lam_obs_geo)


# Modeled versus observed population growth -----------------------------------
df_lambda_mean <- tibble(
  estimate = c(
    "IPM asymptotic",
    "Observed arithmetic",
    "Observed geometric"),
  lambda = c(
    lam_mean,
    lam_obs_mean,
    lam_obs_geo))

df_lambda_mean


# # Save data --------------------------------------------------------------------
# write.csv(df_og, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original.csv')))
# write.csv(df_meta, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_meta.csv')))
# write.csv(df, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata.csv')))
# write.csv(df_su, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_survival.csv')))
# write.csv(df_gr, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_growth.csv')))
# write.csv(df_fl, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_flower.csv')))
# write.csv(df_fr, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_fruit.csv')))
# write.csv(df_re, row.names = F,
#           file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_recruit.csv')))
# 
