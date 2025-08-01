# IPM mean - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.29

# Study organism: Polygala lewtonii
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.318.1
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.318.1
# Citing publication: https://doi.org/10.1071/BT11271
# Time periode: 2001-2017


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
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, MASS, GGally)


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Polygala lewtonii')
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
v_mod_set_su <- c()
v_mod_set_gr <- c()
v_mod_set_fl <- c()
v_mod_set_fe <- c()

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
df_og   <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original.csv')))
df_meta <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_meta.csv')))
df      <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata.csv')))


# Double ID --------------------------------------------------------------------

df %>%
  mutate(
    id_part1 = str_extract(id, "^([^_]+_[^_]+_[^_]+)"),
    id_part2 = sub("^[^_]+_[^_]+_[^_]+_", "", id)
  ) %>%
  group_by(site, quad, year, id_part1) %>%
  filter(n() > 1) %>%  # Keep only groups with more than one row
  ungroup()

# Remove all of them 
# Send the id and other info on them to Aldo


# Variability in the number of individuals per quadrats ------------------------
df %>%
  group_by(site, quad) %>%
  summarise(n_individuals = n_distinct(id), .groups = "drop") %>%
  summarise(
    total_quadrats   = n(),
    min_indivs       = min(n_individuals),
    max_indivs       = max(n_individuals),
    mean_indivs      = mean(n_individuals),
    median_indivs    = median(n_individuals),
    nr_quad_above3   = sum( n_individuals <= 3),
    prop_quad_above3 = mean(n_individuals <= 3))
  

# Survival ---------------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives), !is.na(logvol_t0), !is.na(logvol_t0_2), !is.na(logvol_t0_3),
         !is.na(recruits)) %>%
  filter(size_t0 != 0) %>%
  mutate(recruits = as.factor(recruits)) %>% 
  dplyr::select(id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
                stage, recruits)

fig_su_raw <- ggplot(
  data = plot_binned_prop(df_su, 10, logvol_t0, survives)) +
  geom_jitter(data = df_su, aes(x = logvol_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logvol_t0, y = survives),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logvol_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Survival',
       subtitle = v_ggp_suffix,
       x = expression('log(volume)'[t0]),
       y = expression('Survival to time t1'))
fig_su_raw


# Survival model ---------------------------------------------------------------
mod_su_0  <- glm(survives ~ 1, data = df_su, family = 'binomial') 
mod_su_1  <- glm(survives ~ logvol_t0, data = df_su, family = 'binomial') 
mod_su_2  <- glm(survives ~ logvol_t0 + logvol_t0_2, data = df_su, family = 'binomial')  
mod_su_3  <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_su, family = 'binomial')  
mod_su_01 <- glm(survives ~ recruits, data = df_su, family = 'binomial')
mod_su_11 <- glm(survives ~ logvol_t0 + recruits, data = df_su, family = 'binomial')
mod_su_21 <- glm(survives ~ logvol_t0 + logvol_t0_2 + recruits, data = df_su, family = 'binomial')  
mod_su_31 <- glm(survives ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits, data = df_su, family = 'binomial')
mod_su_12 <- glm(survives ~ logvol_t0 * recruits, data = df_su, family = 'binomial')
mod_su_22 <- glm(survives ~ logvol_t0 * recruits + logvol_t0_2:recruits, data = df_su, family = 'binomial')  
mod_su_32 <- glm(survives ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, data = df_su, family = 'binomial')

mods_su      <- list(
  mod_su_0,  mod_su_1,  mod_su_2,  mod_su_3,
  mod_su_01, mod_su_11, mod_su_21, mod_su_31,
  mod_su_12, mod_su_22, mod_su_32)
mods_su_dAIC <- AICctab(mods_su, weights = T, sort = F)$dAIC

mods_su_sorted       <- order(mods_su_dAIC)
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index       <- mod_su_index_bestfit - 1 
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index       <- v_mod_set_su
  }
mod_su_bestfit       <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef         <- coef(mod_su_bestfit)

df_su_newdata <- df_su %>%
  group_by(recruits) %>%
  summarise(x = list(seq(min(logvol_t0), max(logvol_t0), 
                         length.out = 100)), .groups = "drop") %>%
  unnest(x) %>%
  mutate(
    logvol_t0 = x,
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)) %>% 
  mutate(predicted = predict(mod_su_bestfit, newdata = ., type = "response"))

fig_su_line <- ggplot(
  df_su, aes(x = logvol_t0, y = survives, color = factor(recruits))) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3) +
  geom_line(data = df_su_newdata, aes(y = predicted), size = 1) +
  labs(
    title    = 'Survival probability by volume and recruits',
    subtitle = v_ggp_suffix,
    x        = 'Volumen t0 (log)',
    y        = 'Probability of survival',
    color    = 'Recruits') +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  theme_bw() +
  theme(legend.position = 'none')

df_su_bindata <- df_su %>%
  group_by(recruits) %>%
  group_modify(~ plot_binned_prop(.x, 10, logvol_t0, survives)) %>%
  ungroup() %>%
  mutate(recruits = factor(recruits))

df_su_pred <- df_su %>%
  group_by(recruits) %>%
  reframe(logvol_t0 = seq(min(logvol_t0), max(logvol_t0), length.out = 100)) %>%
  mutate(
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)) %>%
  mutate(survives = predict(mod_su_bestfit, newdata = ., type = 'response'))

fig_su_bin <- ggplot() +
  geom_point(
    data = df_su_bindata, aes(x = logvol_t0, y = survives, color = recruits)) +
  geom_errorbar(aes(x = logvol_t0, ymin = lwr, ymax = upr, color = recruits),
                data = df_su_bindata, width = 0.1) +
  geom_line(aes(x = logvol_t0, y = survives, color = recruits),
            data = df_su_pred, size = 1.2) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  labs(x = 'Volumen t0 (log)', y = '', color  = 'Recruit') +
  theme_bw() +
  ylim(0, 1)

fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  filter(
    !is.na(size_t1),      !is.na(size_t0), 
    !is.na(logsize_t0),   !is.na(logsize_t1), 
    !is.na(logsize_t0_2), !is.na(logsize_t0_3),
    !is.na(volume_t0),    !is.na(volume_t1), 
    !is.na(logvol_t0),    !is.na(logvol_t1), 
    !is.na(logvol_t0_2),  !is.na(logvol_t0_3),
    !is.na(recruits)) %>% 
  mutate(recruits = as.factor(recruits)) %>% 
  dplyr::select(
    id, year, size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
    stage, recruits)

fig_gr_raw <- ggplot(
  data  = df_gr, aes(logvol_t0, logvol_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(volume) ' [t0]),
       y        = expression('log(volume)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))
fig_gr_raw


# Growth model -----------------------------------------------------------------
mod_gr_0   <- lm(logvol_t1 ~ 1, data = df_gr)
mod_gr_1   <- lm(logvol_t1 ~ logvol_t0, data = df_gr)
mod_gr_2   <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr)  
mod_gr_3   <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr)
mod_gr_01  <- lm(logvol_t1 ~ recruits, data = df_gr)
mod_gr_11  <- lm(logvol_t1 ~ logvol_t0 + recruits, data = df_gr)
mod_gr_21  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + recruits, data = df_gr)  
mod_gr_31  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits, data = df_gr)
mod_gr_12  <- lm(logvol_t1 ~ logvol_t0 * recruits, data = df_gr)
mod_gr_221 <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2, data = df_gr)  
mod_gr_22  <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits, data = df_gr)  
mod_gr_32  <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, data = df_gr)

mods_gr      <- list(
  mod_gr_0,  mod_gr_1,  mod_gr_2,   mod_gr_3,
  mod_gr_01, mod_gr_11, mod_gr_21,  mod_gr_31,
  mod_gr_12, mod_gr_22, mod_gr_221, mod_gr_32)
mods_gr_dAIC <- AICctab(mods_gr, weights = T, sort = F)$dAIC

mods_gr_sorted       <- order(mods_gr_dAIC)
if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr
}
mod_gr_bestfit       <- mod_gr_221
mod_gr_ranef         <- coef(mod_gr_bestfit)

df_gr_newdata <- df_gr %>%
  group_by(recruits) %>%
  summarise(range = list(seq(min(logvol_t0), max(logvol_t0), length.out = 100)), .groups = 'drop') %>%
  unnest(range) %>%
  mutate(
    logvol_t0 = range,
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)) %>%
  dplyr::select(-range)

df_gr_newdata <- cbind(df_gr_newdata, as.data.frame(predict(
  mod_gr_bestfit,
  newdata = df_gr_newdata,
  interval = 'confidence')))


fig_gr <- ggplot(df_gr, aes(x = logvol_t0, y = logvol_t1, color = recruits)) +
  geom_point(alpha = 0.4) +
  geom_line(data = df_gr_newdata, aes(y = fit), size = 1) +
  geom_ribbon(
    data = df_gr_newdata,
    aes(y = fit, ymin = lwr, ymax = upr, fill = recruits),
    alpha = 0.2,
    color = NA) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  scale_fill_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  labs(
    title = 'Growth prediction',
    subtitle = v_ggp_suffix,
    x = 'Volumen t0 (log)',
    y = 'Volumen t1 (log)',
    color = 'Recruit',
    fill = 'Recruit') + 
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


# Flower data ------------------------------------------------------------------
# We exclude all recruits form the probability to flower
df_fl <- df %>% 
  filter(!is.na(flower), !is.na(logvol_t0), !is.na(logvol_t0_2), !is.na(logvol_t0_3)) %>%
  filter(recruits == 0) %>% 
  dplyr::select(id, year, size_t0, flower, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
                stage)

fig_fl_raw <- ggplot(
  data = plot_binned_prop(df_fl, 10, logvol_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logvol_t0, y = flower), alpha = 0.1, 
              position = position_jitter(width = 0.1, height = 0.3)) +
  geom_point(aes(x = logvol_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logvol_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text     = element_text(size = 8),
        title         = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('Volume (log)'[t0]),
       y        = expression('Flowering probability t0'))
fig_fl_raw


# Flower model -----------------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ 1,
                data = df_fl, family = 'binomial') 
# Logistic regression
mod_fl_1 <- glm(flower ~ logvol_t0,
                data = df_fl, family = 'binomial') 
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logvol_t0 + logvol_t0_2,
                data = df_fl, family = 'binomial')  
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
                data = df_fl, family = 'binomial')  

# Compare models using AIC
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- AICctab(mods_fl, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fl_sorted <- order(mods_fl_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index       <- mod_fl_index_bestfit - 1 
} else {
  mod_fl_index_bestfit <- v_mod_set_fl +1
  v_mod_fl_index       <- v_mod_set_fl
}

mod_fl_bestfit <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef   <- coef(mod_fl_bestfit)

# Generate predictions for survival across a range of sizes
mod_fl_x <- seq(
  min(df_fl$logvol_t0, na.rm = T),
  max(df_fl$logvol_t0, na.rm = T), length.out = 100)

df_fl_pred <- predictor_fun(mod_fl_x, mod_fl_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logvol_t0 = mod_fl_x, flower = .)

fig_fl_line <- ggplot() +
  geom_jitter(data = df_fl, aes(x = logvol_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred, aes(x = logvol_t0, y = flower),
            color = line_color_pred_fun(mod_fl_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Flowering prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_fl_bin <- ggplot() +
  geom_point(data = plot_binned_prop(df_fl, 10, logvol_t0, flower), 
             aes(x = logvol_t0, y = flower)) +
  geom_errorbar(
    data = plot_binned_prop(df_fl, 10, logvol_t0, flower), 
    aes(x = logvol_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_fl_pred, aes(x = logvol_t0, y = flower),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_fl <- fig_fl_line + fig_fl_bin + plot_layout()
fig_fl


# Fecundity --------------------------------------------------------------------
# Conditional on flowering
df_fec <- df %>%
  filter(flower == 1, !is.na(logvol_t0))

df_fec %>% filter(recruits == 1 & flowering_stems > 0)
'there are no recruits that produce a flowering stock in t0'

# Since there are no 0s in the dataset we go for a truncated nb model
df_fec$flowering_stems %>% summary()
'I couldnt find a functioning truncated nb function'

mod_fe_0 <- glm.nb(flowering_stems ~ 1, data = df_fec)
mod_fe_1 <- glm.nb(flowering_stems ~ logvol_t0, data = df_fec)
mod_fe_2 <- glm.nb(flowering_stems ~ logvol_t0 + logvol_t0_2, data = df_fec)
mod_fe_3 <- glm.nb(flowering_stems ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_fec)

mods_fe      <- list(mod_fe_0, mod_fe_1, mod_fe_2, mod_fe_3)
mods_fe_dAIC <- AICctab(mods_fe, weights = T, sort = F)$dAIC
mods_fe_i_sort <- order(mods_fe_dAIC)

if (length(v_mod_set_fe) == 0) {
  mod_fe_i_best <- mods_fe_i_sort[1]
  v_mod_fe_i    <- mod_fe_i_best - 1 
} else {
  mod_fe_i_best <- v_mod_set_fe +1
  v_mod_fe_i    <- v_mod_set_fe
}

mod_fe_best  <- mods_fe[[mod_fe_i_best]]
mod_fe_ranef <- coef(mod_fe_best)

df_fec_preddata <- tibble(
  logvol_t0 = seq(min(df_fec$logvol_t0, na.rm = TRUE),
                  max(df_fec$logvol_t0, na.rm = TRUE),
                  length.out = 100)) %>%
  mutate(logvol_t0_2 = logvol_t0^2,
         logvol_t0_3 = logvol_t0^3)

mod_fe_pred <- predict(mod_fe_best, newdata = df_fec_preddata, type = 'link', se.fit = TRUE)

df_fec_preddata <- df_fec_preddata %>%
  mutate(
    fit_link = mod_fe_pred$fit,
    se_link = mod_fe_pred$se.fit,
    fit_link_lower = fit_link - 1.96 * se_link,
    fit_link_upper = fit_link + 1.96 * se_link,
    predicted_stems = exp(fit_link),
    predicted_stems_lower = exp(fit_link_lower),
    predicted_stems_upper = exp(fit_link_upper))

fig_fe <- ggplot(df_fec, aes(x = logvol_t0, y = flowering_stems)) +
  geom_jitter(width = 0.1, height = 0.2, alpha = 0.4) +
  geom_line(data = df_fec_preddata, aes(x = logvol_t0, y = predicted_stems), color = 'darkgreen', size = 1.2) +
  geom_ribbon(
    data = df_fec_preddata,
    aes(x = logvol_t0, ymin = predicted_stems_lower, ymax = predicted_stems_upper),
    fill = 'darkgreen', alpha = 0.2,
    inherit.aes = FALSE) +
  labs(
    title    = 'Fecundity',
    subtitle = v_ggp_suffix,
    x        = 'Volume t0 (log)',
    y        = 'Number of Flowering Stems') +
  theme_bw()
fig_fe


# Flowering stock to Recruit transition ----------------------------------------
df_fs2r <- df %>%
  filter(!is.na(flowering_stems)) %>%
  group_by(site, quad, year) %>%
  summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = 'drop') %>%
  mutate(year_recruits = year + 1) %>%
  left_join(
    df %>%
      filter(recruits == 1) %>%
      group_by(site, quad, year) %>%
      summarise(recruit_count = n(), .groups = 'drop'),
    by = c('site', 'quad', 'year_recruits' = 'year')
  ) %>%
  mutate(recruit_count = ifelse(is.na(recruit_count), 0, recruit_count)) %>%
  filter(total_stocks < 100)

df_fs2r %>% 
  group_by(year) %>% 
  summarise(total_stocks = sum(total_stocks))

mod_fs2r <- glm.nb(recruit_count ~ total_stocks, data = df_fs2r)

df_fs2r_newdata <- data.frame(
  total_stocks = seq(min(df_fs2r$total_stocks, na.rm = TRUE),
                     max(df_fs2r$total_stocks, na.rm = TRUE), length.out = 100))

mod_fs2r_preds <- predict(mod_fs2r, newdata = df_fs2r_newdata, type = 'link', se.fit = TRUE)
df_fs2r_newdata$fit <- exp(mod_fs2r_preds$fit)  # inverse link (log)
df_fs2r_newdata$lower <- exp(mod_fs2r_preds$fit - 1.96 * mod_fs2r_preds$se.fit)
df_fs2r_newdata$upper <- exp(mod_fs2r_preds$fit + 1.96 * mod_fs2r_preds$se.fit)

ggplot(df_fs2r, aes(x = total_stocks, y = recruit_count)) +
  geom_jitter(height = 0.2, width = 0.5, alpha = 0.4, color = 'gray40') +
  geom_ribbon(data = df_fs2r_newdata, aes(x = total_stocks, ymin = lower, ymax = upper), inherit.aes = FALSE,
              alpha = 0.2, fill = 'red') +
  geom_line(data = df_fs2r_newdata, aes(x = total_stocks, y = fit), inherit.aes = FALSE,
            color = 'darkred', size = 1.2) +
  labs(
    title = 'Recruits t1 by Flowering Stems t0 - Negative Binomial',
    x = 'Total Flowering Stems t0 (per site/plot/year)',
    y = 'Number of Recruits t1'
  ) +
  theme_bw()


# Flowering stock to Recruit transition - Zero truncated, Bayesian -------------
#install.packages('brms') 
library(brms)

df_fs2r_0t <- df_fs2r %>%
  filter(recruit_count > 0)

# mod_fs2r_0t <- brm(
#   bf(recruit_count | trunc(lb = 1) ~ total_stocks),
#   data = df_fs2r_0t,
#   family = negbinomial(link = 'log'),
#   chains = 4,
#   cores = 4,
#   iter = 2000,
#   control = list(adapt_delta = 0.95)
# )
# saveRDS(mod_fs2r_0t, file = file.path(dir_data, 'bayes_recruit_model.rds'))
mod_fs2r_0t <- readRDS(file.path(dir_data, 'bayes_recruit_model.rds'))
summary(mod_fs2r_0t)
plot(mod_fs2r_0t)
pp_check(mod_fs2r_0t)

df_fs2r_0t_pred <- df_fs2r_0t %>%
  bind_cols(fitted(mod_fs2r_0t, newdata = df_fs2r_0t, re_formula = NA, summary = TRUE))

ggplot(df_fs2r_0t_pred, aes(x = total_stocks, y = recruit_count)) +
  geom_point(alpha = 0.4, color = 'black') + 
  geom_line(aes(y = Estimate), color = 'blue', size = 1) +  
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, fill = 'blue') +
  labs(
    x = 'Total Stocks (year t)',
    y = 'Recruits (year t+1, truncated at >0)',
    title = 'Posterior mean and 95% credible interval'
  ) +
  theme_minimal()


# Recruits by number of individuals per quadrat --------------------------------
ggplot(df %>%
         group_by(site, quad, year) %>%
         summarise(
           n_individuals = n_distinct(id),
           n_recruits = sum(recruits, na.rm = TRUE),
           .groups = 'drop'), aes(x = n_individuals, y = n_recruits)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
  geom_smooth(method = 'glm', method.args = list(family = 'poisson'), color = 'darkred') +
  labs(
    x = 'Number of Individuals per Quadrat-Year',
    y = 'Number of Recruits',
    title = 'Recruitment vs. Density') +
  theme_bw()


# Recruits by total volume per quadrat and fire---------------------------------
ggplot(df %>%
         group_by(site, quad, year) %>%
         summarise(
           total_volume = sum(volume_t0, na.rm = TRUE),
           n_recruits   = sum(recruits, na.rm = TRUE),
           fire = as.factor(
             if (all(is.na(postburn_plant))) {'_NA_ture'} else {
               max(postburn_plant, na.rm = TRUE)}),
           .groups = 'drop'),
       aes(x = total_volume, y = n_recruits, color = fire)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c('0' = 'darkred', '1' = 'forestgreen', '_NA_ture' = 'gray50')) +
  theme_bw() +
  labs(x = 'Total Volume per Quadrat', y = 'Number of Recruits', color = 'Postburn Status')


# Tails from the publication ---------------------------------------------------
# Seedbank
"P. lewtonii plants are usually killed by fire and seedlings recruit
post-fire from a soil seedbank (B. Pace-Aldana, pers. comm.;
C. W. Weekley, unpubl. data)."
"A separate field
experiment demonstrated that buried seeds may retain viability
for at least 2 years (C. W. Weekley and E. S. Menges, unpubl.
data), indicating that P. lewtonii is capable of accumulating a
persistent (sensuThompson and Grime 1979) soil seedbank. The
sudden appearance of large cohorts of seedling recruits suggests
that the persistent seedbanks of P. lewtonii are ecologically
relevant."
" In the March 2002 census,
burned quadrats out-recruited unburned quadrats 3.7–1
(184.0% v. 49.2%; c2 = 48.880, d.f. = 1, P < 0.001). Over the 2
censuses, burned quadrats had proportional recruitment of
between 75.8% and 79.0% of seedlings recruited into burned
quadrats with no surviving reproductive adults, indicating that
most recruitment was from the preburn seedbank rather than
from newly produced seeds."


# Recruitment data -------------------------------------------------------------
df_re <- df %>%
  group_by(year, site, quad) %>%
  summarise(tot_p_volume = sum(volume_t0, na.rm = TRUE), .groups = 'drop') %>%
  {
    df_quad <- .
    df_group <- df_quad %>%
      group_by(year) %>%
      summarise(g_cov = mean(tot_p_volume), .groups = 'drop')
    
    df_cover <- left_join(df_quad, df_group, by = 'year') %>%
      mutate(year = as.integer(year + 1)) %>%
      drop_na()
    
    df_re <- df %>%
      group_by(year, site, quad) %>%
      summarise(nr_recs = sum(recruits, na.rm = TRUE), .groups = 'drop')
    
    left_join(df_cover, df_re, by = c('year', 'site', 'quad'))
  }

ggplot(
  df_re, aes(x = tot_p_volume, y = nr_recs)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant volume '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))

# Density dependency
df_re_qd <- df %>% 
  group_by(site, quad, year) %>%
  dplyr::select(recruits) %>% 
  summarise(rec_qd_t1 = sum(recruits, na.rm = T)) %>%
  left_join(df %>% 
              group_by(site, quad, year) %>% 
              summarise(nr_ind = sum(!is.na(volume_t0))) %>% 
              mutate(year = year - 1),
            by = c('site', 'quad', 'year'))

fig_re_dens <- ggplot(data = df_re_qd) + 
  geom_jitter(aes(y = rec_qd_t1, x = nr_ind)) + 
  geom_smooth(aes(y = rec_qd_t1, x = nr_ind), method = 'lm') + 
  theme_bw() + 
  labs(title    = 'Recruitment - desity dependence: Quad level',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant volume '[t0]),   
       y        = expression('Number of recruits '       [t1])) +
  theme(plot.subtitle = element_text(size = 8))

ggsave(file.path(dir_result, 'mean_rec_density_dependency.png'), 
       plot = fig_re_dens, width = 10, height = 5, dpi = 300)


# Recruitment model ------------------------------------------------------------
df_re_mod <- df_re %>% filter(!is.na(nr_recs))
# Fit a negative binomial model for recruitment
mod_rec <- MASS::glm.nb(nr_recs ~ 1, data = df_re_mod)

# Generate predictions for recruitment
df_re_mod <- df_re_mod %>% 
  mutate(mod_pred = predict(mod_rec, type = 'response')) 


# Per-capita reproduction ------------------------------------------------------
df_repr_pc <- df %>%
  filter(!is.na(volume_t0)) %>% 
  summarize(n_adults = n()) %>%
  bind_cols(
    df_re_mod %>%
      summarize(nr_recs = sum(nr_recs, na.rm = TRUE),
                mod_pred = sum(mod_pred, na.rm = TRUE))) %>%
  mutate(repr_pc_mean = mod_pred / n_adults,
         repr_pc_obs  = nr_recs  / n_adults) %>%
  drop_na()

# overall level
df %>% 
  filter(!is.na(volume_t0)) %>% 
  summarize(n_adults = n()) %>% 
  bind_cols(df %>% 
              filter(!is.na(recruits)) %>%
              summarize(n_rec = n())) %>% 
  mutate(rp_pc_m = n_rec / n_adults)


# site level
df %>% 
  filter(!is.na(volume_t0)) %>%
  group_by(site) %>% 
  summarize(n_adults = n()) %>% 
  left_join(df %>% 
              filter(!is.na(recruits)) %>%
              group_by(site) %>%
              summarize(n_rec = n()),
            by = 'site') %>% 
  mutate(n_rec = ifelse(is.na(n_rec), 0, n_rec)) %>%
  mutate(repr_pc_mean = n_rec / n_adults) %>% 
  summarise(n_adults = sum(n_adults),
            n_rec    = sum(n_rec),
            rp_pc_m  = mean(repr_pc_mean))


# quad level
df %>% 
  filter(!is.na(volume_t0)) %>%
  group_by(quad) %>% 
  summarize(n_adults = n()) %>% 
  left_join(df %>% 
              filter(!is.na(recruits)) %>%
              group_by(quad) %>%
              summarize(n_rec = n()),
            by = 'quad') %>% 
  mutate(n_rec = ifelse(is.na(n_rec), 0, n_rec)) %>%
  mutate(repr_pc_mean = n_rec / n_adults) %>% 
  summarise(n_adults = sum(n_adults),
            n_rec    = sum(n_rec),
            rp_pc_m  = mean(repr_pc_mean))


# Extracting parameter estimates -----------------------------------------------
# Survival
coef_su_fe  <- data.frame(coefficient = names(coef(mod_su_bestfit)),
                          value       =       coef(mod_su_bestfit))

coef_su <- Reduce(function(...) rbind(...), list(coef_su_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Growth
coef_gr_fe  <- data.frame(coefficient = names(coef(mod_gr_bestfit)),
                          value       =       coef(mod_gr_bestfit))
coef_gr_var <- data.frame(coefficient = names(coef(mod_gr_var)),
                          value       =       coef(mod_gr_var))

coef_gr <- Reduce(function(...) rbind(...), list(coef_gr_fe, coef_gr_var)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Flower
coef_fl_fe  <- data.frame(coefficient = names(coef(mod_fl_bestfit)),
                          value       =       coef(mod_fl_bestfit))

coef_fl <- Reduce(function(...) rbind(...), list(coef_fl_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Fruit
coef_fe_fe  <- data.frame(coefficient = names(coef(mod_fe_best)),
                          value       =       coef(mod_fe_best))

coef_fr <- Reduce(function(...) rbind(...), list(coef_fe_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Recruitment 
df_re_size <- df %>% subset(recruits == 1)

# # Miscellany
# coef_misc   <- data.frame(coefficient = c('rec_siz', 'rec_sd',
#                                           'fecu_b0', 
#                                           'max_siz', 'min_siz'),
#                           value       = c(mean(log(df_re_size$size_t0), na.rm = T), 
#                                           sd(  log(df_re_size$size_t0), na.rm = T),
#                                           repr_pc_mean_stock,
#                                           df_gr$logsize_t0 %>% max, 
#                                           df_gr$logsize_t0 %>% min))
# 
# extr_value <- function(x, field){
#   subset(x, coefficient == field)$value
# }
# 
# pars <- Filter(function(x) length(x) > 0, list(
#   prefix      = v_script_prefix,
#   species     = v_species,
#   surv_b0     = extr_value(coef_su, 'b0'),
#   surv_b1     = extr_value(coef_su, 'logvol_t0'),
#   surv_b2     = extr_value(coef_su, 'logvol_t0_2'),
#   surv_b3     = extr_value(coef_su, 'logvol_t0_3'),
#   surv_br     = extr_value(coef_su, 'recruits1'),
#   grow_b0     = extr_value(coef_gr, 'b0'),
#   grow_b1     = extr_value(coef_gr, 'logvol_t0'),
#   grow_b2     = extr_value(coef_gr, 'logvol_t0_2'),
#   grow_b3     = extr_value(coef_gr, 'logvol_t0_3'),
#   grow_br     = extr_value(coef_gr, 'recruits1'),
#   grow_b1_br  = extr_value(coef_gr, 'logvol_t0:recruits1'),
#   grow_b2_br0 = extr_value(coef_gr, 'recruits0:logvol_t0_2'),
#   grow_b2_br1 = extr_value(coef_gr, 'recruits0:logvol_t0_2'),
#   a           = extr_value(coef_gr, 'a'),
#   b           = extr_value(coef_gr, 'b'),
#   fl_b0       = extr_value(coef_fl, 'b0'),
#   fl_b1       = extr_value(coef_fl, 'logvol_t0'),
#   fl_b2       = extr_value(coef_fl, 'logvol_t0_2'),
#   fl_b3       = extr_value(coef_fl, 'logvol_t0_3'),
#   fr_b0       = extr_value(coef_fl, 'b0'),
#   fr_b1       = extr_value(coef_fr, 'logvol_t0'),
#   fr_b2       = extr_value(coef_fr, 'logvol_t0_2'),
#   fr_b3       = extr_value(coef_fr, 'logvol_t0_3'),
#   fecu_b0     = extr_value(coef_misc, 'fecu_b0'),
#   recr_sz     = extr_value(coef_misc, 'rec_siz'),
#   recr_sd     = extr_value(coef_misc, 'rec_sd'),
#   L           = extr_value(coef_misc, 'min_siz'),
#   U           = extr_value(coef_misc, 'max_siz'),
#   mat_siz     = 200,
#   mod_su_index = v_mod_su_index,
#   mod_gr_index = v_mod_gr_index,
#   mod_gr_index = v_mod_fl_index
# ))
# 
# 
# # Building the IPM -------------------------------------------------------------
# # Function describing the invert logit
# inv_logit <- function(x) {exp(x) / (1 + exp(x))}
# 
# # Survival of x-sized individual to time t1
# sx <- function(x, pars, num_pars = v_mod_su_index) {
#   survival_value <- pars$surv_b0
#   for (i in 1:num_pars) {
#     param_name <- paste0('surv_b', i)
#     if (!is.null(pars[[param_name]])) {
#       survival_value <- survival_value + pars[[param_name]] * x^(i)
#     }
#   }
#   return(inv_logit(survival_value))
# }
# 
# # Function describing standard deviation of growth model
# grow_sd <- function(x, pars) {
#   pars$a * (exp(pars$b* x)) %>% sqrt 
# }
# 
# # Growth from size x to size y
# gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
#   mean_value <- 0
#   for (i in 0:num_pars) {
#     param_name <- paste0('grow_b', i)
#     if (!is.null(pars[[param_name]])) {
#       mean_value <- mean_value + pars[[param_name]] * x^i
#     }
#   }
#   sd_value <- grow_sd(x, pars)
#   return(dnorm(y, mean = mean_value, sd = sd_value))
# }
# 
# # Function describing the transition kernel
# pxy <- function(x, y, pars) {
#   return(sx(x, pars) * gxy(x, y, pars))
# }
# 
# # Flowering of x-sized individual at time t0
# fl_x <- function(x, pars, num_pars = v_mod_fl_index) {
#   val <- pars$fl_b0
#   for (i in 1:num_pars) {
#     param <- paste0('fl_b', i)
#     if (!is.null(pars[[param]])) {
#       val <- val + pars[[param]] * x^i
#     }
#   }
#   inv_logit(val)
# }
# 
# # Fruiting of x-sized individuals at time t0
# fr_x <- function(x, pars, num_pars = v_mod_fr_i) {
#   val <- pars$fr_b0
#   for (i in 1:num_pars) {
#     param <- paste0('fr_b', i)
#     if (!is.null(pars[[param]])) {
#       val <- val + pars[[param]] * x^i
#     }
#   }
#   exp(val)  # Negative binomial uses log link
# }
# 
# # Recruitment size distribution at time t1
# re_y_dist <- function(y, pars) {
#   dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
# }
# 
# # F-kernel
# fyx <- function(y, x, pars) {
#   fl_x(x, pars) *
#     fr_x(x, pars) *
#     pars$fecu_b0 *
#     re_y_dist(y, pars)
# }
# 
# # Kernel
# kernel <- function(pars) {
#   
#   # number of bins over which to integrate
#   n   <- pars$mat_siz 
#   # lower limit of integration
#   L   <- pars$L  
#   # upper limit of integration
#   U   <- pars$U       
#   # bin size
#   h   <- (U - L) / n  
#   # lower boundaries of bins
#   b   <- L + c(0:n) * h             
#   # midpoints of bins
#   y   <- 0.5 * (b[1:n] + b[2:(n + 1)]) 
#   
#   # Survival vector
#   Smat   <- c()
#   Smat   <- sx(y, pars)
#   
#   # Growth matrix
#   Gmat   <- matrix(0, n, n)
#   Gmat[] <- t(outer(y, y, gxy, pars)) * h
#   
#   # Growth/survival transition matrix
#   Tmat   <- matrix(0, n, n)
#   
#   # Correct for eviction of offspring
#   for(i in 1:(n / 2)) {
#     Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
#     Tmat[,i]  <- Gmat[,i] * Smat[i]
#   }
#   
#   # Correct eviction of large adults
#   for(i in (n / 2 + 1):n) {
#     Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
#     Tmat[,i]  <- Gmat[,i] * Smat[i]
#   }
#   
#   # Fertility matrix
#   Fmat <- outer(y, y, Vectorize(function(x, y) fyx(x, y, pars))) * h
#   
#   # Full Kernel is simply a summation of fertility and transition matrices
#   k_yx <- Fmat + Tmat
#   
#   return(list(k_yx    = k_yx,
#               Fmat    = Fmat,
#               Tmat    = Tmat,
#               Gmat    = Gmat,
#               meshpts = y))
# }
# 
# lambda_ipm <- function(i) {
#   return(Re(eigen(kernel(i)$k_yx)$value[1]))
# }
# 
# # mean population growth rate
# lam_mean <- lambda_ipm(pars)
# lam_mean
# 
# 
# # Observed population growth ---------------------------------------------------
# df_counts_year <- df %>%
#   group_by(year) %>%
#   filter(!is.na(survives)) %>% 
#   summarise(n = n())
# 
# # Then compute observed lambda
# lam_obs_y <- df_counts_year$n[-1] / df_counts_year$n[-nrow(df_counts_year)]
# lam_obs_mean <- mean(lam_obs_y, na.rm = TRUE)
# 
# 
# # IPM investigation ------------------------------------------------------------
# summary(df$size_t0, na.rm = TRUE)
# hist(df$size_t0, na.rm = TRUE)
# 
# # Use actual observed size limits
# min_x <- min(df$size_t0, na.rm = TRUE)
# max_x <- max(df$size_t0, na.rm = TRUE)
# 
# # Recalculate mesh points
# n_mesh <- 200  # resolution
# x_vals <- seq(min_x, max_x, length.out = n_mesh)
# y_vals <- x_vals  # assuming y follows same range
# 
# 
# fl_vals <- sapply(x_vals, function(x) fl_x(x, pars))
# plot(x_vals, fl_vals, type = "l", main = "Flowering Probability vs Size",
#      xlab = "Size (x)", ylab = "Flowering Probability")
# 
# 
# fr_vals <- sapply(x_vals, function(x) fr_x(x, pars))
# plot(x_vals, fr_vals, type = "l", main = "Fruiting Counts vs Size",
#      xlab = "Size (x)", ylab = "Mean Number of Fruits")
# 
# 
# re_vals <- sapply(y_vals, function(y) re_y_dist(y, pars))
# plot(y_vals, re_vals, type = "l", main = "Recruitment Size Distribution",
#      xlab = "Size (y)", ylab = "Density")
# 
# 
# fyx_vals <- outer(x_vals, y_vals, Vectorize(function(x, y) fyx(x, y, pars)))
# selected_x <- c(2, 5, 8)  # example sizes
# matplot(y_vals, t(fyx_vals[selected_x, ]), type = "l", lty = 1, col = 1:length(selected_x),
#         main = "Kernel fyx across y for selected x",
#         xlab = "Size (y)", ylab = "fyx value")
# legend("topright", legend = paste("x =", selected_x), col = 1:length(selected_x), lty = 1)
# 
# 
# # Approximated lambda
# dx <- x_vals[2] - x_vals[1]
# dy <- y_vals[2] - y_vals[1]
# lambda_approx <- sum(fyx_vals) * dx * dy
# 
# cat("F-only lambda (approx):", lambda_approx, "\n")
# cat("Full IPM lambda (eigen):", lam_mean, "\n")
# 
# 
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
