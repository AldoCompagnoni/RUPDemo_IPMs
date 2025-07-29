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
mod_su_index_bestfit <- mods_su_sorted[1]
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
  data  = df_gr, aes(logvol_t0, logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))
fig_gr_raw


# Growth model -----------------------------------------------------------------
mod_gr_0  <- lm(logvol_t1 ~ 1, data = df_gr)
mod_gr_1  <- lm(logvol_t1 ~ logvol_t0, data = df_gr)
mod_gr_2  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2, data = df_gr)  
mod_gr_3  <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3, data = df_gr)
mod_gr_01 <- lm(logvol_t1 ~ recruits, data = df_gr)
mod_gr_11 <- lm(logvol_t1 ~ logvol_t0 + recruits, data = df_gr)
mod_gr_21 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + recruits, data = df_gr)  
mod_gr_31 <- lm(logvol_t1 ~ logvol_t0 + logvol_t0_2 + logvol_t0_3 + recruits, data = df_gr)
mod_gr_12 <- lm(logvol_t1 ~ logvol_t0 * recruits, data = df_gr)
mod_gr_22 <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits, data = df_gr)  
mod_gr_32 <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, data = df_gr)

mods_gr      <- list(
  mod_gr_0,  mod_gr_1,  mod_gr_2,  mod_gr_3,
  mod_gr_01, mod_gr_11, mod_gr_21, mod_gr_31,
  mod_gr_12, mod_gr_22, mod_gr_32)
mods_gr_dAIC <- AICctab(mods_gr, weights = T, sort = F)$dAIC

mods_gr_sorted       <- order(mods_gr_dAIC)
mod_gr_index_bestfit <- mods_gr_sorted[1]
mod_gr_bestfit       <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef         <- coef(mod_gr_bestfit)

df_gr_newdata <- df_gr %>%
  group_by(recruits) %>%
  summarise(range = list(seq(min(logvol_t0), max(logvol_t0), length.out = 100)), .groups = 'drop') %>%
  unnest(range) %>%
  mutate(
    logvol_t0 = range,
    logvol_t0_2 = logvol_t0^2,
    recruits = factor(recruits)
  ) %>%
  dplyr::select(-range) %>%
  mutate(predicted = predict(mod_gr_bestfit, newdata = .))

fig_gr <- ggplot(df_gr, aes(x = logvol_t0, y = logvol_t1, color = recruits)) +
  geom_point(alpha = 0.4) +
  geom_line(data = df_gr_newdata, aes(y = predicted), size = 1) +
  scale_color_manual(values = c('0' = '#BBB857', '1' = '#3666DC')) +
  labs(
    title    = 'Growth prediction',
    subtitle = v_ggp_suffix,
    x = 'Volumen t0 (log)',
    y = 'Volumen t1 (log)',
    color = 'Recruit') +
  theme_minimal()

fig_gr


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
  filter(flower == 1)

# Since there are no 0s in the dataset we go for a truncated nb model
df_fec$flowering_stems %>% summary()
'I couldnt find a functioning truncated nb function'

mod_fe_nb <- glm.nb(flowering_stems ~ logvol_t0, data = df_fec)

df_fec_preddata <- tibble(
  logvol_t0 = seq(min(df_fec$logvol_t0, na.rm = TRUE),
                  max(df_fec$logvol_t0, na.rm = TRUE),
                  length.out = 100))
df_fec_preddata$predicted_stems <- predict(mod_fe_nb, newdata = df_fec_preddata, type = 'response')

# Plot
fig_fe <- ggplot(df_fec, aes(x = logvol_t0, y = flowering_stems)) +
  geom_jitter(width = 0.1, height = 0.2, alpha = 0.4) +
  geom_line(data = df_fec_preddata, aes(y = predicted_stems), color = 'darkgreen', size = 1.2) +
  labs(
    title    = 'Fecundity',
    subtitle = v_ggp_suffix,
    x        = 'Volume t0 (log)',
    y        = 'Number of Flowering Stems') +
  theme_minimal()
fig_fe

# Flowering stock to Recruit transition ----------------------------------------
stocks_by_plot <- df %>%
  filter(!is.na(flowering_stems)) %>%
  group_by(site, quad, year) %>%
  summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = 'drop') %>%
  mutate(year_recruits = year + 1)

recruits_by_plot <- df %>%
  filter(recruits == 1) %>%
  group_by(site, quad, year) %>%
  summarise(recruits_present = 1, .groups = 'drop')  # binary presence

stocks_to_recruits <- stocks_by_plot %>%
  left_join(recruits_by_plot, by = c('site', 'quad', 'year_recruits' = 'year')) %>%
  mutate(recruits_present = ifelse(is.na(recruits_present), 0, recruits_present))

mod_stocks_to_recruits <- glm(recruits_present ~ total_stocks, data = stocks_to_recruits, family = binomial)
summary(mod_stocks_to_recruits)

ggplot(stocks_to_recruits, aes(x = total_stocks, y = recruits_present)) +
  geom_jitter(height = 0.05, alpha = 0.4) +
  stat_smooth(method = 'glm', method.args = list(family = 'binomial'), color = 'blue', se = TRUE) +
  labs(
    title = 'Probability of Recruitment by Total Flowering Stems',
    x = 'Total Flowering Stems (per site/plot/year)',
    y = 'Recruitment Presence (Next Year)'
  ) +
  theme_minimal()


mod_log <- glm(recruits_present ~ log1p(total_stocks), data = stocks_to_recruits, family = binomial)
summary(mod_log)




# Aggregate flowering stems by year and shift forward
stock_to_recruit_by_year <- {
  
  # Total number of flowering stems per year
  stocks_by_year <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE)) %>%
    mutate(year = year + 1)  # Shift forward: flowering stems in year → recruits in year + 1
  
  # Count recruits by year
  recruits_by_year <- df %>%
    filter(recruit == 1) %>%
    group_by(year) %>%
    summarise(n_recruits = n())
  
  # Combine and compute per-stock recruitment
  recruits_by_year %>%
    left_join(stocks_by_year, by = 'year') %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)}

stock_to_recruit_by_year %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE))

hist(stock_to_recruit_by_year$recruits_per_stock,
     main = 'Recruitment Efficiency (Recruits per Flowering Stem)',
     xlab = 'Recruits per Stock', col = 'lightblue', border = 'white')


repr_pc_mean_stock   <- mean(stock_to_recruit_by_year$recruits_per_stock, na.rm = TRUE)
repr_pc_median_stock <- median(stock_to_recruit_by_year$recruits_per_stock, na.rm = TRUE)


stock_to_recruit_by_year_site <- {
  # Total flowering stems by site and year
  stocks_by_year_site <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(site, year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = 'drop') %>%
    mutate(year = year + 1)
  
  # Recruits by site and year
  recruits_by_year_site <- df %>%
    filter(recruit == 1) %>%
    group_by(site, year) %>%
    summarise(n_recruits = n(), .groups = 'drop')
  
  # Join and compute ratio
  recruits_by_year_site %>%
    left_join(stocks_by_year_site, by = c('site', 'year')) %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)}

stock_to_recruit_by_year_plot <- {
  # Flowering stems per plot and year
  stocks_by_year_plot <- df %>%
    filter(!is.na(flowering_stems)) %>%
    group_by(site, quad, year) %>%
    summarise(total_stocks = sum(flowering_stems, na.rm = TRUE), .groups = 'drop') %>%
    mutate(year = year + 1)
  
  # Recruits per plot and year
  recruits_by_year_plot <- df %>%
    filter(recruit == 1) %>%
    group_by(site, quad, year) %>%
    summarise(n_recruits = n(), .groups = 'drop')
  
  # Join and compute
  recruits_by_year_plot %>%
    left_join(stocks_by_year_plot, by = c('site', 'quad', 'year')) %>%
    mutate(recruits_per_stock = n_recruits / total_stocks)}

stock_to_recruit_by_year_site <- stock_to_recruit_by_year_site %>%
  filter(!is.na(total_stocks) & total_stocks > 0)

stock_to_recruit_by_year_plot <- stock_to_recruit_by_year_plot %>%
  filter(!is.na(total_stocks) & total_stocks > 0)

# Summary for site-level
stock_to_recruit_by_year_site %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE))

# Summary for plot-level
stock_to_recruit_by_year_plot %>%
  summarise(
    mean   = mean(recruits_per_stock, na.rm = TRUE),
    sd     = sd  (recruits_per_stock, na.rm = TRUE),
    median = median(recruits_per_stock, na.rm = TRUE))

