# IPM mean - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.06.02


# Website    : https://portal.edirepository.org/nis/mapbrowse?packageid=edi.219.1
# Publication: https://bioone.org/journals/southeastern-naturalist/volume-15/issue-3/058.015.0318/Ecology-and-Conservation-of-the-Endangered-Legume-Crotalaria-avonensis-in/10.1656/058.015.0318.short

# rm(list = ls())


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(
  # negative binomial modeling
  MASS,
  # load tidyverse after MASS to not mask the select function
  tidyverse,
  # bbmle is for AICtab
  bbmle,
  # patchwork plot alingment
  patchwork,
  # binom.cofint for the survival plot
  binom,
  skimr,
  lubridate) # , skimr, ipmr, binom, janitor, lme4


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Crotalaria avonensis')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()


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
v_mod_set_gr <- c(2)
# fig_gr
v_mod_set_su <- c()
# fig_su
v_mod_set_fl <- c()
# fig_fl
v_mod_set_fr <- c()


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
df_og <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data_v2.csv'), col_types = cols(
  burnA = col_double(),
  burnB = col_double(),
  burnC = col_double(),
  burnD = col_double(),
  burnE = col_double(),
  burnF = col_double())) %>% 
  janitor::clean_names() %>%  
  mutate(
    plant_id = as.factor(paste(site, quad, plant, sep = '_')),
    quad_id  = as.factor(paste(site, quad, sep = '_')),
    year     = as.numeric(substr(date, 1, 4)),  
    month    = as.numeric(substr(date, 6, 7)),
    site     = as.factor(site),
    quad     = as.factor(quad),
    mp       = as.factor(mp),
    plant    = as.factor(plant),
    caged    = as.factor(caged),
    veg      = as.factor(veg)) %>%
  arrange(site, quad, quad_id, plant, plant_id, year, month)


df_meta <- data.frame(variable = colnames(df_og)) %>% 
  mutate(definition = c(
    'study site', 'quadrat number',	'macroplot number',	
    'plant number (within quad)', 'direction within circular quad',	
    'distance from quad center',	'was quad caged from 2012 onward',	
    'vegetation type', 'year quadrat was initiated', 
    'year-month of observation',	'fire severity Dec. 2014 or Jan. 2015',
    'fire severity Aug. 2005',	'fire severity May/June 2009',
    "fire severity B's ridge 2016", 'fire severity Feb. 2017',
    'fire severity Oct. 2017', 'survival code for month', 'number of stems',
    'number of branch tips', 'number of flowers (corolla showing)', 
    'number of developing fruits', 'number of mature fruits', 'herbivory code',
    'plant identification', 'quadrat identification', 'sample year', 
    'sample month'))

df_og_extended <- df_og %>%
  # Survival:
  #  is the plant not appearing but it is 3 years before latest record, then NA
  mutate(date = ymd(paste0(date, "-01"))) %>%
  group_by(plant_id) %>%
  mutate(
    # s = 6 and 5 mean death (I think)
    latest_alive_date      = max(date[(s > 0 & s < 4)], na.rm = TRUE),
    earliest_recorded_date = min(date[ s > 0]         , na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    survival = case_when(
      #  is it before the earliest recorded data, then NA
      date <  earliest_recorded_date ~ NA_real_,
      #  is the plant still alive some time after, then 1
      date <  latest_alive_date & !is.infinite(latest_alive_date) ~ 1,
      #  is the plant not alive some time after, then 0
      date == latest_alive_date ~ 0,
      # is it after the latest date that the plant is alive, then NA
      date >  latest_alive_date ~ NA_real_,
      latest_alive_date == 2017 ~ NA_real_)) %>%
  # Recruits
  mutate(
    recruit = case_when(
      s == 3 ~ 1,
      date == earliest_recorded_date & year(earliest_recorded_date) > 1999 + 4 ~ 1,
      TRUE ~ NA_real_))


# Mean data frame --------------------------------------------------------------
df <- df_og_extended %>%
  group_by(site, quad_id , plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(survival )), NA_real_, min(survival,  na.rm = T)),
    size_t0  = if_else(all(is.na(br)),        NA_real_, max(br,        na.rm = T)),
    flower   = if_else(all(is.na(fl)),        NA_real_, max(fl,        na.rm = T)),  
    fruit    = if_else(all(is.na(fr)),        NA_real_, max(fr,        na.rm = T)),
    recruit  = if_else(all(is.na(recruit)),   NA_real_, min(recruit,   na.rm = T)),
    fire_sev = if_else(
      all(is.na(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f))),
      NA_real_, 
      mean(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f), na.rm = TRUE)),
    .groups  = 'drop'
  ) %>% 
  ungroup() %>% 
  mutate(survives = if_else(survives == 0 & year > 2017 - 4, NA, survives))%>%
  # Define dormancy
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0)   ~ 1,
      size_t0  >  0 & !is.na(survives) ~ 0, 
      TRUE ~ NA_real_ 
    ),
    # Generate a new column 'dormancy_count' that counts consecutive 1s
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 & 
        lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1 & lag(dormancy, 5) == 1 ~ 6,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1                         ~ 5,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1                                                 ~ 4,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1           ~ 3,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 0           ~ 2,
      dormancy == 1                                                           ~ 1,
      TRUE ~ dormancy)) %>%
  group_by(site, quad_id, plant_id) %>%
  # Set size_t1 based on survival; propagate size_t0 if survives, otherwise set to NA
  mutate(
    size_t1 = case_when(
      survives == 1 ~ lead(size_t0), 
      TRUE ~ NA_real_)) %>%
  ungroup() %>% 
  # Compute log-transformed sizes and their powers for modeling
  mutate(
    logsize_t0   = log(size_t0),     
    logsize_t1   = log(size_t1),    
    logsize_t0_2 = logsize_t0^2,     
    logsize_t0_3 = logsize_t0^3)


# Survival data ----------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  dplyr::select(plant_id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

fig_su_overall <- ggplot(
  data = plot_binned_prop(df, 10, logsize_t0, survives)) +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logsize_t0, y = survives),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Survival',
       subtitle = v_ggp_suffix,
       x = expression('log(size)'[t0]),
       y = expression('Survival to time t1'))
fig_su_overall 


# Survival model ---------------------------------------------------------------
# Logistic regression
mod_su_0 <- glm(survives ~ 1,
                data = df_su, family = 'binomial') 
# Logistic regression
mod_su_1 <- glm(survives ~ logsize_t0,
                data = df_su, family = 'binomial') 
# Quadratic logistic model
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2,
                data = df_su, family = 'binomial')  
# Cubic logistic model
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                data = df_su, family = 'binomial')  


# Compare models using AIC
mods_su      <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- AICtab(mods_su, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_su_sorted <- order(mods_su_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index       <- mod_su_index_bestfit - 1 
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index       <- v_mod_set_su
}

mod_su_bestfit   <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef         <- coef(mod_su_bestfit)

# Generate predictions for survival across a range of sizes
mod_su_x <- seq(
  min(df_su$logsize_t0, na.rm = T),
  max(df_su$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_su_pred <- predictor_fun(mod_su_x, mod_su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_su_x, survives = .)

# Survival plots
fig_su_line <- ggplot() +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives),
            color = line_color_pred_fun(mod_su_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_su_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su

ggsave(file.path(dir_result, 'mean_survival.png'), 
       plot = fig_su, width = 10, height = 5, dpi = 300)


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  dplyr::select(plant_id, year, size_t0, size_t1,
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

ggplot(
  data  = df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))


# Growth model -----------------------------------------------------------------
# Intercept model 
mod_gr_0 <- lm(logsize_t1 ~ 1,
               data = df_gr)
# Linear model
mod_gr_1 <- lm(logsize_t1 ~ logsize_t0, 
               data = df_gr)
# Quadratic model
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2, 
               data = df_gr)  
# Cubic model
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
               data = df_gr)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- AICtab(mods_gr, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_gr_sorted <- order(mods_gr_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr
}

mod_gr_bestfit         <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef           <- coef(mod_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_gr$pred <- predict(mod_gr_bestfit, type = 'response')

# Growth plot
fig_gr_line <- ggplot(
  df_gr, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, mod_gr_ranef), 
                color = line_color_pred_fun(mod_gr_ranef), lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_pred <- ggplot(
  df_gr, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1), color = 'red', lwd = 2) + 
  theme_bw()

fig_gr <- fig_gr_line + fig_gr_pred + plot_layout() 
fig_gr

ggsave(file.path(dir_result, 'mean_growth.png'), 
       plot = fig_gr, width = 10, height = 5, dpi = 300)


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
df_fl <- df %>% 
  filter(!is.na(flower)) %>%
  dplyr::select(plant_id, year, size_t0, flower, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>% 
  mutate(flower = if_else(flower > 0, 1, flower))

fig_fl_overall <- ggplot(
  data = plot_binned_prop(df_fl, 10, logsize_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower), alpha = 0.1, 
              position = position_jitter(width = 0.1, height = 0.3)) +
  geom_point(aes(x = logsize_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text     = element_text(size = 8),
        title         = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Flowering probability in t0'))
fig_fl_overall


# Flower model -----------------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ 1,
                data = df_fl, family = 'binomial') 
# Logistic regression
mod_fl_1 <- glm(flower ~ logsize_t0,
                data = df_fl, family = 'binomial') 
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logsize_t0 + logsize_t0_2,
                data = df_fl, family = 'binomial')  
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                data = df_fl, family = 'binomial')  


# Compare models using AIC
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- AICtab(mods_fl, weights = T, sort = F)$dAIC

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
  min(df_fl$logsize_t0, na.rm = T),
  max(df_fl$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_fl_pred <- predictor_fun(mod_fl_x, mod_fl_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_fl_x, flower = .)

# Survival plots
fig_fl_line <- ggplot() +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower),
            color = line_color_pred_fun(mod_fl_ranef), lwd = 2) +  
  theme_bw() + 
  labs(title    = 'Flowering prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_fl_bin <- ggplot() +
  geom_point(data = plot_binned_prop(df_fl, 10, logsize_t0, flower), 
             aes(x = logsize_t0, y = flower)) +
  geom_errorbar(
    data = plot_binned_prop(df_fl, 10, logsize_t0, flower), 
    aes(x = logsize_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_fl <- fig_fl_line + fig_fl_bin + plot_layout()
fig_fl

ggsave(file.path(dir_result, 'mean_flowering.png'), 
       plot = fig_fl, width = 10, height = 5, dpi = 300)


# Fruit data -------------------------------------------------------------------
df_fr <- df %>%
  filter(!is.na(fruit), !is.na(logsize_t0))

# Over dispersed -> negative binomial
# otherwise Poisson
mean(df_fr$fruit)
var(df_fr$fruit)


# Fruit model ------------------------------------------------------------------
fr_mod_0 <- glm.nb(fruit ~ 1, data = df_fr)
fr_mod_1 <- glm.nb(fruit ~ logsize_t0, data = df_fr)
fr_mod_2 <- glm.nb(fruit ~ logsize_t0 + logsize_t0_2, data = df_fr)
fr_mod_3 <- glm.nb(fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, data = df_fr)

fr_mods      <- list(fr_mod_0, fr_mod_1, fr_mod_2, fr_mod_3)
fr_mods_dAIC <- AICtab(fr_mods, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
fr_mods_i_sort <- order(fr_mods_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_fr) == 0) {
  fr_mod_i_best <- fr_mods_i_sort[1]
  v_mod_fr_i    <- fr_mod_i_best - 1 
} else {
  fr_mod_i_best <- v_mod_set_fr +1
  v_mod_fr_i    <- v_mod_set_fr
}

fr_mod_best  <- fr_mods[[fr_mod_i_best]]
fr_mod_ranef <- coef(fr_mod_best)

# Plot
{
  # Create prediction range
  logsize_seq <- seq(min(df_fr$logsize_t0, na.rm = TRUE),
                     max(df_fr$logsize_t0, na.rm = TRUE),
                     length.out = 200)
  
  # Build prediction dataframe dynamically
  pred_data <- data.frame(logsize_t0 = logsize_seq)
  
  if ('logsize_t0_2' %in% names(coef(fr_mod_best))) {
    pred_data$logsize_t0_2 <- logsize_seq^2
  }
  if ('logsize_t0_3' %in% names(coef(fr_mod_best))) {
    pred_data$logsize_t0_3 <- logsize_seq^3
  }
  
  # Predict values
  pred_data$pred <- predict(fr_mod_best, newdata = pred_data, type = 'response')
  
  # Create and assign plot
  fig_fr <- ggplot(df_fr, aes(x = logsize_t0, y = fruit)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_line(data = pred_data, aes(x = logsize_t0, y = pred), color = line_color_pred_fun(fr_mod_ranef), linewidth = 1.2) +
    theme_minimal() +
    labs(
      title = 'Predicted Fruit Count vs. log(Size)',
      x = 'log(Size at t0)',
      y = 'Number of Fruits'
    )
}

ggsave(file.path(dir_result, 'mean_fruiting.png'), 
       plot = fig_fr, width = 10, height = 5, dpi = 300)


# Fruit to recruit -------------------------------------------------------------
repr_pc_by_year <- {
  fruit_by_year <- df %>%
    filter(!is.na(fruit)) %>%
    group_by(year) %>%
    summarise(total_fruit = sum(fruit, na.rm = TRUE)) %>%
    mutate(year = year + 1)
  
  recruits_by_year <- df %>%
    filter(recruit == 1) %>%
    group_by(year) %>%
    summarise(n_recruits = n())
  
  recruits_by_year %>%
    left_join(fruit_by_year, by = 'year') %>%
    mutate(repr_pc_mean = n_recruits / total_fruit)
}

repr_pc_by_year %>% 
  summarise(mean(repr_pc_mean),
            sd  (repr_pc_mean),
            median(repr_pc_mean))

hist(repr_pc_by_year$repr_pc_mean)

repr_pc_mean   <- mean(  repr_pc_by_year$repr_pc_mean, na.rm = T)
repr_pc_median <- median(repr_pc_by_year$repr_pc_mean, na.rm = T)


# Recruitment data -------------------------------------------------------------
df_re <- df %>%
  group_by(year, site, quad_id) %>%
  summarise(tot_p_area = sum(size_t0, na.rm = TRUE), .groups = "drop") %>%
  {
    df_quad <- .
    df_group <- df_quad %>%
      group_by(year) %>%
      summarise(g_cov = mean(tot_p_area), .groups = "drop")
    
    df_cover <- left_join(df_quad, df_group, by = "year") %>%
      mutate(year = as.integer(year + 1)) %>%
      drop_na()
    
    df_re <- df %>%
      group_by(year, site, quad_id) %>%
      summarise(nr_quad = sum(recruit, na.rm = TRUE), .groups = "drop")
    
    left_join(df_cover, df_re, by = c("year", "site", "quad_id"))
  }

ggplot(
  df_re, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))

# Density dependency
df_re_qd <- df %>% 
  group_by(site, quad_id, year) %>%
  dplyr::select(recruit) %>% 
  summarise(rec_qd_t1 = sum(recruit, na.rm = T)) %>%
  left_join(df %>% 
              group_by(site, quad_id, year) %>% 
              summarise(nr_ind = sum(!is.na(size_t0))) %>% 
              mutate(year = year - 1),
            by = c('site', 'quad_id', 'year'))

fig_re_dens <- ggplot(data = df_re_qd) + 
  geom_jitter(aes(y = rec_qd_t1, x = nr_ind)) + 
  geom_smooth(aes(y = rec_qd_t1, x = nr_ind), method = 'lm') + 
  theme_bw() + 
  labs(title    = 'Recruitment - desity dependence: Quad level',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))

ggsave(file.path(dir_result, 'mean_rec_density_dependency.png'), 
       plot = fig_re_dens, width = 10, height = 5, dpi = 300)


# Recruitment model ------------------------------------------------------------
df_re_mod <- df_re %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
mod_rec <- MASS::glm.nb(nr_quad ~ 1, data = df_re_mod)

# Generate predictions for recruitment
df_re_mod <- df_re_mod %>% 
  mutate(mod_pred = predict(mod_rec, type = 'response')) 

# Per-capita reproduction ------------------------------------------------------
df_repr_pc <- df %>%
  filter(!is.na(size_t0)) %>% 
  summarize(n_adults = n()) %>%
  bind_cols(
    df_re_mod %>%
      summarize(nr_quad = sum(nr_quad, na.rm = TRUE),
                mod_pred = sum(mod_pred, na.rm = TRUE))
  ) %>%
  mutate(
    repr_pc_mean = mod_pred / n_adults,
    repr_pc_obs = nr_quad / n_adults
  ) %>%
  drop_na()

# overall level
df %>% 
  filter(!is.na(size_t0)) %>% 
  summarize(n_adults = n()) %>% 
  bind_cols(df %>% 
              filter(!is.na(recruit)) %>%
              summarize(n_rec = n())) %>% 
  mutate(rp_pc_m = n_rec / n_adults)


# site level
df %>% 
  filter(!is.na(size_t0)) %>%
  group_by(site) %>% 
  summarize(n_adults = n()) %>% 
  left_join(df %>% 
              filter(!is.na(recruit)) %>%
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
  filter(!is.na(size_t0)) %>%
  group_by(quad_id) %>% 
  summarize(n_adults = n()) %>% 
  left_join(df %>% 
              filter(!is.na(recruit)) %>%
              group_by(quad_id) %>%
              summarize(n_rec = n()),
            by = 'quad_id') %>% 
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
coef_fr_fe  <- data.frame(coefficient = names(coef(fr_mod_best)),
                          value       =       coef(fr_mod_best))

coef_fr <- Reduce(function(...) rbind(...), list(coef_fl_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Recruitment 
df_re_size <- df %>% subset(recruit == 1)

# Miscellany
coef_misc   <- data.frame(coefficient = c('rec_siz', 'rec_sd',
                                          'fecu_b0', 
                                          'max_siz', 'min_siz'),
                          value       = c(mean(log(df_re_size$size_t0), na.rm = T), 
                                          sd(  log(df_re_size$size_t0), na.rm = T),
                                          repr_pc_median,
                                          df_gr$logsize_t0 %>% max, 
                                          df_gr$logsize_t0 %>% min))

extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = v_script_prefix,
  species = v_species,
  surv_b0 = extr_value(coef_su, 'b0'),
  surv_b1 = extr_value(coef_su, 'logsize_t0'),
  surv_b2 = extr_value(coef_su, 'logsize_t0_2'),
  surv_b3 = extr_value(coef_su, 'logsize_t0_3'),
  grow_b0 = extr_value(coef_gr, 'b0'),
  grow_b1 = extr_value(coef_gr, 'logsize_t0'),
  grow_b2 = extr_value(coef_gr, 'logsize_t0_2'),
  grow_b3 = extr_value(coef_gr, 'logsize_t0_3'),
  a       = extr_value(coef_gr, 'a'),
  b       = extr_value(coef_gr, 'b'),
  fl_b0   = extr_value(coef_fl, 'b0'),
  fl_b1   = extr_value(coef_fl, 'logsize_t0'),
  fl_b2   = extr_value(coef_fl, 'logsize_t0_2'),
  fl_b3   = extr_value(coef_fl, 'logsize_t0_3'),
  fr_b0   = extr_value(coef_fl, 'b0'),
  fr_b1   = extr_value(coef_fr, 'logsize_t0'),
  fr_b2   = extr_value(coef_fr, 'logsize_t0_2'),
  fr_b3   = extr_value(coef_fr, 'logsize_t0_3'),
  fecu_b0 = extr_value(coef_misc, 'fecu_b0'),
  recr_sz = extr_value(coef_misc, 'rec_siz'),
  recr_sd = extr_value(coef_misc, 'rec_sd'),
  L       = extr_value(coef_misc, 'min_siz'),
  U       = extr_value(coef_misc, 'max_siz'),
  mat_siz = 200,
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index
))


# Building the IPM -------------------------------------------------------------
# Function describing the invert logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

# Survival of x-sized individual to time t1
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  for (i in 1:num_pars) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^(i)
    }
  }
  return(inv_logit(survival_value))
}

# Function describing standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- 0
  for (i in 0:num_pars) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}

# Function describing the transition kernel
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Flowering of x-sized individual at time t0
fl_x <- function(x, pars, num_pars = v_mod_fl_index) {
  val <- pars$fl_b0
  for (i in 1:num_pars) {
    param <- paste0('fl_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  inv_logit(val)
}

# Fruiting of x-sized individuals at time t0
fr_x <- function(x, pars, num_pars = v_mod_fr_i) {
  val <- pars$fr_b0
  for (i in 1:num_pars) {
    param <- paste0('fr_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  exp(val)  # Negative binomial uses log link
}

# Recruitment size distribution at time t1
re_y_dist <- function(y, pars) {
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
}

# F-kernel
fyx <- function(y, x, pars) {
  fl_x(x, pars) *
    fr_x(x, pars) *
    pars$fecu_b0 *
    re_y_dist(y, pars)
}

# Kernel
kernel <- function(pars) {
  
  # number of bins over which to integrate
  n   <- pars$mat_siz 
  # lower limit of integration
  L   <- pars$L  
  # upper limit of integration
  U   <- pars$U       
  # bin size
  h   <- (U - L) / n  
  # lower boundaries of bins
  b   <- L + c(0:n) * h             
  # midpoints of bins
  y   <- 0.5 * (b[1:n] + b[2:(n + 1)]) 
  
  # Survival vector
  Smat   <- c()
  Smat   <- sx(y, pars)
  
  # Growth matrix
  Gmat   <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Growth/survival transition matrix
  Tmat   <- matrix(0, n, n)
  
  # Correct for eviction of offspring
  for(i in 1:(n / 2)) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for(i in (n / 2 + 1):n) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Fertility matrix
  Fmat <- outer(y, y, Vectorize(function(x, y) fyx(x, y, pars))) * h
  
  # Full Kernel is simply a summation of fertility and transition matrices
  k_yx <- Fmat + Tmat
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y))
}

lambda_ipm <- function(i) {
  return(Re(eigen(kernel(i)$k_yx)$value[1]))
}

# mean population growth rate
lam_mean <- lambda_ipm(pars)
lam_mean


# Observed population growth ---------------------------------------------------
df_counts_year <- df %>%
  group_by(year) %>%
  filter(!is.na(size_t0)) %>% 
  summarise(n = n())

# Then compute observed lambda
lam_obs_y <- df_counts_year$n[-1] / df_counts_year$n[-nrow(df_counts_year)]
lam_obs_mean <- mean(lam_obs_y, na.rm = TRUE)


# IPM investigation ------------------------------------------------------------
summary(df$size_t0, na.rm = TRUE)
hist(df$size_t0, na.rm = TRUE)

# Use actual observed size limits
min_x <- min(df$size_t0, na.rm = TRUE)
max_x <- max(df$size_t0, na.rm = TRUE)

# Recalculate mesh points
n_mesh <- 200  # resolution
x_vals <- seq(min_x, max_x, length.out = n_mesh)
y_vals <- x_vals  # assuming y follows same range


fl_vals <- sapply(x_vals, function(x) fl_x(x, pars))
plot(x_vals, fl_vals, type = "l", main = "Flowering Probability vs Size",
     xlab = "Size (x)", ylab = "Flowering Probability")


fr_vals <- sapply(x_vals, function(x) fr_x(x, pars))
plot(x_vals, fr_vals, type = "l", main = "Fruiting Counts vs Size",
     xlab = "Size (x)", ylab = "Mean Number of Fruits")


re_vals <- sapply(y_vals, function(y) re_y_dist(y, pars))
plot(y_vals, re_vals, type = "l", main = "Recruitment Size Distribution",
     xlab = "Size (y)", ylab = "Density")


fyx_vals <- outer(x_vals, y_vals, Vectorize(function(x, y) fyx(x, y, pars)))
selected_x <- c(2, 5, 8)  # example sizes
matplot(y_vals, t(fyx_vals[selected_x, ]), type = "l", lty = 1, col = 1:length(selected_x),
        main = "Kernel fyx across y for selected x",
        xlab = "Size (y)", ylab = "fyx value")
legend("topright", legend = paste("x =", selected_x), col = 1:length(selected_x), lty = 1)


# Approximated lambda
dx <- x_vals[2] - x_vals[1]
dy <- y_vals[2] - y_vals[1]
lambda_approx <- sum(fyx_vals) * dx * dy

cat("F-only lambda (approx):", lambda_approx, "\n")
cat("Full IPM lambda (eigen):", lam_mean, "\n")


# Save data --------------------------------------------------------------------
write.csv(df_og, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original.csv')))
write.csv(df_meta, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_meta.csv')))
write.csv(df, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata.csv')))
write.csv(df_su, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_survival.csv')))
write.csv(df_gr, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_growth.csv')))
write.csv(df_fl, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_flower.csv')))
write.csv(df_fr, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_fruit.csv')))
write.csv(df_re, row.names = F,
          file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_recruit.csv')))
