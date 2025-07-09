# IPM mean with fire - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.07.09 (y-m-d)


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
    logsize_t0_3 = logsize_t0^3,
    # FIRE - if there was a fire in this year then 1, else 0
    fire         = case_when(
      fire_sev >= 0 ~ 1,
      (is.na(fire_sev) & survives >= 0) ~ 0, 
      TRUE ~ NA_real_))


# Survival data ----------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  dplyr::select(plant_id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, 
                fire)

df_su_f0 <- df_su %>% filter(fire == 0)
df_su_f1 <- df_su %>% filter(fire == 1)

fig_su_overall_f0 <- ggplot(
  data = plot_binned_prop(df_su_f0, 10, logsize_t0, survives)) +
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
       y = expression('Survival without fire event'))

fig_su_overall_f1 <- ggplot(
  data = plot_binned_prop(df_su_f1, 10, logsize_t0, survives)) +
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
  labs(title = 'Survival with fire event',
       subtitle = v_ggp_suffix,
       x = expression('log(size)'[t0]),
       y = expression(''))

fig_su_overall <- fig_su_overall_f0 + fig_su_overall_f1 + plot_layout()


# Survival model ---------------------------------------------------------------
# Logistic regression
mod_su_0 <- glm(survives ~ fire,
                data = df_su, family = 'binomial') 
# Logistic regression
mod_su_1 <- glm(survives ~ logsize_t0 + fire,
                data = df_su, family = 'binomial') 
# Quadratic logistic model
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2 + fire,
                data = df_su, family = 'binomial')  
# Cubic logistic model
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire,
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

mod_su_bestfit <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef   <- coef(mod_su_bestfit)

# Create prediction data frames for fire = 0 and fire = 1
df_su_pred_0 <- data.frame(logsize_t0 = seq(min(df_su$logsize_t0), max(df_su$logsize_t0), length.out = 100), fire = 0)
df_su_pred_1 <- data.frame(logsize_t0 = seq(min(df_su$logsize_t0), max(df_su$logsize_t0), length.out = 100), fire = 1)

# Add predictions using model
df_su_pred_0$survives <- predict(mod_su_bestfit, newdata = df_su_pred_0, type = "response")
df_su_pred_1$survives <- predict(mod_su_bestfit, newdata = df_su_pred_1, type = "response")

fig_su_line_0 <- ggplot() +
  geom_jitter(data = df_su_f0, aes(x = logsize_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred_0, aes(x = logsize_t0, y = survives),
            color = "red", lwd = 2) +  
  theme_bw() + 
  labs(title = 'Survival (No Fire)',
       y = 'Survival to time t1', x = '')

fig_su_bin_0 <- ggplot() +
  geom_point(data = plot_binned_prop(df_su_f0, 10, logsize_t0, survives),
             aes(x = logsize_t0, y = survives)) +
  geom_errorbar(data = plot_binned_prop(df_su_f0, 10, logsize_t0, survives),
                aes(x = logsize_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_su_pred_0, aes(x = logsize_t0, y = survives),
            color = 'red', lwd = 2) +
  theme_bw() +
  ylim(0, 1) + 
  labs(y = 'Survival to time t1', x = 'Size at time t0 (log())')

fig_su_line_1 <- ggplot() +
  geom_jitter(data = df_su_f1, aes(x = logsize_t0, y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred_1, aes(x = logsize_t0, y = survives),
            color = "red", lwd = 2) +  
  theme_bw() + 
  labs(title = 'Survival (Fire)',
       y = '', x = '')

fig_su_bin_1 <- ggplot() +
  geom_point(data = plot_binned_prop(df_su_f1, 10, logsize_t0, survives),
             aes(x = logsize_t0, y = survives)) +
  geom_errorbar(data = plot_binned_prop(df_su_f1, 10, logsize_t0, survives),
                aes(x = logsize_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_su_pred_1, aes(x = logsize_t0, y = survives),
            color = 'red', lwd = 2) +
  theme_bw() +
  ylim(0, 1) +
  labs(y = '', x = 'Size at time t0 (log())')

fig_su_all <- (fig_su_line_0 + fig_su_line_1) /
  (fig_su_bin_0 + fig_su_bin_1)  + 
  plot_annotation(title = "Survival",
                  subtitle = v_ggp_suffix,
                  theme = theme(
                    plot.title =    element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10, face = "italic")))

fig_su_all


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  dplyr::select(plant_id, year, size_t0, size_t1,
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, 
                fire)

df_gr_f0 <- df_gr %>% filter(fire == 0)
df_gr_f1 <- df_gr %>% filter(fire == 1)

fig_gr_overall_f0 <- ggplot(
  data  = df_gr_f0, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_overall_f1 <- ggplot(
  data  = df_gr_f1, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_overall <- fig_gr_overall_f0 + fig_gr_overall_f1 + plot_layout()
fig_gr_overall


# Growth model -----------------------------------------------------------------
# Intercept model 
mod_gr_0 <- lm(logsize_t1 ~ fire,
               data = df_gr)
# Linear model
mod_gr_1 <- lm(logsize_t1 ~ logsize_t0 + fire, 
               data = df_gr)
# Quadratic model
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + fire, 
               data = df_gr)  
# Cubic model
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire, 
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

# Create prediction data for fire = 0 and fire = 1
df_gr_pred_0 <- data.frame(logsize_t0 = seq(min(df_gr$logsize_t0, na.rm = TRUE),
                                            max(df_gr$logsize_t0, na.rm = TRUE), 
                                            length.out = 100)) %>%
  mutate(logsize_t0_2 = logsize_t0^2,
         fire = 0)

df_gr_pred_1 <- data.frame(logsize_t0 = seq(min(df_gr$logsize_t0, na.rm = TRUE),
                                            max(df_gr$logsize_t0, na.rm = TRUE), 
                                            length.out = 100)) %>%
  mutate(logsize_t0_2 = logsize_t0^2,
         fire = 1)

# Predict from model
df_gr_pred_0$logsize_t1 <- predict(mod_gr_bestfit, newdata = df_gr_pred_0)
df_gr_pred_1$logsize_t1 <- predict(mod_gr_bestfit, newdata = df_gr_pred_1)

# Plot lines
fig_gr_line_0 <- ggplot() +
  geom_point(data = df_gr_f0, aes(x = logsize_t0, y = logsize_t1), alpha = 0.4) +
  geom_line(data = df_gr_pred_0, aes(x = logsize_t0, y = logsize_t1), color = 'green', lwd = 2) +
  theme_bw() +
  labs(title = 'Growth (No Fire)', y = 'Size at time t1 (log)', 
       x = 'Size at time t0 (log)')

fig_gr_line_1 <- ggplot() +
  geom_point(data = df_gr_f1, aes(x = logsize_t0, y = logsize_t1), alpha = 0.4) +
  geom_line(data = df_gr_pred_1, aes(x = logsize_t0, y = logsize_t1), color = 'green', lwd = 2) +
  theme_bw() +
  labs(title = 'Growth (Fire)', y = '', x = 'Size at time t0 (log)')

# Plot predicted vs observed (optional but recommended)
fig_gr_pred_0 <- ggplot(df_gr_f0, aes(x = predict(mod_gr_bestfit, newdata = df_gr_f0), y = logsize_t1)) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, color = 'red', lwd = 1.5) +
  theme_bw() +
  labs(x = 'Predicted', y = 'Observed')

fig_gr_pred_1 <- ggplot(df_gr_f1, aes(x = predict(mod_gr_bestfit, newdata = df_gr_f1), y = logsize_t1)) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, color = 'red', lwd = 1.5) +
  theme_bw() +
  labs(x = 'Predicted', y = '')

# Combine
fig_gr_all <- (fig_gr_line_0 + fig_gr_line_1) /
  (fig_gr_pred_0 + fig_gr_pred_1) +
  plot_annotation(title = "Growth Prediction",
                  subtitle = v_ggp_suffix,
                  theme = theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10, face = "italic"))) &
  theme(axis.title.y = element_text(size = 12))

fig_gr_all


# Growth variance --------------------------------------------------------------
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)  
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2  
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE)) 


# Flowering data ----------------------------------------------------------------
df_fl <- df %>% 
  filter(!is.na(flower)) %>%
  dplyr::select(plant_id, year, size_t0, flower, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, 
                fire) %>% 
  mutate(flower = if_else(flower > 0, 1, flower))

df_fl_f0 <- df_fl %>% filter(fire == 0)
df_fl_f1 <- df_fl %>% filter(fire == 1)

fig_fl_overall_f0 <- ggplot(
  data = plot_binned_prop(df_fl_f0, 10, logsize_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logsize_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Flowering without fire event',
       subtitle = v_ggp_suffix,
       x = expression('log(size)'[t0]),
       y = expression('Flowering probability'))

fig_fl_overall_f1 <- ggplot(
  data = plot_binned_prop(df_fl_f1, 10, logsize_t0, flower)) +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower), 
              position = position_jitter(width = 0.1, height = 0.3)
              , alpha = .1) +
  geom_point(aes(x = logsize_t0, y = flower),
             alpha = 1, pch = 16, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                linewidth = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title = 'Flowering with fire event',
       subtitle = v_ggp_suffix,
       x = expression('log(size)'[t0]),
       y = expression(''))

fig_fl_overall <- fig_fl_overall_f0 + fig_fl_overall_f1 + plot_layout()
fig_fl_overall


# Flower model -----------------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ fire,
                data = df_fl, family = 'binomial') 
# Logistic regression
mod_fl_1 <- glm(flower ~ logsize_t0 + fire,
                data = df_fl, family = 'binomial') 
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logsize_t0 + logsize_t0_2 + fire,
                data = df_fl, family = 'binomial')  
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire,
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

# Create prediction data frames for fire = 0 and fire = 1
df_fl_pred_0 <- data.frame(logsize_t0 = seq(min(df_fl$logsize_t0, na.rm = T), 
                                            max(df_fl$logsize_t0, na.rm = T), 
                                            length.out = 100), 
                           fire = 0) %>% mutate(logsize_t0_2 = logsize_t0^2)
df_fl_pred_1 <- data.frame(logsize_t0 = seq(min(df_fl$logsize_t0, na.rm = T), 
                                            max(df_fl$logsize_t0, na.rm = T), 
                                            length.out = 100), 
                           fire = 1) %>% mutate(logsize_t0_2 = logsize_t0^2)

# Add predictions using model
df_fl_pred_0$flower <- predict(mod_fl_bestfit, newdata = df_fl_pred_0, type = "response")
df_fl_pred_1$flower <- predict(mod_fl_bestfit, newdata = df_fl_pred_1, type = "response")

fig_fl_line_0 <- ggplot() +
  geom_jitter(data = df_fl_f0, aes(x = logsize_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred_0, aes(x = logsize_t0, y = flower),
            color = "red", lwd = 2) +  
  theme_bw() + 
  labs(title = 'Flowering (No Fire)',
       y = 'Flowering at time t0', x = '')

fig_fl_bin_0 <- ggplot() +
  geom_point(data = plot_binned_prop(df_fl_f0, 10, logsize_t0, flower),
             aes(x = logsize_t0, y = flower)) +
  geom_errorbar(data = plot_binned_prop(df_fl_f0, 10, logsize_t0, flower),
                aes(x = logsize_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_fl_pred_0, aes(x = logsize_t0, y = flower),
            color = 'red', lwd = 2) +
  theme_bw() +
  ylim(0, 1) + 
  labs(y = 'Flowering at time t0', x = 'Size at time t0 (log())')

fig_fl_line_1 <- ggplot() +
  geom_jitter(data = df_fl_f1, aes(x = logsize_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred_1, aes(x = logsize_t0, y = flower),
            color = "red", lwd = 2) +  
  theme_bw() + 
  labs(title = 'Flowering (Fire)',
       y = '', x = '')

fig_fl_bin_1 <- ggplot() +
  geom_point(data = plot_binned_prop(df_fl_f1, 10, logsize_t0, flower),
             aes(x = logsize_t0, y = flower)) +
  geom_errorbar(data = plot_binned_prop(df_fl_f1, 10, logsize_t0, flower),
                aes(x = logsize_t0, ymin = lwr, ymax = upr)) +
  geom_line(data = df_fl_pred_1, aes(x = logsize_t0, y = flower),
            color = 'red', lwd = 2) +
  theme_bw() +
  ylim(0, 1) +
  labs(y = '', x = 'Size at time t0 (log())')

fig_fl_all <- (fig_fl_line_0 + fig_fl_line_1) /
  (fig_fl_bin_0 + fig_fl_bin_1)  + 
  plot_annotation(title = "Flowering",
                  subtitle = v_ggp_suffix,
                  theme = theme(
                    plot.title =    element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10, face = "italic")))

fig_fl_all


# # Fruit data -------------------------------------------------------------------
# df_fr <- df %>%
#   filter(!is.na(fruit), !is.na(logsize_t0))
# 
# # Over dispersed -> negative binomial
# # otherwise Poisson
# mean(df_fr$fruit)
# var(df_fr$fruit)
# 
# 
# # Fruit model ------------------------------------------------------------------
# fr_mod_0 <- glm.nb(fruit ~ fire, data = df_fr)
# fr_mod_1 <- glm.nb(fruit ~ logsize_t0 + fire, data = df_fr)
# fr_mod_2 <- glm.nb(fruit ~ logsize_t0 + logsize_t0_2 + fire, data = df_fr)
# fr_mod_3 <- glm.nb(fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire, data = df_fr)
# 
# fr_mods      <- list(fr_mod_0, fr_mod_1, fr_mod_2, fr_mod_3)
# fr_mods_dAIC <- AICtab(fr_mods, weights = T, sort = F)$dAIC
# 
# # Get the sorted indices of dAIC values
# fr_mods_i_sort <- order(fr_mods_dAIC)
# 
# # Establish the index of model complexity
# if (length(v_mod_set_fr) == 0) {
#   fr_mod_i_best <- fr_mods_i_sort[1]
#   v_mod_fr_i    <- fr_mod_i_best - 1 
# } else {
#   fr_mod_i_best <- v_mod_set_fr +1
#   v_mod_fr_i    <- v_mod_set_fr
# }
# 
# fr_mod_best  <- fr_mods[[fr_mod_i_best]]
# fr_mod_ranef <- coef(fr_mod_best)
# 
# # Prediction data for fire = 0
# df_pred_0 <- data.frame(
#   logsize_t0   = seq(min(df_fr$logsize_t0, na.rm = TRUE),
#                      max(df_fr$logsize_t0, na.rm = TRUE),
#                      length.out = 200),
#   fire         = 0
# )
# df_pred_0$logsize_t0_2 <- df_pred_0$logsize_t0^2
# df_pred_0$pred         <- predict(fr_mod_best, newdata = df_pred_0, type = "response")
# 
# # Prediction data for fire = 1
# df_pred_1 <- data.frame(
#   logsize_t0   = df_pred_0$logsize_t0,  # same x values
#   fire         = 1
# )
# df_pred_1$logsize_t0_2 <- df_pred_1$logsize_t0^2
# df_pred_1$pred         <- predict(fr_mod_best, newdata = df_pred_1, type = "response")
# 
# # Plot
# fig_fr <- ggplot(df_fr, aes(x = logsize_t0, y = fruit)) +
#   geom_point(alpha = 0.3, size = 1) +
#   geom_line(data = df_pred_0, aes(x = logsize_t0, y = pred), 
#             color = "#1b9e77", linewidth = 1.2, linetype = "solid") +
#   geom_line(data = df_pred_1, aes(x = logsize_t0, y = pred), 
#             color = "#d95f02", linewidth = 1.2, linetype = "dashed") +
#   theme_minimal() +
#   labs(
#     title    = 'Fruit Production by log(Size) and Fire',
#     subtitle = v_ggp_suffix,
#     x        = 'log(Size at t0)',
#     y        = 'Number of Fruits'
#   ) +
#   annotate("text", x = Inf, y = Inf, label = "Solid = no fire; Dashed = fire",
#            hjust = 1.1, vjust = 1.5, size = 3, color = "gray40")
# 
# fig_fr
