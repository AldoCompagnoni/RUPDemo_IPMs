# IPM mean - Archbold -  - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.04.02


# Website    : 
# Publication: 


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
  lubridate,
  ipmr,
  janitor) # , skimr, binom, lme4


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Chrysopsis highlandsensis')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()


# Create a unique species abbreviation for file naming
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Define script prefix
v_script_prefix <- str_c(v_head)

v_suffix <- ''

# Plot subtitle
v_ggp_suffix    <- paste(
  tools::toTitleCase(v_head), '-', v_species)


# Models
v_mod_set_gr <- c()
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
df_og <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    plant_id = identifier,
    year     = year0,
    survival = survival_1) %>%  
  mutate(
    plant_id = as.factor(plant_id)) %>%
  arrange(site, plant_id, year, survival) # quad, quad_id, plant, , month

# df_meta <- data.frame(variable = colnames(df_og)) %>% 
#   mutate(definition = c(
#     'study site', 'quadrat number',	'macroplot number',	
#     'plant number (within quad)', 'direction within circular quad',	
#     'distance from quad center',	'was quad caged from 2012 onward',	
#     'vegetation type', 'year quadrat was initiated', 
#     'year-month of observation',	'fire severity Dec. 2014 or Jan. 2015',
#     'fire severity Aug. 2005',	'fire severity May/June 2009',
#     "fire severity B's ridge 2016", 'fire severity Feb. 2017',
#     'fire severity Oct. 2017', 'survival code for month', 'number of stems',
#     'number of branch tips', 'number of flowers (corolla showing)', 
#     'number of developing fruits', 'number of mature fruits', 'herbivory code',
#     'plant identification', 'quadrat identification', 'sample year', 
#     'sample month'))

df_og_extended <- df_og %>%
    # Recruits
  mutate(
    recruit = ifelse(astg == 1, 1, 0))


# Fire data
df_fire <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_fire.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    year = burn_yr,
    fire = treatment) %>% 
  mutate(fire = as.factor(fire)) %>%  
  select(!c(year0, notes))

# Mean data frame --------------------------------------------------------------
df <- df_og_extended %>%
  group_by(site , plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(survival )), NA_real_, min(survival,  na.rm = T)),
    size_t0  = if_else(all(is.na(dia)),       NA_real_, max(dia,        na.rm = T)), # height of tallest scape (cm)
    size_t1  = if_else(all(is.na(dia_1)),     NA_real_, max(dia_1,      na.rm = T)),
    flower   = if_else(all(is.na(hd)),        NA_real_, max(hd,        na.rm = T)),  
    recruit  = if_else(all(is.na(recruit)),   NA_real_, min(recruit,   na.rm = T)),
  ) %>% 
  ungroup() %>% 
  mutate(
    logsize_t0   = log(size_t0), # avoiding size_t0 = 0    
    logsize_t1   = log(size_t1),    
    logsize_t0_2 = logsize_t0^2,     
    logsize_t0_3 = logsize_t0^3) %>%
  full_join(df_fire, by = c("site", "year")) %>%
  mutate(
    fire = case_when(
      is.na(fire)      ~ 0,
      fire == "burn"   ~ 1,
      TRUE             ~ NA)) # since there is also mechanical return we have them as NA
  


# Prepare data frames for analysis ---------------------------------------------
# Survival data frame
surv_df <- subset(df, !is.na(survives)) %>%
  subset(size_t0 != 0) %>%
  select(site, plant_id, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         fire) %>% 
  filter(!is.na(fire)) %>%
  mutate(
    fire = factor(fire,
                  levels = c(0, 1),
                  labels = c("No fire", "Fire")))


# Growth data frame
grow_df <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(site, plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Recruit data frame
# Parent area cover data at year t0
# Total area data frame
site_df <- df %>%  
  group_by (site, year) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

group_df <- site_df %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

cover_df <- left_join(site_df, group_df) %>%
  mutate(year = year + 1) %>%
  mutate(year = as.integer(year)) %>% 
  drop_na()

# Recruitment data at year t1
recr_df <- df %>% 
  group_by (year, site) %>%
  summarise(nr_rec = sum(recruit, na.rm = T)) %>% 
  ungroup

recr_df <- full_join(cover_df, recr_df, by = c('site', 'year'))


# Plotting the data --------------------------------------------------------
# Survival analysis
# Load custom function for plotting binned proportions
source('helper_functions/plot_binned_prop.R')

# Generate a plot for overall survival based on size at t0
df_su_binned <- surv_df %>%
  group_split(fire) %>%
  purrr::map_df(~ plot_binned_prop(.x, 10, logsize_t0, survives) %>%
                  mutate(fire = unique(.x$fire)))

fig_su_overall <- ggplot(df_su_binned, aes(x = logsize_t0, y = survives, color = fire)) +
  geom_jitter(data = surv_df, aes(x = logsize_t0, y = survives, color = fire),
              position = position_jitter(width = 0.1, height = 0.3), alpha = 0.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "top") +
  labs(
    title = "Survival",
    subtitle = v_ggp_suffix,
    x = expression('log(size)'[t0]),
    y = "Survival Probability")

fig_su_overall

# Growth analysis
# Create a scatter plot of size at t0 versus size at t1
g_gr_overall <- ggplot(
  data  = grow_df, aes( x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))

ggsave(paste0(dir_result, '/1.2_overall_gr', v_suffix,'.png'), 
       plot = g_gr_overall, 
       width = 4, height = 3, units = 'in', dpi = 150)

# Recruitment analysis
# Create a scatter plot showing the relationship between 
# total parent plant area and number of recruits
g_rec_overall <- ggplot(
  recr_df, aes(x = tot_p_area, y = nr_rec)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))

# Save the recruitment plot as a PNG file
ggsave(paste0(dir_result, '/1.3_overall_rec', v_suffix,'.png'), 
       plot = g_rec_overall, 
       width = 4, height = 3, units = 'in', dpi = 150)


# Fit vital rate models for the mean IPM ---------------------------------------
# Survival
# Fit models to predict survival based on size at time t0
# Logistic regression
su_mod_mean_0 <- glm(survives ~ 1, 
                     data = surv_df, family = 'binomial') 
# Logistic regression
su_mod_mean   <- glm(survives ~ logsize_t0, 
                     data = surv_df, family = 'binomial') 
# Quadratic logistic model
su_mod_mean_2 <- glm(survives ~ logsize_t0 + logsize_t0_2, 
                     data = surv_df, family = 'binomial')  
# Cubic logistic model
su_mod_mean_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
                     data = surv_df, family = 'binomial')  

# Compare models using AIC
su_mods              <- list(su_mod_mean_0, su_mod_mean, su_mod_mean_2, su_mod_mean_3)
su_dAIC_values       <- bbmle::AICctab(su_mods, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
su_sorted_indices    <- order(su_dAIC_values)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- su_sorted_indices[1]
  v_mod_su_index <- mod_su_index_bestfit - 1 
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index <- v_mod_set_su
}

su_mod_bestfit   <- su_mods[[mod_su_index_bestfit]]
su_ranef         <- coef(su_mod_bestfit)

# Generate predictions for survival across a range of sizes
surv_x <- seq(min(surv_df$logsize_t0, na.rm = T), 
              max(surv_df$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
surv_pred_df <- predictor_fun(surv_x, su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = surv_x, survives = .)

# Plot observed survival with fitted line
g_surv_line <- ggplot() +
  geom_jitter(data = surv_df, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0, height = 0.25) + 
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = line_color_pred_fun(su_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

# Plot binned survival proportions with error bars
g_surv_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        y = survives) ) +
  geom_errorbar(data =  plot_binned_prop(
    df, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        ymin = lwr,
        ymax = upr) ) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
g_surv_overall_pred <- g_surv_line + g_surv_bin + plot_layout()
g_surv_overall_pred

# Save the survival prediction plot
ggsave(paste0(dir_result, '/2.2_overall_surv_pred_logs', v_suffix, '.png'),
       plot = g_surv_overall_pred,
       width = 8, height = 3, units = 'in', dpi = 150)


# Growth
# Fit growth models to predict size at time t1 based on size at time t0
# Intercept model
gr_mod_mean_0 <- lm(logsize_t1 ~ 
                      1, data = grow_df)
# Linear model
gr_mod_mean   <- lm(logsize_t1 ~ 
                      logsize_t0, data = grow_df)
# Quadratic model
gr_mod_mean_2 <- lm(logsize_t1 ~ 
                      logsize_t0 + logsize_t0_2, data = grow_df)  
# Cubic model
gr_mod_mean_3 <- lm(logsize_t1 ~ 
                      logsize_t0 + logsize_t0_2 + logsize_t0_3, data = grow_df)

# Compare models using AIC
gr_mods              <- list(gr_mod_mean_0, gr_mod_mean, gr_mod_mean_2, gr_mod_mean_3)
gr_dAIC_values       <- bbmle::AICctab(gr_mods, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
gr_sorted_indices    <- order(gr_dAIC_values)

# Establish the index of model complexity
if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- gr_sorted_indices[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr
}

gr_mod_bestfit       <- gr_mods[[mod_gr_index_bestfit]]
gr_ranef             <- coef(gr_mod_bestfit)

# Predict size at time t1 using the mean growth model
grow_df$pred <- predict(gr_mod_bestfit, type = 'response')

# Plot observed size at time t1 against size at time t0 with the fitted line
# Define functions
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

g_grow_line <- ggplot(
  grow_df, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, gr_ranef), 
                color = line_color_pred_fun(gr_ranef), 
                lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

# Plot predicted versus observed size at time t1
g_grow_pred <- ggplot(
  grow_df, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw()

# Combine growth line and prediction plots
g_grow_overall_pred <- g_grow_line + g_grow_pred + plot_layout() 
g_grow_overall_pred

# Save the growth prediction plot
ggsave(paste0(dir_result, '/2.1_overall_grow_pred_logs', v_suffix, '.png'),
       plot = g_grow_overall_pred,
       width = 8, height = 4, units = 'in', dpi = 150)

# Fit a model to assess variance in growth
# Fitted values from growth model
x         <- fitted(gr_mod_bestfit)  
# Squared residuals
y         <- resid(gr_mod_bestfit)^2  
# Non-linear model for variance
gr_var_m  <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0),
                 control = nls.control(
                   maxiter = 1000, tol = 1e-6, warnOnly = TRUE))


# Recruitment
# Filter recruitment data to exclude NAs
recr_nona_nr_rec <- recr_df %>% filter(!is.na(nr_rec))
# Fit a negative binomial model for recruitment
rec_mod_mean <- MASS::glm.nb(nr_rec ~ 1, data = recr_nona_nr_rec)

# Generate predictions for recruitment
recr_nona_nr_rec <- recr_nona_nr_rec %>% 
  mutate(pred_mod_mean = predict(rec_mod_mean, type = 'response')) 

# Summarize total number of recruits and predictions
rec_sums_df_m <- recr_nona_nr_rec %>%
  summarize(nr_rec = sum(nr_rec),
            pred_mod_mean = sum(pred_mod_mean))

# Count number of adult individuals
indiv_m <- surv_df %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
repr_pc_m <- indiv_m %>%
  bind_cols(rec_sums_df_m) %>%
  mutate(repr_pc_mean = pred_mod_mean / n_adults) %>%
  mutate(repr_pc_obs = nr_rec / n_adults) %>%
  drop_na 


# Exporting parameter estimates ------------------------------------------------
# Growth
grow_fe <- data.frame(coefficient = names(coef(gr_mod_bestfit)),
                      value       = coef(gr_mod_bestfit))
var_m   <- data.frame(coefficient = names(coef(gr_var_m)),
                      value       = coef(gr_var_m))

grow_out <- Reduce(function(...) rbind(...), list(grow_fe, var_m)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

write.csv(grow_out, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_', v_sp_abb, '_grow_pars_mean.csv'))


# Survival
surv_fe <- data.frame(coefficient = names(coef(su_mod_bestfit)),
                      value       = coef(su_mod_bestfit))

surv_out<- Reduce(function(...) rbind(...), list(surv_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

write.csv(surv_out, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_', v_sp_abb, '_surv_pars_mean.csv'))

# Recruitment 
rec_size <- df %>% subset(recruit == 1)

others   <- data.frame(coefficient = c('rec_siz', 'rec_sd', 
                                       'max_siz', 'min_siz',
                                       'fecu_b0'),
                       value       = c(mean(log(rec_size$size_t0), na.rm = TRUE), 
                                       sd(  log(rec_size$size_t0), na.rm = TRUE),
                                       grow_df$logsize_t0 %>% max, 
                                       grow_df$logsize_t0 %>% min,
                                       repr_pc_m$repr_pc_mean))

write.csv(others, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_', v_sp_abb, '_other_pars_mean.csv'))


# Building the IPM from scratch ------------------------------------------------
extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = v_script_prefix,
  species = v_species,
  surv_b0 = extr_value(surv_out, 'b0'),
  surv_b1 = extr_value(surv_out, 'logsize_t0'),
  surv_b2 = extr_value(surv_out, 'logsize_t0_2'),
  surv_b3 = extr_value(surv_out, 'logsize_t0_3'),
  grow_b0 = extr_value(grow_out, 'b0'),
  grow_b1 = extr_value(grow_out, 'logsize_t0'),
  grow_b2 = extr_value(grow_out, 'logsize_t0_2'),
  grow_b3 = extr_value(grow_out, 'logsize_t0_3'),
  a       = extr_value(grow_out, 'a'),
  b       = extr_value(grow_out, 'b'),
  fecu_b0 = extr_value(others, 'fecu_b0'),
  recr_sz = extr_value(others, 'rec_siz'),
  recr_sd = extr_value(others, 'rec_sd'),
  L       = extr_value(others, 'min_siz'),
  U       = extr_value(others, 'max_siz'),
  mat_siz = 200,
  mod_gr_index = v_mod_gr_index,
  mod_su_index = v_mod_su_index
))

write.csv(pars, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_', v_sp_abb, '_pars.csv'))

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

# Function describing the transition kernel
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Function describing the recruitment 
fy <- function(y, pars, h){
  n_recr  <- pars$fecu_b0
  recr_y  <- dnorm(y, pars$recr_sz, max(h/10, pars$recr_sd)) * h
  recr_y  <- recr_y / sum(recr_y)
  f       <- n_recr * recr_y
  return(f)
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
  
  # Fertility matrix
  Fmat        <- matrix(0, n, n)
  Fmat[]      <- matrix(fy(y, pars, h), n, n)
  
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

# observed population growth rate
# Population counts at time t0
pop_counts_t0 <- df %>%
  group_by(year, site) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  group_by(year, site) %>%
  summarize(n_t1 = n()) %>% 
  ungroup 

# Calculate observed population growth rates, 
# accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
  # by dropping NAs, we remove gaps in sampling!
  drop_na %>% 
  group_by(year) %>% 
  summarise(n_t0 = sum(n_t0),
            n_t1 = sum(n_t1)) %>% 
  ungroup %>% 
  mutate(obs_pgr = n_t1 / n_t0)

# Geometric mean of yearly population growth rates
lam_mean_count <- exp(mean(log(pop_counts$obs_pgr), na.rm = T))

# Overall (aggregated) population growth rate
lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)


# Building the IPM with ipmr ---------------------------------------------------
proto_ipm_p <- init_ipm(sim_gen   = 'simple',
                        di_dd     = 'di',
                        det_stoch = 'det') %>% 
  define_kernel(
    name      = 'P',
    family    = 'CC',
    formula   = s * g,
    s         = plogis(
      surv_b0 + 
        (if (mod_su_index >= 1) surv_b1 * size_1   else 0) +
        (if (mod_su_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3 * size_1^3 else 0)),
    
    mu_g      = grow_b0 + 
      (if (mod_gr_index >= 1) grow_b1 * size_1   else 0) +
      (if (mod_gr_index >= 2) grow_b2 * size_1^2 else 0) +
      (if (mod_gr_index >= 3) grow_b3 * size_1^3 else 0),
    
    g         = dnorm(size_2, mu_g, grow_sig),
    grow_sig  = sqrt(a * exp(b * size_1)),
    data_list = pars,
    states    = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun = 'norm', target = 'g')
  ) %>% 
  
  define_kernel(
    name      = 'F',
    family    = 'CC',
    formula   = fecu_b0 * r_d,
    r_d       = dnorm(size_2, recr_sz, recr_sd),
    data_list = pars,
    states    = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'r_d')
  ) %>% 
  
  define_impl(
    make_impl_args_list(
      kernel_names = c(  'P', 'F'),
      int_rule     = rep('midpoint', 2),
      state_start  = rep('size', 2),
      state_end    = rep('size', 2))
  ) %>%
  
  define_domains(
    size = c(pars$L,
             pars$U,
             pars$mat_siz)
  ) %>%
  
  define_pop_state(
    n_size = rep(1 / 200, 200)
  )

ipmr_p <- make_ipm(proto_ipm  = proto_ipm_p, iterations = 200)

lam_mean_ipmr <- lambda(ipmr_p)

lam_out       <- data.frame(coefficient = names(lam_mean_ipmr), 
                            value       = lam_mean_ipmr)
write.csv(lam_out, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_',v_sp_abb, '_lambda_vec.csv'))

lam_out_wide  <- as.list(pivot_wider(lam_out, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value'))

write.csv(lam_out_wide, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_',v_sp_abb, '_lambda.csv'))
