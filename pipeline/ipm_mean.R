# IPM mean to padrino pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date: 2024.11.08

# reading in, and cleaning the data
#  exploring the overall-years rates
#  setting up the vital rate data-frames for the year specific 

# Setting the stage ------------------------------------------------------------

# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source( 'helper_functions/load_packages.R' )
load_packages( tidyverse, patchwork, skimr, ipmr, binom, bbmle )


# Data -------------------------------------------------------------------------
# Create a unique species abbreviation for file naming
sp_abb  <- tolower(gsub(' ', '', 
                        paste(substr(unlist(strsplit(species, ' ')), 1, 2),
                              collapse = '')))

# Directory 
pub_dir    <- file.path(paste0(author_year, '_', region_abb))
R_dir      <- file.path(pub_dir, 'R',       sp_abb)
data_dir   <- file.path(pub_dir, 'data',    sp_abb)
result_dir <- file.path(pub_dir, 'results', sp_abb)

# Define script prefix
script_prefix <- str_c( str_extract(author_year, "^[^_]+"), 
                        str_sub(str_extract(author_year, "_\\d+$"), -2, -1) )

# anderson16 et al requires state identifier
if( grepl('anderson16',script_prefix) ){
  script_prefix <- str_c(script_prefix, region_abb )
}

# Plant tracker if its not already exists
if (!file.exists(paste0(data_dir, '/', script_prefix, '_', sp_abb, '.csv'))) {
  source(paste0(R_dir, '/', script_prefix, '_', sp_abb, '_tracker.R'))
}

## Read and clean the species data
df <- read.csv(paste0(data_dir, '/', script_prefix, '_', sp_abb, '.csv')) %>% 
  filter(Species == species) %>%
  select(-c(Suspect, nearEdge, Site)) %>%
  mutate(across(c(Quad), as.factor)) %>%
  rename(species  = Species, 
         size_t0  = basalArea_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1,
         year     = Year,
         quad     = Quad,
         track_id = trackID) %>%
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),
         logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3)


# Data exploration -------------------------------------------------------------
# Analyze the quadrat inventory by year
# Group data by quadrant and year, counting the number of individuals
inv_plot_per_year <- df %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>%
  # Create a scatter plot of quadrat counts over the years
  ggplot() + 
  geom_point(aes(x = year, y = quad)) +
  theme_bw() +
  labs(title    = 'Sampling inventory',
       subtitle = paste(script_prefix, species)) +
  theme(axis.text.y = element_text(size = 5))
inv_plot_per_year

if (!dir.exists(paste0(pub_dir, '/results'))) {
  dir.create(paste0(pub_dir, '/results'))
}
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

ggsave(paste0(result_dir, '/0.0_quad_per_year.png'), 
       plot = inv_plot_per_year, 
       width = 8, height = 3, units = 'in', dpi = 150)


# Prepare data frames for analysis ---------------------------------------------
# Survival data frame
surv_df <- subset(df, !is.na(survives)) %>%
  subset(size_t0 != 0) %>%
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Growth data frame
grow_df <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Total area data frame
quad_df <- df %>%  
  group_by (quad, year) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

group_df <- quad_df %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

cover_df <- left_join(quad_df, group_df) %>%
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

# Recruitment data frame
recr_df <- df %>%
  group_by (year, quad) %>% 
  summarise(nr_quad = sum(recruit, na.rm = T)) %>% 
  ungroup

recr_df <- left_join(cover_df, recr_df)

# Save data
write.csv(df, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_data_df.csv'))
write.csv(surv_df,
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_survival_df.csv'))
write.csv(grow_df,
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_growth_df.csv'))
write.csv(recr_df,
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_recruitment_df.csv'))


# Plotting the data --------------------------------------------------------
# Create histograms for log-transformed sizes at t0 and t1
hist_t0 <- 
  ggplot(df, aes(x = logsize_t0)) +
  geom_histogram(binwidth = 0.2, fill = 'grey', color = 'black') +
  labs(title    = 'Histogram',
       subtitle = paste(script_prefix, species), 
       x        = 'Size at time t0') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
hist_t1 <- 
  ggplot(df, aes(x = logsize_t1)) +
  geom_histogram(binwidth = 0.2, fill = 'white', color = 'black') +
  labs(x = 'Size at time t1') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

hist_sizes_log <- hist_t0 + hist_t1

ggsave(paste0(result_dir, '/1.0_overall_hist_sizes_log.png'), 
       plot = hist_sizes_log, 
       width = 8, height = 3, units = 'in', dpi = 150)

# Survival analysis
# Load custom function for plotting binned proportions
source('helper_functions/plot_binned_prop.R')

# Generate a plot for overall survival based on size at t0
surv_overall <- ggplot(data = plot_binned_prop(df, 10, logsize_t0, survives)) +
  geom_point(aes(x = logsize_t0, 
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Sampling inventory',
       subtitle = paste(script_prefix, species),
       x        = expression('log(size)'[t0]),
       y        = expression('Survival to time t1'))

# Save the survival plot to a file
ggsave(paste0(result_dir, '/1.1_overall_surv.png'), 
       plot = surv_overall, 
       width = 4, height = 3, units = 'in', dpi = 150)

# Growth analysis
# Create a scatter plot of size at t0 versus size at t1
gr_overall <-
  ggplot(data  = grow_df, aes( x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Sampling inventory',
       subtitle = paste(script_prefix, species),
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1]))

ggsave(paste0(result_dir, '/1.2_overall_gr.png'), 
       plot = gr_overall, 
       width = 4, height = 3, units = 'in', dpi = 150)

# Recruitment analysis
# Create a scatter plot showing the relationship between 
# total parent plant area and number of recruits
rec_overall <- 
  ggplot(recr_df, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Sampling inventory',
       subtitle = paste(script_prefix, species),
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1]))  

# Save the recruitment plot as a PNG file
ggsave(paste0(result_dir, '/1.3_overall_rec.png'), 
       plot = rec_overall, 
       width = 4, height = 3, units = 'in', dpi = 150)


# Fit vital rate models for the mean IPM -----------------------------------
# Fit growth models to predict size at time t1 based on size at time t0
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
gr_mods              <- list(gr_mod_mean, gr_mod_mean_2, gr_mod_mean_3)
gr_dAIC_values       <- AICtab(gr_mods, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
gr_sorted_indices    <- order(gr_dAIC_values)
gr_mod_bestfit_index <- gr_sorted_indices[1 + gr_complex]
gr_mod_bestfit       <- gr_mods[[gr_mod_bestfit_index]]
gr_ranef             <- coef(gr_mod_bestfit)

# Predict size at time t1 using the mean growth model
grow_df$pred <- predict(gr_mod_bestfit, type = 'response')

# Plot observed size at time t1 against size at time t0 with the fitted line
# Define functions
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

grow_line <- 
  ggplot(grow_df, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, gr_ranef), 
                color = line_color_pred_fun(gr_ranef), 
                lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = paste(script_prefix, species))

# Plot predicted versus observed size at time t1
grow_pred <- 
  ggplot(grow_df, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw()

# Combine growth line and prediction plots
grow_overall_pred <- grow_line + grow_pred + plot_layout() 

# Save the growth prediction plot
ggsave(paste0(result_dir, '/2.1_overall_grow_pred_logs', 
              gr_mod_bestfit_index, '.png'), 
       plot = grow_overall_pred, 
       width = 8, height = 4, units = 'in', dpi = 150)

# Fit a model to assess variance in growth
# Fitted values from growth model
x         <- fitted(gr_mod_mean)  
# Squared residuals
y         <- resid(gr_mod_mean)^2  
# Non-linear model for variance
gr_var_m  <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0),
                 control = nls.control(maxiter = 1000, 
                                       tol = 1e-6, 
                                       warnOnly = TRUE) )  

# Survival
# Fit models to predict survival based on size at time t0
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
su_mods              <- list(su_mod_mean, su_mod_mean_2, su_mod_mean_3)
su_dAIC_values       <- AICtab(su_mods, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
su_sorted_indices    <- order(su_dAIC_values)
su_mod_bestfit_index <- su_sorted_indices[1 + su_complex]
su_mod_bestfit       <- su_mods[[su_mod_bestfit_index]]
su_ranef             <- coef(su_mod_bestfit)

# Generate predictions for survival across a range of sizes
surv_x <- seq(min(surv_df$logsize_t0, na.rm = T), 
              max(surv_df$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
surv_pred_df <- predictor_fun(surv_x, su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = surv_x, survives = .)

# Plot observed survival with fitted line
surv_line <- 
  ggplot() +
  geom_jitter(data = surv_df, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0, height = 0.25) + 
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = line_color_pred_fun(su_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = paste(script_prefix, species))

# Plot binned survival proportions with error bars
surv_bin <- 
  ggplot() +
  geom_point(data =  plot_binned_prop(df, 10, 
                                      logsize_t0, survives), 
             aes(x = logsize_t0, 
                 y = survives) ) +
  geom_errorbar(data =  plot_binned_prop(df, 10, 
                                         logsize_t0, survives), 
                aes(x = logsize_t0, 
                    ymin = lwr,
                    ymax = upr) ) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw()

# Combine survival plots
surv_overall_pred <- surv_line + surv_bin + plot_layout()

# Save the survival prediction plot
ggsave(paste0(result_dir, '/2.2_overall_surv_pred_logs', 
              su_mod_bestfit_index, '.png'), 
       plot = surv_overall_pred, 
       width = 8, height = 3, units = 'in', dpi = 150) 

## Recruitment
# Filter recruitment data to exclude NAs
recr_nona_nr_quad <- recr_df %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
rec_mod_mean <- MASS::glm.nb(nr_quad ~ 1, data = recr_nona_nr_quad)

# Generate predictions for recruitment
recr_nona_nr_quad <- 
  recr_nona_nr_quad %>% 
  mutate(pred_mod_mean = predict(rec_mod_mean, type = 'response')) 

# Summarize total number of recruits and predictions
rec_sums_df_m <- 
  recr_nona_nr_quad %>%
  summarize(nr_quad = sum(nr_quad),
            pred_mod_mean = sum(pred_mod_mean))

# Count number of adult individuals
indiv_m <- surv_df %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
repr_pc_m <- indiv_m %>%
  bind_cols(rec_sums_df_m) %>%
  mutate(repr_pc_mean = pred_mod_mean / n_adults) %>%
  mutate(repr_pc_obs = nr_quad / n_adults) %>%
  drop_na 


# Exporting parameter estimates ------------------------------------------------
# Growth
grow_fe <- data.frame(coefficient = names(coef(gr_mod_bestfit)),
                      value       = coef(gr_mod_bestfit))
var_m   <- data.frame(coefficient = names(coef(gr_var_m)),
                      value       = coef(gr_var_m))

grow_out <- Reduce(function(...) rbind(...), list(grow_fe, var_m)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(coefficient, 
                               grepl('Intercept', coefficient), 'b0'))

write.csv(grow_out, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_grow_pars_mean.csv'), 
          row.names = F)


# Survival
surv_fe <- data.frame(coefficient = names(coef(su_mod_bestfit)),
                      value       = coef(su_mod_bestfit))

surv_out<- Reduce(function(...) rbind(...), list(surv_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(coefficient, 
                               grepl('Intercept', coefficient), 'b0'))
write.csv(surv_out, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_surv_pars_mean.csv'), 
          row.names = F)

# Recruitment 
rec_size <- df %>% subset(recruit == 1)

others   <- data.frame(coefficient = c('rec_siz', 'rec_sd', 
                                       'max_siz', 'min_siz',
                                       'fecu_b0'),
                       value       = c(mean(log(rec_size$size_t0)), 
                                       sd(  log(rec_size$size_t0)),
                                       grow_df$logsize_t0 %>% max, 
                                       grow_df$logsize_t0 %>% min,
                                       repr_pc_m$repr_pc_mean))

write.csv(others, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_other_pars_mean.csv'), 
          row.names = F)


# Building the IPM from scratch ------------------------------------------------
extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = script_prefix,
  species = species,
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
  su_mod_bestfit_index = su_mod_bestfit_index,
  gr_mod_bestfit_index = gr_mod_bestfit_index
))

write.csv(pars, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_pars.csv'), 
          row.names = F)


# Function describing standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = gr_mod_bestfit_index) {
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
sx <- function(x, pars, num_pars = su_mod_bestfit_index) {
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
  recr_y  <- dnorm(y, pars$recr_sz, pars$recr_sd) * h
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
  group_by(year, quad) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  group_by(year, quad) %>%
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
        (if (su_mod_bestfit_index >= 1) surv_b1 * size_1   else 0) +
        (if (su_mod_bestfit_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (su_mod_bestfit_index >= 3) surv_b3 * size_1^3 else 0)),
    
    mu_g      = grow_b0 + 
      (if (gr_mod_bestfit_index >= 1) grow_b1 * size_1   else 0) +
      (if (gr_mod_bestfit_index >= 2) grow_b2 * size_1^2 else 0) +
      (if (gr_mod_bestfit_index >= 3) grow_b3 * size_1^3 else 0),
     
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

ipmr_p <- make_ipm(proto_ipm  = proto_ipm_p, 
                   iterations = 200)

lam_mean_ipmr <- lambda(ipmr_p)

lam_out       <- data.frame(coefficient = names(lam_mean_ipmr), 
                            value       = lam_mean_ipmr)
write.csv(lam_out, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_lambda_vec.csv'), 
          row.names = F)

lam_out_wide  <- as.list(pivot_wider(lam_out, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value'))
write.csv(lam_out_wide, 
          paste0(data_dir, '/', 
                 script_prefix, '_', sp_abb, '_lambda.csv'), 
          row.names = F)
