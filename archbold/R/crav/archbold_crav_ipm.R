# Building IPM - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.28


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
  tidyverse,
  # bbmle is for AICtab
  bbmle,
  # patchwork plot alingment
  patchwork,
  # binom.cofint for the survival plot
  binom) # , skimr, ipmr, binom, janitor, lme4


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
v_mod_set_gr <- c()
v_mod_set_su <- c()
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
df <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data_v2.csv')) %>% 
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


df_meta <- data.frame(variable = colnames(df)) %>% 
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


# Original year mean dataframe -------------------------------------------------
df_mean_og <- df %>%
  filter(s < 6 | is.na(s)) %>%
  group_by(site, quad, quad_id, plant, plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(s )), NA_real_, max (s,  na.rm = TRUE)),
    size_t0  = if_else(all(is.na(br)), NA_real_, max (br, na.rm = TRUE)),
    fruit    = if_else(all(is.na(fr)), NA_real_, max (fr, na.rm = TRUE)),   
    flower   = if_else(all(is.na(fl)), NA_real_, max (fl, na.rm = TRUE)),  
    fire_sev = if_else(
      all(is.na(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f))),
      NA_real_, 
      mean(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f), na.rm = TRUE)),
    .groups  = 'drop'
  ) %>% 
  ungroup()


# Base mean dataframe ----------------------------------------------------------
df_mean <- df_mean_og %>% 
  group_by(site, quad, quad_id, plant, plant_id) %>% 
  # Handle survival based on previous dormancy status
  # If the current status is dead or NA
  # Check if survived in the previous 1, 2, or 3 years
  # Check if survives in the next 1, 2, or 3 years
  mutate(
    survives = if_else(
      (survives == 0 | is.na(survives)) & 
        (lag(survives, 1) %in% c(1, 3, 5) | 
           lag(survives, 2) %in% c(1, 3, 5) | 
           lag(survives, 3) %in% c(1, 3, 5)) & 
        (lead(survives, 1) %in% c(1, 3, 5) | 
           lead(survives, 2) %in% c(1, 3, 5) | 
           lead(survives, 3) %in% c(1, 3, 5)),
      1,  
      survives
    ),
    
    # Set survival to 0 if plant dies after 3 consecutive years of dormancy
    survives = if_else(
      survives == 1 & 
        lead(survives, 1) == 0 & 
        lead(survives, 2) == 0 & 
        lead(survives, 3) == 0, 
      0,  
      survives
    )
  ) %>%
  
  # Define dormancy
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0) ~ 1,
      size_t0  >  0                  ~ 0, 
      TRUE ~ NA_real_ 
    ),
    # Generate a new column 'dormancy_count' that counts consecutive 1s
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 0 ~ 2,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 ~ 3,
      dormancy == 1 ~ 1,
      TRUE ~ dormancy)
  ) %>% 
  
  # Define recruits
  mutate(
    recruit = case_when(
      (survives == 3 | survives == 5) ~ 1, 
      TRUE ~ NA_real_  
    )
  ) %>% 
  
  mutate(
    # Make sure that no survival value exceeds 1
    survives = if_else(survives > 1, 1, survives), 
    
    # Handle missing size_t0 for dead plants by setting NA in survival column
    survives = if_else(survives == 0 & is.na(size_t0), NA, survives)
  ) %>% 
  
  ungroup() %>% 
  arrange(site, quad, plant, year) %>%
  group_by(site, quad, plant, plant_id) %>%
  
  # Set size_t1 based on survival; propagate size_t0 if survives, otherwise set to NA
  mutate(
    size_t1 = case_when(
      survives == 1 ~ lead(size_t0), 
      TRUE ~ NA_real_  
    )
  ) %>%
  
  ungroup() %>% 
  
  # Compute log-transformed sizes and their powers for modeling
  mutate(
    logsize_t0   = log(size_t0),     
    logsize_t1   = log(size_t1),    
    logsize_t0_2 = logsize_t0^2,     
    logsize_t0_3 = logsize_t0^3      
  ) %>%
  group_by(site, quad, plant, plant_id) %>%
  mutate(age = NA_real_) %>%  # Initialize age as NA for all
  
  # Apply for loop to calculate age
  mutate(age = {
    # Initialize age vector
    age_vector <- rep(NA_real_, n())
    
    # Iterate over each row in the group (each site, quad, and plant combination)
    for (i in 1:n()) {
      if (!is.na(recruit[i]) && recruit[i] == 1) {
        age_vector[i] <- 1  # Set age to 1 for recruits
      } else if (!is.na(dormancy[i]) && dormancy[i] >= 0) {
        # If plant survives, increment the age based on previous value
        if (i > 1 && !is.na(age_vector[i - 1])) {
          age_vector[i] <- age_vector[i - 1] + 1
        }
      }
    }
    
    # Return the computed age vector
    age_vector
  }) %>%
  ungroup() %>%
  arrange(site, quad, plant, plant_id, year) %>%
  group_by(site, quad, plant, plant_id) %>%
  mutate(
    fire_event = ifelse(!is.na(fire_sev) & fire_sev > 0, 1, 0),
    fire_gap = {
      gap <- numeric(n())
      counter <- NA_real_
      for (i in seq_along(fire_event)) {
        if (is.na(fire_event[i])) {
          # Keep NA if we haven't seen a fire yet
          gap[i] <- counter
        } else if (fire_event[i] == 1) {
          counter <- 0
          gap[i] <- counter
        } else {
          counter <- ifelse(is.na(counter), NA_real_, counter + 1)
          gap[i] <- counter
        }
      }
      gap
    }
  ) %>%
  ungroup()


# Survival data ----------------------------------------------------------------
df_su <- df_mean %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

fig_su_overall <- ggplot(
  data = plot_binned_prop(df_mean, 10, logsize_t0, survives)) +
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

ggplot(data = df_su) +
  geom_jitter(aes(x = logsize_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('Survival to time t1')) +
  theme(plot.subtitle = element_text(size = 8))


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
  geom_jitter(data = df_su, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, 
                                     y = survives),
            color = line_color_pred_fun(mod_su_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_su_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_mean, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df_mean, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        ymin = lwr,
        ymax = upr) ) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
fig_su <- fig_su_line + fig_su_bin + plot_layout()
fig_su


# Growth data ------------------------------------------------------------------
df_gr <- df_mean %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(plant_id, year, size_t0, size_t1, age,
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
mod_gr_1   <- lm(logsize_t1 ~ logsize_t0, 
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
                color = line_color_pred_fun(mod_gr_ranef), 
                lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

fig_gr_pred <- ggplot(
  df_gr, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw()

fig_gr <- fig_gr_line + fig_gr_pred + plot_layout() 
fig_gr


# Fit a model to assess variance in growth
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)  
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2  
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE) ) 


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


coef_misc   <- data.frame(coefficient = c('max_siz', 'min_siz'),
                          value       = c(df_gr$logsize_t0 %>% max, 
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
  L       = extr_value(coef_misc, 'min_siz'),
  U       = extr_value(coef_misc, 'max_siz'),
  mat_siz = 200,
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index
))


# Building the IPM -------------------------------------------------------------
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
  
  # Fertility matrix -----------------------------------------------------------
  Fmat        <- matrix(0, n, n)

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
