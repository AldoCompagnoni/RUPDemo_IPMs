# IPM mean - Archbold - Polygala lewtonii

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.05.20

# Study organism: Polygala lewtonii
# Link: 
# Meta data link: 
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
load_packages(MASS, tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4, GGally)


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
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
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
source('helper_functions/plot_binned_prop.R')


# Data -------------------------------------------------------------------------
df_og   <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_original_260519.csv')))
df_site <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_sitehist_260519.csv')))

df_meta <- data.frame(
  variable    = c('unique_ID', 'quad', 'patch', 'year', 'month', 'date', 'qsurv', 'stg', 'ht', 'mcd', 'st', 'flst'),
  description = c(
    'Unique identifier for each individual plant; Numeric',
    'Name of circular quadrat with radius of 0.25m; Numeric',
    'Name of an aggregate of several quads; Numeric',
    'Year data were collected; Numeric',
    'Month data were collected for March, June, September, and December. April was also used in the first sampling; Numeric',
    'Combination of year and month; MMYY',
    'Quarterly survival code. 0 or 1 refer to survival from the previous quarter. If a plant was dropped from monitoring, survival is NA and these plants should NOT be considered dead. 0 = dead; 1 = alive; 2 = missing; 3 = new adult; 5 = seedling',
    'Stage class from March census. If plant was a seedling during the previous year it is considered a seedling here. 1 = putative seedling; 2 = vegetative; 3 = reproductive adult',
    'Height from March census in cm; Numeric',
    'Maximum crown diameter from March census in cm; Numeric',
    'Total number of stems; Numeric',
    'Number of flowering stems; Numeric'),
  stringsAsFactors = FALSE)

df <- read.csv(file.path(dir_data,  paste0('ab_', v_sp_abb, '_df_workdata_260519.csv'))) %>% 
  # Define Disturbance
  mutate(disturbance = as.factor(if_else(burn_mean > 0, 1, 0)))


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

mod_su_bestfit <- glm(survives ~ logvol_t0 * recruits + logvol_t0_2:recruits + logvol_t0_3:recruits, 
                      data = df_su, family = 'binomial')

fig_su


# Growth -----------------------------------------------------------------------
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

mod_gr_bestfit       <- lm(logvol_t1 ~ logvol_t0 * recruits + logvol_t0_2, 
                           data = df_gr)

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


# Flower probability -----------------------------------------------------------
# We exclude all recruits form the probability to flower
df_fl <- df %>% 
  filter(!is.na(flower), !is.na(logvol_t0), !is.na(logvol_t0_2), !is.na(logvol_t0_3)) %>%
  filter(recruits == 0) %>% 
  dplyr::select(id, year, size_t0, flower, fl_nr, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                volume_t0, volume_t1, logvol_t0, logvol_t1, logvol_t0_2, logvol_t0_3,
                stage)

mod_fl_bestfit <-  glm(flower ~ logvol_t0 + logvol_t0_2,
                       data = df_fl, family = 'binomial')

fig_fl_all


# Number of flowers conditional on flowering -----------------------------------
df_fl_cond <- df_fl %>%
  filter(flower == 1) %>%
  filter(fl_nr %% 1 == 0)

mod_fl_n_bestfit <- glm.nb(fl_nr ~ logvol_t0 + logvol_t0_2 + logvol_t0_3,
                           data = df_fl_cond)

fig_fl_n_all


# Flowering stock to Recruit transition - Zero truncated, Bayesian -------------
#install.packages('brms') 
library(brms)

df_fs2r_0t <- df_fs2r %>%
  filter(recruit_count > 0)

# mod_fs2r_0t <- brm(
#   bf(recruit_count | trunc(lb = 1) ~ total_stocks + disturbance),
#   data = df_fs2r_0t,
#   family = negbinomial(link = "log"),
#   chains = 4,
#   cores = 4,
#   iter = 2000,
#   control = list(adapt_delta = 0.95))
# saveRDS(mod_fs2r_0t, file = file.path(dir_data, 'bayes_recruit_model_disturbance_260519.rds'))

mod_fs2r_0t <- readRDS(file.path(dir_data, 'bayes_recruit_model_disturbance_260519.rds'))

fig_fl2r_bay


# Recruitment data -------------------------------------------------------------
# Recruitment 
df_re_size <- df %>% subset(recruits == 1)


df_re <- df %>%
  group_by(year, site) %>%
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
      group_by(year, site) %>%
      summarise(nr_recs = sum(recruits, na.rm = TRUE), .groups = 'drop')
    
    left_join(df_cover, df_re, by = c('year', 'site'))
  }

df_re_mod <- df_re %>% filter(!is.na(nr_recs))
# Fit a negative binomial model for recruitment
mod_rec <- MASS::glm.nb(nr_recs ~ 1, data = df_re_mod)

# Generate predictions for recruitment
df_re_mod <- df_re_mod %>% 
  mutate(mod_pred = predict(mod_rec, type = 'response')) 

# Per-capita reproduction
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

# Flower probability
coef_fl_fe  <- data.frame(coefficient = names(coef(mod_fl_bestfit)),
                          value       =       coef(mod_fl_bestfit))

coef_fl <- Reduce(function(...) rbind(...), list(coef_fl_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Flower number conditional
coef_fln <- data.frame(coefficient = names(coef(mod_fl_n_bestfit)),
                       value       = coef(mod_fl_n_bestfit)) %>%
  mutate(coefficient = as.character(coefficient),
         coefficient = replace(coefficient, grepl('Intercept', coefficient), 'b0'))

# Miscellany
coef_misc   <- data.frame(coefficient = c('rec_siz', 'rec_sd',
                                          'fecu_b0',
                                          'max_siz', 'min_siz'),
                          value       = c(mean(log(df_re_size$size_t0), na.rm = T),
                                          sd(  log(df_re_size$size_t0), na.rm = T),
                                          df_fs2r_0t_pred$Estimate %>% mean(),
                                          df_gr$logvol_t0 %>% max,
                                          df_gr$logvol_t0 %>% min))

extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix      = v_script_prefix,
  species     = v_species,
  surv_b0     = extr_value(coef_su, 'b0'),
  surv_b1     = extr_value(coef_su, 'logvol_t0'),
  surv_b2     = extr_value(coef_su, 'logvol_t0_2'),
  surv_b3     = extr_value(coef_su, 'logvol_t0_3'),
  surv_br     = extr_value(coef_su, 'recruits1'),
  grow_b0     = extr_value(coef_gr, 'b0'),
  grow_b1     = extr_value(coef_gr, 'logvol_t0'),
  grow_b2     = extr_value(coef_gr, 'logvol_t0_2'),
  grow_b3     = extr_value(coef_gr, 'logvol_t0_3'),
  grow_br     = extr_value(coef_gr, 'recruits1'),
  grow_b1_br  = extr_value(coef_gr, 'logvol_t0:recruits1'),
  grow_b2_br0 = extr_value(coef_gr, 'recruits0:logvol_t0_2'),
  grow_b2_br1 = extr_value(coef_gr, 'recruits0:logvol_t0_2'),
  a           = extr_value(coef_gr, 'a'),
  b           = extr_value(coef_gr, 'b'),
  fl_b0       = extr_value(coef_fl, 'b0'),
  fl_b1       = extr_value(coef_fl, 'logvol_t0'),
  fl_b2       = extr_value(coef_fl, 'logvol_t0_2'),
  fl_b3       = extr_value(coef_fl, 'logvol_t0_3'),
  fr_b0       = extr_value(coef_fl, 'b0'),
  fln_b0 = extr_value(coef_fln, 'b0'),
  fln_b1 = extr_value(coef_fln, 'logvol_t0'),
  fln_b2 = extr_value(coef_fln, 'logvol_t0_2'),
  fln_b3 = extr_value(coef_fln, 'logvol_t0_3'),
  fecu_b0     = extr_value(coef_misc, 'fecu_b0'),
  recr_sz     = extr_value(coef_misc, 'rec_siz'),
  recr_sd     = extr_value(coef_misc, 'rec_sd'),
  L           = extr_value(coef_misc, 'min_siz'),
  U           = extr_value(coef_misc, 'max_siz'),
  mat_siz     = 200,
  mod_su_index   = v_mod_su_index,
  mod_gr_index   = v_mod_gr_index,
  mod_fl_index   = v_mod_fl_index,
  mod_fl_n_index = v_mod_fl_n_index))


# Function describing standard deviation of growth model -----------------------
# Function describing the invert logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

# Survival of x-sized individual to time t1
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  for (i in 1:num_pars) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^i
    }
  }
  inv_logit(survival_value)
}

# Growth variation
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- pars$grow_b0
  for (i in 1:num_pars) {
    param_name <- paste0("grow_b", i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  if (!is.null(pars$grow_b1_br)) {
    mean_value <- mean_value
  }
  sd_value <- grow_sd(x, pars)
  dnorm(y, mean = mean_value, sd = sd_value)
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


# Flowering of x-sized individual at time t0
fl_n_x <- function(x, pars, num_pars = v_mod_fl_n_index) {
  val <- pars$fln_b0
  for (i in 1:num_pars) {
    param <- paste0('fln_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }
  exp(val)
}

# Recruitment size distribution at time t1
re_y_dist <- function(y, pars) {
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
}

# F-kernel
fyx <- function(y, x, pars) {
  fl_x(x, pars) *
    fl_n_x(x, pars) *
    pars$fecu_b0 *
    re_y_dist(y, pars)
}


# Kernel -----------------------------------------------------------------------
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
  Smat <- sx(y, pars)
  
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

lambda_ipm <- function(pars) {
  Re(eigen(kernel(pars)$k_yx)$value[1])
}

# mean population growth rate
lam_mean <- lambda_ipm(pars)
lam_mean


# Observed population growth ---------------------------------------------------
df_counts_year <- df %>%
  group_by(year) %>%
  filter(!is.na(survives)) %>%
  summarise(n = n())

# Then compute observed lambda
lam_obs_y <- df_counts_year$n[-1] / df_counts_year$n[-nrow(df_counts_year)]
lam_obs_mean <- mean(lam_obs_y, na.rm = TRUE)
lam_obs_mean