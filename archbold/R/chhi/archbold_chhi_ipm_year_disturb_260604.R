# IPM year-specific with fire - Archbold -  - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.05.28


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
  # GLMM: glmer function
  lme4, 
  janitor) # , skimr, ipmr, binom


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head             <- c('archbold')
# Define species
v_species          <- c('Chrysopsis highlandsensis')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c()
# Years that we want to remove from the analysis
v_years_re         <- c()


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
# fig_su
v_mod_set_gr <- c()
# fig_gr
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
# function to plot your survival data 'binned' per year (instead of 'jittered')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df_og <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>% 
  janitor::clean_names()


# Fire data --------------------------------------------------------------------
df_fire <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_fire.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    year = burn_yr,
    fire = treatment) %>% 
  select(!c(year0, notes))

# Mean data frame --------------------------------------------------------------
df <- df_og %>% 
  rename(
    plant_id = identifier,
    year     = year0,
    survival = survival_1) %>%  
  mutate(
    plant_id = as.factor(plant_id)) %>%
  arrange(site, plant_id, year, survival) %>%
  mutate(recruit = ifelse(astg == 1, 1, 0)) %>%
  group_by(site , plant_id, year) %>%
  summarise(
    survives = if (all(is.na(survival))) NA_real_ else min(survival, na.rm = TRUE),
    size_t0  = if (all(is.na(dia)))      NA_real_ else max(dia,      na.rm = TRUE),
    size_t1  = if (all(is.na(dia_1)))    NA_real_ else max(dia_1,    na.rm = TRUE),
    flower   = if (all(is.na(hd)))       NA_real_ else max(hd,       na.rm = TRUE),
    recruit  = if (all(is.na(recruit)))  NA_real_ else min(recruit,  na.rm = TRUE),
    .groups = 'drop') %>%
  mutate(
    fl_nr  = flower,
    flower = if_else(flower > 0, 1, flower)) %>% 
  mutate(
    logsize_t0   = log(size_t0),
    logsize_t1   = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  full_join(df_fire, by = c('site', 'year')) %>%
  mutate(
    fire = case_when(
      is.na(fire)    ~ 'No fire',
      fire == 'burn' ~ 'Fire',
      TRUE           ~ NA_character_),
    fire = factor(fire, levels = c('No fire', 'Fire')))


# Survival data ----------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives)) %>%
  #  filter(size_t0 != 0) %>% # what do I need this for and why does it give me a different length of data????
  dplyr::select(plant_id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
                fire) %>%
  # Complete cases only for the AICctab
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year)


# Survival model ---------------------------------------------------------------
# GLMM; binomial
# Intercept model
mod_su_0 <- glmer(
  survives ~ fire + (1 | year),
  data = df_su, family = binomial) 
# Linear size
mod_su_1 <- glmer(
  survives ~ logsize_t0 + fire + (1 | year),
  data = df_su, family = binomial) 
# Quadratic size
mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + fire + (1 | year),
  data = df_su, family = binomial)  
# Cubic size
mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire + (1 | year),
  data = df_su, family = binomial)  

# Compare models using AIC
mods_su      <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- bbmle::AICctab(mods_su, weights = T, sort = F)$dAIC

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

mod_su_best <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef   <- coef(mod_su_best)$year


# FIX YEAR TYPE
df_su$year <- factor(df_su$year)

# remake split object
surv_bin_yrs <- split(df_su, df_su$year)

v <- names(surv_bin_yrs)

surv_yr_plots <- function(i) {
  
  surv_temp <- as.data.frame(surv_bin_yrs[[i]])
  
  x_temp <- seq(
    min(surv_temp$logsize_t0, na.rm = TRUE),
    max(surv_temp$logsize_t0, na.rm = TRUE),
    length.out = 100)
  
  # NO FIRE
  pred_no_fire <- data.frame(
    logsize_t0   = x_temp,
    logsize_t0_2 = x_temp^2,
    logsize_t0_3 = x_temp^3,
    fire         = factor("No fire", levels = levels(df_su$fire)),
    year         = factor(v[i], levels = levels(df_su$year)))
  
  pred_no_fire$surv <- predict(
    mod_su_best,
    newdata = pred_no_fire,
    type = "response",
    re.form = NULL)
  
  # FIRE
  pred_fire <- data.frame(
    logsize_t0   = x_temp,
    logsize_t0_2 = x_temp^2,
    logsize_t0_3 = x_temp^3,
    fire         = factor("Fire", levels = levels(df_su$fire)),
    year         = factor(v[i], levels = levels(df_su$year)))
  
  pred_fire$surv <- predict(
    mod_su_best,
    newdata = pred_fire,
    type = "response",
    re.form = NULL)
  
  # BINNED MEAN SURVIVAL POINTS
  surv_pts <- surv_temp %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, fire) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      survives   = mean(survives, na.rm = TRUE),
      n           = n(),
      .groups = "drop") %>%
    filter(n > 0)
  
  ggplot() +
    # regular points
    geom_point(data = surv_pts %>% filter(fire == "No fire"),
               aes(logsize_t0, survives),
               size = 1, color = "black") +
    # fire points
    geom_point(data = surv_pts %>% filter(fire == "Fire"),
               aes(logsize_t0, survives),
               size = 1, color = "red") +
    # no fire prediction
    geom_line(data = pred_no_fire,
              aes(logsize_t0, surv),
              color = "black", lwd = 1) +
    # fire prediction
    geom_line(aes(logsize_t0, surv),
              data = pred_fire, color = "red", lwd = 1) +
    labs(title = v[i],
         x     = expression('log(size)'[t0]),
         y     = expression('Survival probability '[t1])) +
    ylim(0,1) +
    theme_bw()
}

surv_yrs <- lapply(1:length(surv_bin_yrs), surv_yr_plots)

g_surv_years <- wrap_plots(surv_yrs) +
  plot_layout(ncol = 4)

g_surv_years


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  filter(size_t0 != 0) %>%
  filter(size_t1 != 0) %>%
  dplyr::select(
    plant_id, year, size_t0, size_t1,
    logsize_t0, logsize_t1,
    logsize_t0_2, logsize_t0_3,
    fire) %>%
  filter(is.finite(logsize_t1)) %>%
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year) %>% 
  mutate(year = as.factor(year))


# Growth model -----------------------------------------------------------------
# Intercept model 
mod_gr_0 <- lmer(
  logsize_t1 ~ fire + (logsize_t0 | year), 
  data = df_gr)
# Linear size
mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + fire + (logsize_t0 | year), 
  data = df_gr)
# Quadratic size
mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + fire + (logsize_t0 | year),
  data = df_gr)
# Cubic size
mod_gr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire + (logsize_t0 | year),
  data = df_gr)
# Cubic size
mod_gr_4 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (logsize_t0 | year),
  data = df_gr)

# Compare models using AIC
mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- bbmle::AICctab(mods_gr, weights = T, sort = F)$dAIC

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

mod_gr_best <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef   <- coef(mod_gr_best)$year


grow_yr_plots <- function(i){
  
  grow_temp <- df_gr %>%
    filter(year == i)
  
  x_temp <- seq(
    min(grow_temp$logsize_t0, na.rm = TRUE),
    max(grow_temp$logsize_t0, na.rm = TRUE),
    length.out = 100)
  
  # NO FIRE
  pred_no_fire <- data.frame(
    logsize_t0   = x_temp,
    logsize_t0_2 = x_temp^2,
    logsize_t0_3 = x_temp^3,
    fire         = factor("No fire", levels = levels(df_gr$fire)),
    year         = factor(i, levels = levels(df_gr$year)))
  
  pred_no_fire$pred <- predict(
    mod_gr_best,
    newdata = pred_no_fire,
    re.form = NULL)
  
  # FIRE
  pred_fire <- data.frame(
    logsize_t0   = x_temp,
    logsize_t0_2 = x_temp^2,
    logsize_t0_3 = x_temp^3,
    fire         = factor("Fire", levels = levels(df_gr$fire)),
    year         = factor(i, levels = levels(df_gr$year)))
  
  pred_fire$pred <- predict(
    mod_gr_best,
    newdata = pred_fire,
    re.form = NULL)
  
  # BINNED MEANS
  grow_pts <- grow_temp %>%
    mutate(bin = cut(logsize_t0, breaks = 8)) %>%
    group_by(bin, fire) %>%
    summarise(
      logsize_t0 = mean(logsize_t0, na.rm = TRUE),
      logsize_t1 = mean(logsize_t1, na.rm = TRUE),
      n = n(),
      .groups = "drop") %>%
    filter(n > 0)
  
  ggplot() +
    # no fire points
    geom_point(
      data = grow_pts %>% filter(fire == "No fire"),
      aes(logsize_t0, logsize_t1),
      color = "black",
      size = 1) +
    # fire points
    geom_point(
      data = grow_pts %>% filter(fire == "Fire"),
      aes(logsize_t0, logsize_t1),
      color = "red",
      size = 1) +
    # no fire line
    geom_line(aes(logsize_t0, pred),
              data = pred_no_fire, color = "black", lwd = 1) +
    # fire line
    geom_line(aes(logsize_t0, pred), 
              data = pred_fire, color = "red", lwd = 1) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "blue",
      lty = 2) +
    labs(
      title = i,
      x = expression('log(size)'[t0]),
      y = expression('log(size)'[t1])) +
    theme_bw() +
    theme(text = element_text(size = 5))
}

grow_yrs <- lapply(sort(unique(df_gr$year)), grow_yr_plots)

g_grow_years <- wrap_plots(grow_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Growth - year specific",
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9)))

g_grow_years


# Growth variance year specific ------------------------------------------------
x      <- fitted(mod_gr_best)
y      <- resid( mod_gr_best)^2
gr_var <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0))


# Flower data ------------------------------------------------------------------
df_fl <- df %>%
  filter(!is.na(flower)) %>%
  dplyr::select(plant_id, year, size_t0, flower, fl_nr, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire) %>%
  mutate(flower = if_else(flower > 0, 1, flower)) %>%
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year)


# Flower model -----------------------------------------------------------------
# GLMM; binomial
# Intercept model
mod_fl_0 <- glmer(
  flower ~ fire + (1 | year),
  data = df_fl, family = binomial) 
# Linear size term
mod_fl_1 <- glmer(
  flower ~ logsize_t0 + fire + (1 | year),
  data = df_fl, family = binomial) 
# Quadratic size term
mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + fire + (1 | year),
  data = df_fl, family = binomial)  
# Cubic size term
mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire + (1 | year),
  data = df_fl, family = binomial)  

# Compare models using AIC
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- bbmle::AICctab(mods_fl, weights = T, sort = F)$dAIC

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

mod_fl_best <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef   <- coef(mod_fl_best)$year


# Flower number conditional on flowering ---------------------------------------
df_fl_cond <- df_fl %>%
  filter(flower == 1) %>%
  filter(!is.na(fl_nr), fl_nr > 0)

mod_fl_n_0 <- MASS::glm.nb(
  fl_nr ~ fire,
  data = df_fl_cond)

mod_fl_n_1 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + fire,
  data = df_fl_cond)

mod_fl_n_2 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + fire,
  data = df_fl_cond)

mod_fl_n_3 <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire,
  data = df_fl_cond)

mods_fl_n <- list(mod_fl_n_0, mod_fl_n_1, mod_fl_n_2, mod_fl_n_3)

mods_fl_n_dAIC <- bbmle::AICctab(mods_fl_n, weights = TRUE, sort = FALSE)$dAIC

mods_fl_n_sorted <- order(mods_fl_n_dAIC)

mod_fl_n_bestfit <- mods_fl_n[[mods_fl_n_sorted[1]]]
v_mod_fl_n_index <- mods_fl_n_sorted[1] - 1


# Recruit size distribution -----------------------------------------------

df_recr <- df %>%
  filter(recruit == 1) %>%
  filter(size_t0 > 0) %>%
  mutate(logsize_t0 = log(size_t0))

rc_sz <- data.frame(
  mean = mean(df_recr$logsize_t0, na.rm = TRUE),
  sd   = sd(df_recr$logsize_t0, na.rm = TRUE))

# Exporting parameter estimates ------------------------------------------------
# Survival
# Get the coefficients matrix
su_coef_matrix <- coef(mod_su_best)$year

# Initialize a list to store data frames
su_data_frames <- list(
  data.frame(coefficient = paste0('year_', rownames(su_coef_matrix)), 
             value = su_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0_', rownames(su_coef_matrix)), 
             value = su_coef_matrix[, 'logsize_t0']))

# Loop to create additional data frames if needed
for (i in 2:ncol(su_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(su_coef_matrix)) {
    su_data_frames[[length(su_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, rownames(su_coef_matrix)),
      value = su_coef_matrix[, column_name])
  }
}

su_data_frames[[length(su_data_frames) + 1]] <- data.frame(
  coefficient = paste0('fire'),
  value = su_coef_matrix[, grep("^fire", colnames(su_coef_matrix), value = TRUE)])

# Combine data frames and mutate the coefficient column
su_out_yr <- Reduce(rbind, su_data_frames) %>%
  mutate(coefficient = as.character(coefficient))


# Growth
# Get coefficients matrix
gr_coef_matrix <- coef(mod_gr_best)$year

# Initialize a list for year and logsize data frames, including var_fe
gr_data_frames <- list(
  data.frame(coefficient = names(coef(gr_var)), value = coef(gr_var)),
  data.frame(coefficient = paste0('year_', rownames(gr_coef_matrix)), 
             value = gr_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0_', rownames(gr_coef_matrix)), 
             value = gr_coef_matrix[, 'logsize_t0']))

# Loop to create additional logsize data frames
for (i in 2:ncol(gr_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(gr_coef_matrix)) {
    gr_data_frames[[length(gr_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, '_', rownames(gr_coef_matrix)),
      value = gr_coef_matrix[, column_name])
  }
}

gr_data_frames[[length(gr_data_frames) + 1]] <- data.frame(
  coefficient = paste0('fire'),
  value = gr_coef_matrix[, grep("^fire", colnames(gr_coef_matrix), value = TRUE)])

# Combine all data frames using Reduce and mutate the coefficient column
gr_out_yr <- Reduce(function(...) rbind(...), gr_data_frames) %>%
  mutate(coefficient = as.character(coefficient))





# df constant parameters, fixed effects estimates, and mean parameter estimates
# Seedbank constants
mean_rec <- df %>%
  group_by(year) %>%
  summarise(rec = sum(recruit, na.rm = TRUE), .groups = "drop") %>%
  summarise(mean_rec = mean(rec)) %>%
  pull(mean_rec)

seedbank_size <- 1e6
seed_input <- 200

emerg_rate <- mean_rec / seedbank_size
seed_surv <- seedbank_size / (seedbank_size - mean_rec + seed_input)

constants <- data.frame(
  coefficient = c(
    'recr_sz',
    'recr_sd',
    'a',
    'b',
    'L',
    'U',
    'mat_siz',
    'seedbank_size',
    'emerg_rate',
    'seed_surv'
  ),
  value = c(
    rc_sz$mean,
    rc_sz$sd,
    as.numeric(coef(gr_var)[1]),
    as.numeric(coef(gr_var)[2]),
    min(df_gr$logsize_t0),
    max(df_gr$logsize_t0),
    200,
    seedbank_size,
    emerg_rate,
    seed_surv
  )
)

# Create the data frame dynamically based on the number of fixed effects
su_fe <- data.frame(
  coefficient = c(paste0('surv_b', 0:(length(fixef(mod_su_best)) - 2)), 'surv_bf'),
  value       = fixef(mod_su_best))

gr_fe <- data.frame(
  coefficient = c(paste0('grow_b', 0:(length(fixef(mod_gr_best)) - 2)), 'grow_bf'),
  value       = fixef(mod_gr_best))

fl_fe <- data.frame(
  coefficient = c(paste0('fl_b', 0:(length(fixef(mod_fl_best)) - 2)), 'fl_bf'),
  value       = fixef(mod_fl_best))

fln_fe <- data.frame(
  coefficient = c(paste0('fln_b', 0:(length(coef(mod_fl_n_bestfit)) - 2)), 'fln_bf'),
  value       = coef(mod_fl_n_bestfit))


pars_cons <- Reduce(function(...) rbind(...), 
                    list(su_fe, gr_fe, fl_fe, fln_fe, constants)) %>%
  mutate(coefficient = as.character(coefficient))

pars_cons_wide <- as.list(
  pivot_wider(pars_cons,
              names_from = 'coefficient',
              values_from = 'value')
)

pars_mean <- pars_cons_wide



# DF varying parameters
# Function to create coefficient data frames dynamically
create_coef_df <- function(model, prefix) {
  coef_matrix <- coef(model)$year
  
  column_map <- c(
    '(Intercept)'   = '0',
    'logsize_t0'    = '1',
    'logsize_t0_2'  = '2',
    'logsize_t0_3'  = '3',
    'fire'          = 'f' )
  
  # Only keep column_map entries that exist in coef_matrix
  present_cols <- intersect(names(column_map), colnames(coef_matrix))
  column_map <- column_map[present_cols]
  
  lapply(names(column_map), function(col_name) {
    suffix <- column_map[[col_name]]
    data.frame(
      coefficient = paste0(prefix, suffix, '_', rownames(coef_matrix)),
      value = coef_matrix[, col_name]
    )
  })
}



# Create data frames for survival and growth models
su_data_frames <- create_coef_df(mod_su_best, 'surv_b')
gr_data_frames <- create_coef_df(mod_gr_best, 'grow_b')

# Combine all data frames into one
pars_var <- Reduce(rbind, c(
  su_data_frames, gr_data_frames))

pars_var_wide <- as.list(pivot_wider(pars_var, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value'))


# Building the year-specific IPMs from scratch ---------------------------------
# Functions
# Inverse logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

# Survival of x-sized individual to time t1
sx <- function(x,
               pars,
               fire = 0,
               num_params = mod_su_index_bestfit) {
  
  survival_value <- pars$surv_b0
  
  for (i in 1:num_params) {
    
    param_name <- paste0('surv_b', i)
    
    if (!is.null(pars[[param_name]])) {
      
      survival_value <- survival_value +
        pars[[param_name]] * x^i
    }
  }
  
  survival_value <- survival_value +
    pars$surv_bf * fire
  
  inv_logit(survival_value)
}

# Standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b * x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x,
                y,
                pars,
                fire = 0,
                num_params = mod_gr_index_bestfit) {
  mean_value <- 0 # Why is the mean value 0?
  for (i in 0:num_params) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  mean_value <- mean_value + pars$grow_bf * fire
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}

# Transition of x-sized individual to y-sized individual at time t1
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

fl_x <- function(x, pars, fire = 0, num_params = v_mod_fl_index) {
  val <- pars$fl_b0
  for (i in 1:num_params) {
    param_name <- paste0("fl_b", i)
    if (!is.null(pars[[param_name]])) {
      val <- val + pars[[param_name]] * x^i
    }
  }
  val <- val + pars$fl_bf * fire
  inv_logit(val)
}

fl_n_x <- function(x, pars, fire = 0, num_params = v_mod_fl_n_index) {
  val <- pars$fln_b0
  for (i in 1:num_params) {
    param_name <- paste0("fln_b", i)
    if (!is.null(pars[[param_name]])) {
      val <- val + pars[[param_name]] * x^i
    }
  }
  val <- val + pars$fln_bf * fire
  exp(val)
}

seed_prod_x <- function(x, pars, fire = 0) {
  fl_x(x, pars, fire) *
    fl_n_x(x, pars, fire)
}

# Recruit size distribution
re_y_dist <- function(y, pars) {
  dnorm(y,
        mean = pars$recr_sz,
        sd   = pars$recr_sd)
}

# Kernel
kernel <- function(pars, fire = 0) {
  
  n      <- pars$mat_siz
  L      <- pars$L
  U      <- pars$U
  h      <- (U - L) / n
  b      <- L + c(0:n) * h
  y      <- 0.5 * (b[1:n] + b[2:( n + 1 )])
  
  Smat   <- c()
  Smat <- sx(y, pars, fire = fire)
  
  Gmat   <- matrix(0, n, n)
  Gmat[] <- t(
    outer(y, y,
          Vectorize(function(x, y)
            gxy(x, y, pars, fire = fire)))) * h
  
  Tmat   <- matrix(0, n, n)
  
  for(i in 1:(n / 2)) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  for(i in (n / 2 + 1):n) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Seedbank blocks
  S11 <- matrix(pars$seed_surv, 1, 1)
  S12 <- matrix(
    seed_prod_x(y, pars, fire) * h,
    nrow = 1)
  S21 <- matrix(
    pars$emerg_rate * re_y_dist(y, pars),
    ncol = 1)
  S22 <- Tmat
  
  k_yx <- rbind(
    cbind(S11, S12),
    cbind(S21, S22))
  
  return(list(
    k_yx    = k_yx,
    Tmat    = Tmat,
    Gmat    = Gmat,
    meshpts = y))
  
}


# Mean population growth rate --------------------------------------------------
pars_mean <- pars_cons_wide

lambda_ipm <- function(i, fire = 0) {
  Re(eigen(kernel(i, fire = fire)$k_yx)$value[1])
}

lambda_ipm(pars_mean, fire = 0)
lambda_ipm(pars_mean, fire = 1)



# population growth rates for each year -------------------------------------

years_v <- sort(unique(df$year))

pars_yr <- vector(mode = 'list', length = length(years_v))

extr_value_list <- function(x, field) {
  as.numeric(x[paste0(field)] %>% unlist())
}

safe_get <- function(x, field) {
  val <- x[[field]]
  if (is.null(val)) return(NULL)
  as.numeric(val)
}

prep_pars <- function(i, num_surv_params, num_grow_params) {
  
  yr_now <- years_v[i]
  
  pars_year <- list()
  
  # add all constant / fixed parameters
  for(nm in names(pars_cons_wide)) {
    pars_year[[nm]] <- safe_get(pars_cons_wide, nm)
  }
  
  # replace survival year-specific params
  for (j in 0:num_surv_params) {
    nm <- paste0("surv_b", j)
    yr_nm <- paste(nm, yr_now, sep = "_")
    val <- safe_get(pars_var_wide, yr_nm)
    if (!is.null(val)) pars_year[[nm]] <- val
  }
  
  # replace growth year-specific params
  for (j in 0:num_grow_params) {
    nm <- paste0("grow_b", j)
    yr_nm <- paste(nm, yr_now, sep = "_")
    val <- safe_get(pars_var_wide, yr_nm)
    if (!is.null(val)) pars_year[[nm]] <- val
  }
  
  pars_year
}


# remove bad years
years_v <- sort(unique(as.numeric(df$year)))

pars_yr <- lapply(
  1:length(years_v),
  num_surv_params = v_mod_su_index,
  num_grow_params = v_mod_gr_index,
  prep_pars)

contains_numeric0 <- sapply(pars_yr, function(x) {
  any(sapply(x, function(y) identical(y, numeric(0))))
})

pars_yr <- pars_yr[!contains_numeric0]
years_v <- years_v[!contains_numeric0]

# fire vector ---------------------------------------------------------------

fire_lookup <- df %>%
  distinct(year, fire) %>%
  mutate(fire_num = ifelse(fire == "Fire", 1, 0))

fire_v <- fire_lookup$fire_num[
  match(years_v, fire_lookup$year)
]

# lambda per year -----------------------------------------------------------

calc_lambda <- function(i) {
  
  yr <- years_v[i]
  
  fire_y <- fire_lookup %>%
    filter(year == yr) %>%
    summarise(fire_num = max(fire_num, na.rm = TRUE)) %>%
    pull(fire_num)
  
  Re(eigen(kernel(pars_yr[[i]], fire = fire_y)$k_yx)$value[1])
}

lambdas_yr <- sapply(
  seq_along(pars_yr),
  calc_lambda)

names(lambdas_yr) <- years_v
lambdas_yr


# observed population growth rates -----------------------------------------

pop_counts_t0 <- df %>%
  group_by(year, site) %>%
  summarise(n_t0 = n(), .groups = 'drop') %>%
  mutate(year = year + 1)

pop_counts_t1 <- df %>%
  group_by(year, site) %>%
  summarise(n_t1 = n(), .groups = 'drop')

pop_counts <- left_join(
  pop_counts_t0,
  pop_counts_t1,
  by = c("year", "site")) %>%
  mutate(year = year - 1) %>%
  drop_na() %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0),
    n_t1 = sum(n_t1),
    .groups = 'drop') %>%
  mutate(
    obs_pgr = n_t1 / n_t0) %>%
  full_join(
    data.frame(
      year   = as.numeric(names(lambdas_yr)),
      lambda = unlist(lambdas_yr)),
    by = 'year') %>%
  drop_na()


pop_counts <- pop_counts %>%
  left_join(
    fire_lookup %>%
      select(year, fire),
    by = "year")


g_mod_vs_obs <- ggplot(pop_counts) +
  geom_point(
    aes(x = lambda, y = obs_pgr, color = fire),
    size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  labs(
    title = 'Observed vs modeled population growth',
    x = expression('Modeled '*lambda),
    y = 'Observed population growth rate',
    color = 'Fire') +
  theme_classic()

g_mod_vs_obs

pop_counts %>% summarize(mean(obs_pgr))

pop_counts %>%
  summarise(
    geometric_mean_pgr = exp(mean(log(obs_pgr), na.rm = TRUE)))


# Observed PGR vs asymptotic lambda vs projected lambda ------------------------

make_initial_n <- function(year0, pars){
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  df0 <- df %>%
    filter(year == year0,
           !is.na(logsize_t0),
           is.finite(logsize_t0))
  
  adult_counts <- hist(
    pmin(pmax(df0$logsize_t0, L), U),
    breaks = seq(L, U, length.out = n + 1),
    plot = FALSE,
    include.lowest = TRUE
  )$counts
  
  adult_density <- adult_counts / h
  
  c(seedbank = pars$seedbank_size, adult_density)
}

get_fire_y <- function(yr) {
  fire_lookup %>%
    filter(year == yr) %>%
    summarise(fire_num = max(fire_num, na.rm = TRUE)) %>%
    pull(fire_num)
}

project_one_year <- function(yr) {
  
  pars_y <- pars_yr[[match(yr, years_v)]]
  fire_y <- get_fire_y(yr)
  
  n_obs <- make_initial_n(yr, pars_y)
  K     <- kernel(pars_y, fire = fire_y)$k_yx
  
  h <- (pars_y$U - pars_y$L) / pars_y$mat_siz
  
  n_proj <- K %*% n_obs
  
  data.frame(
    year = yr,
    asym_lambda = Re(eigen(K)$values[1]),
    proj_lambda = (sum(n_proj[-1]) * h) / (sum(n_obs[-1]) * h)
  )
}

df_lams <- bind_rows(lapply(1999:2023, project_one_year))

df_obs_pgr <- df %>%
  filter(!is.na(logsize_t0),
         is.finite(logsize_t0)) %>%
  group_by(year) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(obs_pgr = lead(n) / n) %>%
  filter(year %in% 1999:2023)

df_compare <- df_obs_pgr %>%
  left_join(df_lams, by = "year") %>%
  left_join(
    fire_lookup %>%
      group_by(year) %>%
      summarise(
        fire = if_else(any(fire == "Fire", na.rm = TRUE), "Fire", "No fire"),
        .groups = "drop"),
    by = "year")

df_plot <- df_compare %>%
  select(year, obs_pgr, asym_lambda, proj_lambda, fire) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = "lambda_type",
    values_to = "lambda") %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = "Asymptotic lambda",
      proj_lambda = "Projected lambda from observed size distribution"))

g_mod_vs_obs <- ggplot(df_plot, aes(x = lambda, y = obs_pgr, color = fire)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = "free_x") +
  labs(
    title = "Observed population growth vs modeled lambda",
    x = expression("Modeled " * lambda),
    y = "Observed population growth rate",
    color = "Fire") +
  theme_classic()

g_mod_vs_obs