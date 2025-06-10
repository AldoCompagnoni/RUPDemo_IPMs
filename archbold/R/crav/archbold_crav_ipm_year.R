# IPM year-specific - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.06.04


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
  # GLMM: glmer function
  lme4) # , skimr, ipmr, binom, janitor


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head             <- c('archbold')
# Define species
v_species          <- c('Crotalaria avonensis')
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
# function to plot your survival data 'binned' per year (instead of 'jittered')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
# Original data
df_org  <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df_org.csv'))) %>%  
  mutate(
    plant_id = as.factor(plant_id),
    quad_id  = as.factor(quad_id),
    site     = as.factor(site),
    quad     = as.factor(quad),
    mp       = as.factor(mp),
    plant    = as.factor(plant),
    caged    = as.factor(caged),
    veg      = as.factor(veg))

# Meta data
df_meta <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df_meta.csv')))

# Aggrigated annual data 
df_agg  <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df_agg.csv'))) %>%  
  mutate(
    plant_id = as.factor(plant_id),
    quad_id  = as.factor(quad_id),
    site     = as.factor(site),
    quad     = as.factor(quad),
    plant    = as.factor(plant))


# Working data -----------------------------------------------------------------
df <- read.csv(
  file.path(dir_data, paste0('ab_', v_sp_abb, '_df.csv'))) %>%  
  mutate(
    plant_id = as.factor(plant_id),
    quad_id  = as.factor(quad_id),
    site     = as.factor(site),
    quad     = as.factor(quad),
    plant    = as.factor(plant))


# Survival data ----------------------------------------------------------------
df_su <- df %>% 
  filter(!is.na(survives)) %>%
#  filter(size_t0 != 0) %>% # what do I need this for and why does it give me a different length of data????
  dplyr::select(plant_id, year, size_t0, survives, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>%
  # Complete cases only for the AICctab
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year)


# Growth data ------------------------------------------------------------------
df_gr <- df %>% 
  filter(size_t0 != 0) %>%
  #    filter(size_t1 != 0) %>% # what do I need this for and why does it give me a different length of data????
  dplyr::select(plant_id, year, size_t0, size_t1, age,
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>%
  # Complete cases only for the AICctab
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year)


# Flower data ------------------------------------------------------------------
df_fl <- df %>% 
  filter(!is.na(flower)) %>%
  dplyr::select(plant_id, year, size_t0, flower, size_t1, 
                logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>% 
  mutate(flower = if_else(flower > 0, 1, flower)) %>%
  # Complete cases only for the AICctab
  drop_na(logsize_t0, logsize_t0_2, logsize_t0_3, year)


# Fruit data -------------------------------------------------------------------
df_fr <- df %>%
  filter(!is.na(fruit), !is.na(size_t0))


# Recruitment data -------------------------------------------------------------
df_re <- df %>%
  group_by(year, site, quad) %>%
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
      group_by(year, site, quad) %>%
      summarise(nr_quad = sum(recruit, na.rm = TRUE), .groups = "drop")
    
    left_join(df_cover, df_re, by = c("year", "site", "quad"))
  }


# Removing year with too few data ----------------------------------------------
# Years original
v_years_og <- sort(unique(df$year))

df    <- df    %>% filter(!is.na(year) & !(year %in% v_years_re))
df_su <- df_su %>% filter(!is.na(year) & !(year %in% v_years_re))
df_gr <- df_gr %>% filter(!is.na(year) & !(year %in% v_years_re))
df_fl <- df_fl %>% filter(!is.na(year) & !(year %in% v_years_re))
df_fr <- df_fr %>% filter(!is.na(year) & !(year %in% v_years_re))
df_re <- df_re %>% filter(!is.na(year) & !(year %in% v_years_re))

# Years analysed
v_years      <- sort(unique(df$year))


# Survival model ---------------------------------------------------------------
# GLMM; binomial
# Intercept model
mod_su_0 <- glmer(
  survives ~ 1 + (1 | year),
  data = df_su, family = binomial) 
# Linear size
mod_su_1 <- glmer(
  survives ~ logsize_t0 + (1 | year),
  data = df_su, family = binomial) 
# Quadratic size
mod_su_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_su, family = binomial)  
# Cubic size
mod_su_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
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


# Survival plot ----------------------------------------------------------------
su_yrs     <- data.frame(year = df_su$year %>% unique %>% sort)
su_bin_yrs <- lapply(1:nrow(su_yrs), df_binned_prop_year, df, 15,
                       logsize_t0, survives, su_yrs)
su_bin_yrs <- Filter(function(df) nrow(df) > 0, su_bin_yrs)


v <- rep(NA,length(su_bin_yrs))
for (ii in 1:length(su_bin_yrs)) {
  v[ii] <- su_bin_yrs[[ii]][1,'year']
}

surv_yr_plots <- function(i) {
  
  surv_temp    <- as.data.frame(su_bin_yrs[[i]])
  x_temp       <- seq(min(surv_temp$logsize_t0, na.rm = TRUE), 
                      max(surv_temp$logsize_t0, na.rm = TRUE), 
                      length.out = 100)
  
  linear_predictor <- mod_su_ranef[i, 1] 
  
  if (ncol(mod_su_ranef) >= 2) {
    linear_predictor <- linear_predictor + mod_su_ranef[i, 2] * x_temp}
  if (ncol(mod_su_ranef) >= 3) {
    linear_predictor <- linear_predictor + mod_su_ranef[i, 3] * x_temp^2}
  if (ncol(mod_su_ranef) >= 4) {
    linear_predictor <- linear_predictor + mod_su_ranef[i, 4] * x_temp^3}
  pred_temp <- boot::inv.logit(linear_predictor)
  if (ncol(mod_su_ranef) == 2) {line_color <- 'red'
  } else if (ncol(mod_su_ranef) == 3) {line_color <- 'green'
  } else if (ncol(mod_su_ranef) == 4) {line_color <- 'blue'
  } else {line_color <- 'black'  # Default color 
  } 
  
  pred_temp_df <- data.frame(logsize_t0 = x_temp, survives = pred_temp)
  temp_plot <- surv_temp %>% 
    ggplot() +
    geom_point(aes(x = logsize_t0, y = survives), size = 0.5) +
    geom_line(data = pred_temp_df, 
              aes(x = logsize_t0, y = survives), color = line_color, lwd = 1) +
    labs(title = v[i],
         x = expression('log(size)'[t0]),
         y = expression('Survival probability '[t1])) +
    ylim( 0, 1 ) +
    theme_bw() +
    theme(text         = element_text(size = 5),
          axis.title.y = element_text(
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(  
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text(
            margin  = margin(t = 1, r = 0, b = 0, l = 0)),
          axis.text.y  = element_text(
            margin  = margin(t = 0, r = 1, b = 0, l = 0)),
          plot.title   = element_text(
            margin  = margin(t = 2, r = 0, b = 1, l = 0), hjust  = 0.5),
          plot.margin  = margin(t = 0, r = 0, b = 0, l = 5))  
  
  return(temp_plot)
}

surv_yrs   <- lapply(1:length(su_bin_yrs), surv_yr_plots)

fig_surv_years <- wrap_plots(surv_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Survival year specific",  # Set the grand title
    subtitle = v_ggp_suffix,  # Set the grand subtitle
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 8)))


# Growth model -----------------------------------------------------------------
# Intercept model 
mod_gr_0 <- lmer(
  logsize_t1 ~ 1 + (logsize_t0 | year), 
  data = df_gr)
# Linear size
mod_gr_1 <- lmer(
  logsize_t1 ~ logsize_t0 + (logsize_t0 | year), 
  data = df_gr)
# Quadratic size
mod_gr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + (logsize_t0 | year),
  data = df_gr)
# Cubic size
mod_gr_3 <- lmer(
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



# Growth plot ------------------------------------------------------------------
grow_yr_plots <- function(i){
  temp_f <- function(x) {
    linear_predictor <- mod_gr_ranef[which(rownames(mod_gr_ranef) == i), 1]
    if (ncol(mod_gr_ranef) >= 2) {
      linear_predictor <- linear_predictor + 
        mod_gr_ranef[which(rownames(mod_gr_ranef) == i), 2] * x}
    if (ncol(mod_gr_ranef) >= 3) {
      linear_predictor <- linear_predictor + 
        mod_gr_ranef[which(rownames(mod_gr_ranef) == i), 3] * x^2}
    if (ncol(mod_gr_ranef) >= 4) {
      linear_predictor <- linear_predictor + 
        mod_gr_ranef[which(rownames(mod_gr_ranef) == i), 4] * x^3}
    return(linear_predictor)
  }
  if (ncol(mod_gr_ranef) == 2) {line_color <- 'red' } 
  else if (ncol(mod_gr_ranef) == 3) {line_color <- 'green'} 
  else if (ncol(mod_gr_ranef) == 4) {line_color <- 'blue'} 
  else {line_color <- 'black'}
  temp_plot <- df_gr %>% 
    filter(year == i) %>% 
    ggplot() +
    geom_jitter(aes(x = logsize_t0, y = logsize_t1), 
                size = 0.5, alpha = 0.5, width = .05, height = .1) +
    geom_function(fun = temp_f, color = line_color, lwd = 1) +
    geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2) +
    labs(title = i,
         x = expression('log(size) '[ t0]),
         y = expression('log(size) '[ t1])) +
    theme_bw() +
    theme(text         = element_text(size = 5),
          axis.title.y = element_text(
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text( 
            margin = margin(t = 1, r = 0, b = 0, l = 0)),
          axis.text.y  = element_text( 
            margin = margin(t = 0, r = 1, b = 0, l = 0)),
          plot.title   = element_text( 
            margin = margin(t = 2, r = 0, b = 1, l = 0), hjust  = 0.5 ), 
          plot.margin  = margin(t = 0, r = 2, b = 0, l = 0))
  return(temp_plot)
}

grow_yrs     <- lapply(sort(unique(df_gr$year)), grow_yr_plots)
fig_grow_years <- wrap_plots(grow_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Growth - year specific",  
    subtitle = v_ggp_suffix,  
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9)))


# Growth variance year specific ------------------------------------------------
x      <- fitted(mod_gr_best)
y      <- resid( mod_gr_best)^2
gr_var <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0))


# Flower model -----------------------------------------------------------------
# GLMM; binomial
# Intercept model
mod_fl_0 <- glmer(
  flower ~ 1 + (1 | year),
  data = df_fl, family = binomial) 
# Linear size term
mod_fl_1 <- glmer(
  flower ~ logsize_t0 + (1 | year),
  data = df_fl, family = binomial) 
# Quadratic size term
mod_fl_2 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + (1 | year),
  data = df_fl, family = binomial)  
# Cubic size term
mod_fl_3 <- glmer(
  flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year),
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


# Flower plot ------------------------------------------------------------------
fl_yrs     <- data.frame(year = df_fl$year %>% unique %>% sort)
fl_bin_yrs <- lapply(1:nrow(fl_yrs), df_binned_prop_year, df_fl, 15,
                     logsize_t0, flower, fl_yrs)
fl_bin_yrs <- Filter(function(df) nrow(df) > 0, fl_bin_yrs)


v <- rep(NA,length(fl_bin_yrs))
for (ii in 1:length(fl_bin_yrs)) {
  v[ii] <- fl_bin_yrs[[ii]][1,'year']
}

fl_yr_plots <- function(i) {
  
  fl_temp    <- as.data.frame(fl_bin_yrs[[i]])
  x_temp       <- seq(min(fl_temp$logsize_t0, na.rm = TRUE), 
                      max(fl_temp$logsize_t0, na.rm = TRUE), 
                      length.out = 100)
  
  linear_predictor <- mod_fl_ranef[i, 1] 
  
  if (ncol(mod_fl_ranef) >= 2) {
    linear_predictor <- linear_predictor + mod_fl_ranef[i, 2] * x_temp}
  if (ncol(mod_fl_ranef) >= 3) {
    linear_predictor <- linear_predictor + mod_fl_ranef[i, 3] * x_temp^2}
  if (ncol(mod_fl_ranef) >= 4) {
    linear_predictor <- linear_predictor + mod_fl_ranef[i, 4] * x_temp^3}
  pred_temp <- boot::inv.logit(linear_predictor)
  if (ncol(mod_fl_ranef) == 2) {line_color <- 'red'
  } else if (ncol(mod_fl_ranef) == 3) {line_color <- 'green'
  } else if (ncol(mod_fl_ranef) == 4) {line_color <- 'blue'
  } else {line_color <- 'black'  # Default color 
  } 
  
  pred_temp_df <- data.frame(logsize_t0 = x_temp, flower = pred_temp)
  temp_plot <- fl_temp %>% 
    ggplot() +
    geom_point(aes(x = logsize_t0, y = flower), size = 0.5) +
    geom_line(data = pred_temp_df, 
              aes(x = logsize_t0, y = flower), color = line_color, lwd = 1) +
    labs(title = v[i],
         x = expression('log(size)'[t0]),
         y = expression('flival probability '[t1])) +
    ylim( 0, 1 ) +
    theme_bw() +
    theme(text         = element_text(size = 5),
          axis.title.y = element_text(
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(  
            margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text(
            margin  = margin(t = 1, r = 0, b = 0, l = 0)),
          axis.text.y  = element_text(
            margin  = margin(t = 0, r = 1, b = 0, l = 0)),
          plot.title   = element_text(
            margin  = margin(t = 2, r = 0, b = 1, l = 0), hjust  = 0.5),
          plot.margin  = margin(t = 0, r = 0, b = 0, l = 5))  
  
  return(temp_plot)
}

fl_yrs   <- lapply(1:length(fl_bin_yrs), fl_yr_plots)

fig_fl_years <- wrap_plots(fl_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Flowering year specific",  # Set the grand title
    subtitle = v_ggp_suffix,  # Set the grand subtitle
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 8)))


# Fruit model ------------------------------------------------------------------
mod_fr_0 <- glmer.nb(
  fruit ~ 1 + (1 | year), data = df_fr)
mod_fr_1 <- glmer.nb(
  fruit ~ logsize_t0 + (1 | year), data = df_fr)
mod_fr_2 <- glmer.nb(
  fruit ~ logsize_t0 + logsize_t0_2 + (1 | year), data = df_fr)

# Model failed to converge with max|grad| = 0.0024376
# mod_fr_3 <- glmer.nb(
#   fruit ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (1 | year), data = df_fr)

# Compare models using AIC
mods_fr      <- list(mod_fr_0, mod_fr_1, mod_fr_2) #, mod_fr_3
mods_fr_dAIC <- bbmle::AICctab(mods_fr, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fr_i_sort <- order(mods_fr_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_fr) == 0) {
  mod_fr_i_best <- mods_fr_i_sort[1]
  v_mod_fr_i    <- mod_fr_i_best - 1 
} else {
  mod_fr_i_best <- v_mod_set_fr +1
  v_mod_fr_i    <- v_mod_set_fr
}

mod_fr_best       <- mods_fr[[mod_fr_i_best]]
mod_fr_ranef      <- coef(mod_fr_best)$year
mod_fr_ranef$year <- rownames(mod_fr_ranef)


# Fruit plot -------------------------------------------------------------------
# Plotting function for each year
fruit_yr_plot <- function(yr) {
  df_fr_sub <- df_fr %>% filter(year == yr)
  coefs     <- mod_fr_ranef %>% filter(year == yr)
  
  intercept <- coefs[1, "(Intercept)"]
  slope     <- if ("logsize_t0"   %in% names(coefs)) coefs[1, "logsize_t0"]   else 0
  quad      <- if ("logsize_t0_2" %in% names(coefs)) coefs[1, "logsize_t0_2"] else 0
  cubic     <- if ("logsize_t0_3" %in% names(coefs)) coefs[1, "logsize_t0_3"] else 0
  
  pred_fun <- function(x) {
    eta <- intercept + slope * x + quad * x^2 + cubic * x^3
    mu <- exp(eta)  # Negative binomial with log link
    return(mu)
  }
  
  ggplot(df_fr_sub, aes(x = logsize_t0, y = fruit)) +
    geom_jitter(width = 0.05, height = 0.2, alpha = 0.5, size = 0.6) +
    stat_function(fun = pred_fun, color = if_else(
      cubic > 0, 'blue', if_else(quad > 0, 'green', 'red')), linewidth = 1) +
    labs(title = paste("Year:", yr),
         x = expression(log(size)[t[0]]),
         y = "Fruit count") +
    theme_bw(base_size = 8)
}

# Generate one plot per year
fruit_plots <- lapply(sort(unique(df_fr$year)), fruit_yr_plot)

# Combine and annotate
fig_fruit_years <- wrap_plots(fruit_plots) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = "Fruit Production - year specific",
    theme = theme(plot.title = element_text(size = 14, face = "bold")))


# Fruit to recruit -------------------------------------------------------------
# Calculate fruit produced each year (optionally shift year if needed)
fruit_by_year <- df %>%
  filter(!is.na(fruit)) %>%
  group_by(year) %>%
  summarise(total_fruit_t0 = sum(fruit, na.rm = TRUE), .groups = 'drop')

# Calculate number of recruits each year
recruits_by_year <- df %>%
  filter(recruit == 1) %>%
  group_by(year) %>%
  # I am modeling the transition of fruits of this year to the recruits of next year 
  mutate(year = year - 1) %>% 
  summarise(n_recruits_t1 = n(), .groups = 'drop')

# Combine and calculate year-specific ratio
fruit_recruit_ratio_by_year <- recruits_by_year %>%
  left_join(fruit_by_year, by = "year") %>%
  mutate(repr_pc = n_recruits_t1 / total_fruit_t0)



# Recruitment model ------------------------------------------------------------
df_re_mod <- df_re %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
mod_rec <- glmer.nb(nr_quad ~ (1 | year), data = df_re_mod)

# predict the number of recruits per year per quad
df_re_mod <- df_re_mod %>% 
  mutate(pred_mod = predict(mod_rec, type = 'response'))

# sum up the observed and predicted number of recruits per year across all quads
rec_sums_df <- df_re_mod %>% 
  group_by(year) %>% 
  summarise(nr_quad  = sum(nr_quad),
            pred_mod = sum(pred_mod)) %>% 
  ungroup

# number of adults present in each year
indiv_yr <- df %>%
  filter(!is.na(survives)) %>% 
  count(year) %>% 
  rename(n_adults = n) %>% 
  mutate(year = year + 1)

# calculate per-capita recruitment rate
repr_pc_yr <- indiv_yr %>% 
  left_join( rec_sums_df ) %>%
  mutate(repr_percapita = pred_mod / n_adults,
         repr_pc_obs    = nr_quad / n_adults,
         year = year - 1 ) %>% 
  drop_na


# Recruitment plot -------------------------------------------------------------
fig_re <- repr_pc_yr %>% 
  ggplot() +
  geom_point( aes(x = repr_pc_obs,
                  y = repr_percapita)) +
  geom_abline(aes(intercept = 0, slope = 1),
              color = 'red', lwd = 2, alpha = 0.5) +
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = 'Observed per capita recruitment',
       y        = 'Predicted per capita recruitment') +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 8))


# Exporting parameter estimates ------------------------------------------------
# Survival
# Get the coefficients matrix
su_coef_matrix <- coef(mod_su_best)$year

# Initialize a list to store data frames
su_data_frames <- list(
  data.frame(coefficient = paste0('year_', rownames(su_coef_matrix)), 
             value = su_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0_', rownames(su_coef_matrix)), 
             value = su_coef_matrix[, 'logsize_t0'])
)

# Loop to create additional data frames if needed
for (i in 2:ncol(su_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(su_coef_matrix)) {
    su_data_frames[[length(su_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, rownames(su_coef_matrix)),
      value = su_coef_matrix[, column_name]
    )
  }
}

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
             value = gr_coef_matrix[, 'logsize_t0'])
)

# Loop to create additional logsize data frames
for (i in 2:ncol(gr_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(gr_coef_matrix)) {
    gr_data_frames[[length(gr_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, '_', rownames(gr_coef_matrix)),
      value = gr_coef_matrix[, column_name]
    )
  }
}

# Combine all data frames using Reduce and mutate the coefficient column
gr_out_yr <- Reduce(function(...) rbind(...), gr_data_frames) %>%
  mutate(coefficient = as.character(coefficient))

 
# Flower
# Get the coefficients matrix
fl_coef_matrix <- coef(mod_fl_best)$year

# Initialize a list to store data frames
fl_data_frames <- list(
  data.frame(coefficient = paste0('year_', rownames(fl_coef_matrix)), 
             value = fl_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0_', rownames(fl_coef_matrix)), 
             value = fl_coef_matrix[, 'logsize_t0'])
)

# Loop to create additional data frames if needed
for (i in 2:ncol(fl_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(fl_coef_matrix)) {
    fl_data_frames[[length(fl_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, rownames(fl_coef_matrix)),
      value = fl_coef_matrix[, column_name]
    )
  }
}

# Combine data frames and mutate the coefficient column
fl_out_yr <- Reduce(rbind, fl_data_frames) %>%
  mutate(coefficient = as.character(coefficient))



# Fruit
# Get the coefficients matrix
fr_coef_matrix <- coef(mod_fr_best)$year

# Initialize a list to store data frames
fr_data_frames <- list(
  data.frame(coefficient = paste0('year_', rownames(fr_coef_matrix)), 
             value = fr_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0_', rownames(fr_coef_matrix)), 
             value = fr_coef_matrix[, 'logsize_t0'])
)

# Loop to create additional data frames if needed
for (i in 2:ncol(fr_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(fr_coef_matrix)) {
    fr_data_frames[[length(fr_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, rownames(fr_coef_matrix)),
      value = fr_coef_matrix[, column_name]
    )
  }
}

# Combine data frames and mutate the coefficient column
fr_out_yr <- Reduce(rbind, fr_data_frames) %>%
  mutate(coefficient = as.character(coefficient))



# Fruit to recruit
ftr_out_yr <- fruit_recruit_ratio_by_year %>%
  mutate(coefficient = paste0('year_', year),
         value = repr_pc) %>% 
  dplyr::select(coefficient, value)




# Recruitment
rc_pc <- data.frame(coefficient = paste0('rec_pc_',repr_pc_yr$year),
                    value = repr_pc_yr$repr_percapita)

rec_size          <- df %>% subset(recruit == 1)

rc_sz <- data.frame(coefficient = c('rec_siz', 'rec_sd'),
                    value = c(mean(rec_size$logsize_t0, na.rm = TRUE),
                              sd(  rec_size$logsize_t0, na.rm = TRUE)))
mean()

recr_out_yr <- Reduce(function(...) rbind(...), list(rc_pc, rc_sz)) %>%
  mutate(coefficient = as.character(coefficient))



# df constant parameters, fixed effects estimates, and mean parameter estimates
constants <- data.frame(coefficient = c('recr_sz',
                                        'recr_sd',
                                        'a',
                                        'b',
                                        'L',
                                        'U',
                                        'mat_siz'),
                        value = c(mean(df_re_mod$logsize_t0),
                                  sd(  df_re_mod$logsize_t0),
                                  as.numeric(coef(gr_var)[1]),
                                  as.numeric(coef(gr_var)[2]),
                                  grow_df$logsize_t0 %>% min,
                                  grow_df$logsize_t0 %>% max,
                                  200))

# Create the data frame dynamically based on the number of fixed effects
surv_fe <- data.frame(
  coefficient = paste0('surv_b', 0:(length(fixef(mod_su_best)) - 1)),
  value       = fixef(mod_su_best))

grow_fe <- data.frame(
  coefficient = paste0('grow_b', 0:(length(fixef(mod_gr_best)) - 1)),
  value       = fixef(mod_gr_best))

rec_fe  <- data.frame(coefficient = 'fecu_b0',
                      value       = mean(repr_pc_yr$repr_percapita))

pars_cons <- Reduce(function(...) rbind(...), 
                    list(surv_fe, grow_fe, rec_fe, constants)) %>%
  mutate(coefficient = as.character(coefficient))

rownames(pars_cons) <- 1:nrow(pars_cons)

pars_cons_wide <- as.list(pivot_wider(pars_cons, names_from = 'coefficient', 
                                      values_from = 'value'))




# DF varying parameters
# Function to create coefficient data frames dynamically
create_coef_df <- function(model, prefix) {
  coef_matrix <- coef(model)$year
  lapply(0:(ncol(coef_matrix) - 1), function(i) {
    column_name <- if (i == 0) '(Intercept)' 
    else paste0('logsize_t0', ifelse(i == 1, '', paste0('_', i)))
    data.frame(coefficient = paste0(prefix, i, '_', rownames(coef_matrix)),
               value = coef_matrix[, column_name])
  })
}

# Create data frames for survival and growth models
su_data_frames <- create_coef_df(mod_su_best, 'surv_b')
gr_data_frames <- create_coef_df(mod_gr_best, 'grow_b')

# Create the fecundity data frame
fecu_b0 <- data.frame(coefficient = paste0('fecu_b0_', repr_pc_yr$year),
                      value = repr_pc_yr$repr_percapita)

# Combine all data frames into one
pars_var <- Reduce(rbind, c(su_data_frames, gr_data_frames, list(fecu_b0)))

pars_var_wide <- as.list(pivot_wider(pars_var, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value') )



# ----------------------
# 1. Fruit per plot per year
fruit_by_plot_year <- df %>%
  filter(!is.na(fruit)) %>%
  group_by(year, quad_id) %>%
  summarise(total_fruit_t0 = sum(fruit, na.rm = TRUE), .groups = 'drop')

# 2. Recruits per plot per year (shifted back to match fruit of previous year)
recruits_by_plot_year <- df %>%
  filter(recruit == 1) %>%
  mutate(year = year - 1) %>%  # shift recruits back one year
  group_by(year, quad_id) %>%
  summarise(n_recruits_t1 = n(), .groups = 'drop')

# 3. Combine fruit and recruit data at plot-year level
plot_level_ratios <- recruits_by_plot_year %>%
  left_join(fruit_by_plot_year, by = c("year", "quad_id")) %>%
  mutate(total_fruit_t0 = ifelse(
    total_fruit_t0 < n_recruits_t1, n_recruits_t1, total_fruit_t0),
    repr_pc = n_recruits_t1 / total_fruit_t0)

# 4. Mean ratio per year across all plots
yearly_mean_ratios <- plot_level_ratios %>%
  group_by(year) %>%
  summarise(mean_repr_pc = mean(repr_pc, na.rm = TRUE),
            sd_repr_pc   = sd(repr_pc, na.rm = TRUE),
            n_plots      = n(),
            .groups = 'drop')

# Output
print(yearly_mean_ratios)



