# IPM year specific pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.11.11


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)


# Packages ---------------------------------------------------------------------
# load packages
source('helper_functions/load_packages.R' )
load_packages( 
  tidyverse, patchwork, skimr, lme4, bbmle, ipmr, readxl, binom)


# Specification ----------------------------------------------------------------
# Create a unique species abbreviation for file naming
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))
# Define script prefix
v_script_prefix <- str_c(
  str_extract(v_author_year, "^[^_]+"),
  str_sub(str_extract(v_author_year, "_\\d+$"), -2, -1))
# Define prefix for two of the same author and year
if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0("^", v_author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

# Define suffix for plot outputs
v_suffix     <- c("")
if (length(v_years_re)       > 0) {v_suffix <- paste0(v_suffix, "_yr1")}
if (length(v_size_threshold) > 0) {v_suffix <- paste0(v_suffix, "_st1")}
if (length(v_mod_set_gr)     > 0) {v_suffix <- paste0(v_suffix, "_gr1")}
if (length(v_mod_set_su)     > 0) {v_suffix <- paste0(v_suffix, "_su1")}

# Define graph subtitle
v_ggp_suffix <- paste(
  paste0(toupper(substr(v_script_prefix, 1, 1)), 
         substr(v_script_prefix, 2, nchar(v_script_prefix))), '/', 
  v_species, 
  '\n Size threshold:', 
  ifelse(is.null(v_size_threshold), !is.null(v_size_threshold), v_size_threshold),
  '\n Model complexity altered in growth / survival:', 
  ifelse(is.null(v_mod_set_gr), !is.null(v_mod_set_gr), v_mod_set_gr), '/', 
  ifelse(is.null(v_mod_set_su), !is.null(v_mod_set_su), v_mod_set_su),
  '\n Years removed:',
  ifelse(is.null(v_years_re), !is.null(v_years_re), paste(v_years_re, collapse = ", ")))


# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_author_year, '_', v_region_abb))
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


# Save the suffix --------------------------------------------------------------
write.csv(v_suffix, file.path(dir_data, 'v_suffix.csv'), row.names = F)


# Ipm mean and plant tracker if they not already exists ------------------------
if (!file.exists(
  paste0(dir_data, '/', v_script_prefix, '_', v_sp_abb, '_data_df.csv'))) {
  source(paste0(dir_R, '/', v_script_prefix, '_', v_sp_abb, '_ipm_mean.R'))
}


# Data -------------------------------------------------------------------------
df      <- read.csv(
  paste0(dir_data, '/', v_script_prefix, '_', v_sp_abb, '_data_df.csv'))
grow_df <- read.csv(
  paste0(dir_data, '/', v_script_prefix, '_', v_sp_abb, '_growth_df.csv'))
surv_df <- read.csv(
  paste0(dir_data, '/', v_script_prefix, '_', v_sp_abb, '_survival_df.csv'))
recr_df <- read.csv(
  paste0(dir_data, '/', v_script_prefix, '_', v_sp_abb, '_recruitment_df.csv'))

# Implement size threshold
if (length(v_size_threshold) > 0) {
  df <- df %>% 
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0)) %>% 
    filter(logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

df_long <- pivot_longer(
  df, cols = c(logsize_t0, logsize_t1 ), 
  names_to = 'size', values_to = 'size_value' ) %>% 
  select(c(size_t0, year, size,  size_value)) %>% 
  mutate(size = as.factor(size),
         year_fac = as.factor(year))

size_labs        <- c('at time t0', 'at time t1')
names(size_labs) <- c('logsize_t0', 'logsize_t1')

# Histogram of sizes, t0 and t1
g_hist_logsizes_years <- df_long %>% 
  ggplot(aes(x = size_value)) +
  geom_histogram(binwidth = 1) +
  facet_grid(year_fac ~ size, 
             scales = 'free_y',
             labeller = labeller(size = size_labs)) +
  labs(title    = 'Histogram',
       subtitle = v_ggp_suffix,
       x        = 'log(size)',
       y        = 'Frequency') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 5)) +
  theme(plot.subtitle = element_text(size = 8))

ggsave(paste0(dir_result, '/3.1_years_hist_logsizes_years', v_suffix, '.png'),
       plot = g_hist_logsizes_years, width = 6, height = 12, dpi = 150)


# Survival 
# function to plot your survival data 'binned' per year (instead of 'jittered')
source('helper_functions/plot_binned_prop_year.R')

surv_yrs       <- data.frame(year = surv_df$year %>% unique %>% sort)
surv_bin_yrs   <- lapply(1:nrow(surv_yrs), df_binned_prop_year, df, 15, 
                         logsize_t0, survives, surv_yrs)
surv_bin_yrs   <- Filter(function(df) nrow(df) > 0, surv_bin_yrs)

surv_yr_pan_df <- bind_rows(surv_bin_yrs) %>% 
  mutate(transition = paste(
    paste0(year), substr(paste0(year + 1), 3, 4), sep = '-')) %>% 
  mutate(year = as.integer(year - surv_yrs[1,]))

g_survival <- ggplot(
  data = surv_yr_pan_df, aes(x = logsize_t0, y = survives)) +
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr), 
                width = 0.25, linewidth = 0.1) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  # split in panels
  facet_wrap(.~ transition, ncol = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = 8), title = element_text(size = 10),
        strip.text.y = element_text(
          size = 5, margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.text.x = element_text(
          size = 5, margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
        panel.spacing         = unit( '0.5', unit = 'mm' ),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Survival to time t1')) 

ggsave(paste0(dir_result, '/3.2_years_survival', v_suffix, '.png'), 
       plot = g_survival,
       width = 4, height = 9, dpi = 150)


# Growth
grow_yr_pan_df <- grow_df %>%
  mutate(transition = paste(
    paste0(year), substr(
      paste0(year + 1), nchar(paste0(year + 1)) - 1, 
      nchar(paste0(year + 1))), sep = '-')) %>% 
  mutate(year       = as.integer(year - min(year)))

g_growth <- ggplot(
  data = grow_yr_pan_df, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  # split in panels
  facet_wrap(.~ transition, ncol = 4) +
  theme_bw() +
  theme(axis.text    = element_text(size = 8),
        title        = element_text(size = 10),
        strip.text.y = element_text(size = 8, 
                                    margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.text.x = element_text(size = 8, 
                                    margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.switch.pad.wrap = unit('0.5', unit = 'mm'),
        panel.spacing         = unit('0.5', unit = 'mm'),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        =  expression('log(size) '[t0]),
       y        = expression('log(size) '[t1]))
  
ggsave(paste0(dir_result, '/3.3_years_growth', v_suffix, '.png'),
       plot = g_growth,
       width = 4, height = 9, dpi = 150)


# Recruits
indiv_qd <- surv_df %>%
  group_by(quad) %>%
  count(year) %>% 
  rename(n_adults = n) %>% 
  mutate(year = year + 1)

repr_yr <- indiv_qd %>% 
  left_join(recr_df) %>%
  mutate(repr_pc = nr_quad / n_adults) %>% 
  mutate(year = year - 1) %>% 
  drop_na

g_recruits <- repr_yr %>% 
  filter(nr_quad  != max(repr_yr$nr_quad)) %>% 
  filter(n_adults != max(repr_yr$n_adults)) %>% 
  ggplot(aes(x = n_adults, y = nr_quad ) ) +
  geom_point(alpha = 1, pch = 16, size = 1, color = 'red') +
  facet_wrap(.~ year, ncol = 4) +
  theme_bw() +
  theme(axis.text    = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title        = element_text(size = 10),
        strip.text.y = element_text(size = 8, 
                                    margin = margin(0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.text.x = element_text(size = 8, 
                                    margin = margin(0.5, 0.5, 0.5, 0.5, 'mm')),
        strip.switch.pad.wrap = unit('0.5', unit = 'mm'),
        panel.spacing         = unit('0.5', unit = 'mm'),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Recruits',
       subtitle = v_ggp_suffix,
       x        = expression('Number of adults '[ t0]),
       y        = expression('Number of recruits '[ t1]))

ggsave(paste0(dir_result, '/3.4_v_years_recruits', v_suffix, '.png'), 
       plot = g_recruits,
       width = 4, height = 9, dpi = 150)

# Recruitment size
rec_size          <- df %>% subset(recruit == 1)
rec_size$year_fac <- as.factor(rec_size$year)

g_recruitment_size <- rec_size %>% 
  ggplot(aes(x = logsize_t0)) +
  geom_histogram() +
  facet_wrap(year_fac ~ ., scales = 'free_y', 
             ncol = 4) +
  labs(title    = 'Recruitment size',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = 'Frequency') +
  theme(plot.subtitle = element_text(size = 8))

ggsave(paste0(dir_result, '/3.4_v_years_recruitment_size', v_suffix, '.png'), 
       plot = g_recruitment_size,
       width = 4, height = 9, dpi = 150)


# Removing year with too few data ----------------------------------------------
df      <- df      %>% filter(!is.na(year) & !(year %in% v_years_re))
surv_df <- surv_df %>% filter(!is.na(year) & !(year %in% v_years_re))
grow_df <- grow_df %>% filter(!is.na(year) & !(year %in% v_years_re))
recr_df <- recr_df %>% filter(!is.na(year) & !(year %in% v_years_re))
surv_yrs     <- data.frame(year = surv_df$year %>% unique %>% sort)
surv_bin_yrs <- lapply(1:nrow(surv_yrs), df_binned_prop_year, df, 15,
                       logsize_t0, survives, surv_yrs)
surv_bin_yrs <- Filter(function(df) nrow(df) > 0, surv_bin_yrs)

years_v      <- sort(unique(df$year))


# Fitting vital rate models with the random effect of year ---------------------
# Survival year specific -------------------------------------------------------
su_mod_yr_0 <- glmer(
  survives ~ 1 + (1 | year),
  data = surv_df, family = binomial)
su_mod_yr   <- glmer(
  survives ~ logsize_t0 + (1 | year), 
  data = surv_df, family = binomial )
su_mod_yr_2 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 + (1 | year), 
  data = surv_df, family = binomial)
su_mod_yr_3 <- glmer(
  survives ~ logsize_t0 + logsize_t0_2 +  logsize_t0_3 + (1 | year), 
  data = surv_df, family = binomial)
su_mods     <- list(su_mod_yr_0, su_mod_yr, su_mod_yr_2, su_mod_yr_3)

# Assign the best model to the variable
su_dAIC_values    <- AICtab(su_mods, weights = T, sort = F)$dAIC
su_sorted_indices <- order(su_dAIC_values)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  su_mod_ys_index_bestfit <- su_sorted_indices[1]
  v_mod_su_index          <- su_mod_ys_index_bestfit - 1 
} else {
  su_mod_ys_index_bestfit <- v_mod_set_su +1
  v_mod_su_index          <- v_mod_set_su
}

# Specify the model
su_mod_yr_bestfit <- su_mods[[su_mod_ys_index_bestfit]]
su_ranef          <- data.frame(coef(su_mod_yr_bestfit)[1])


v <- rep(NA,length(surv_bin_yrs))
for (ii in 1:length(surv_bin_yrs)) {
  v[ii] <- surv_bin_yrs[[ii]][1,'year']
}

surv_yr_plots <- function(i) {
  
  surv_temp    <- as.data.frame(surv_bin_yrs[[i]])
  x_temp       <- seq(min(surv_temp$logsize_t0, na.rm = TRUE), 
                      max(surv_temp$logsize_t0, na.rm = TRUE), 
                      length.out = 100)
  
  linear_predictor <- su_ranef[i, 1] 
  
  if (ncol(su_ranef) >= 2) {
    linear_predictor <- linear_predictor + su_ranef[i, 2] * x_temp}
  if (ncol(su_ranef) >= 3) {
    linear_predictor <- linear_predictor + su_ranef[i, 3] * x_temp^2}
  if (ncol(su_ranef) >= 4) {
    linear_predictor <- linear_predictor + su_ranef[i, 4] * x_temp^3}
  pred_temp <- boot::inv.logit(linear_predictor)
  if (ncol(su_ranef) == 2) {line_color <- 'red'
  } else if (ncol(su_ranef) == 3) {line_color <- 'green'
  } else if (ncol(su_ranef) == 4) {line_color <- 'blue'
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

surv_yrs   <- lapply(1:length(surv_bin_yrs), surv_yr_plots)
g_surv_years <- wrap_plots(surv_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Survival year specific",  # Set the grand title
    subtitle = v_ggp_suffix,  # Set the grand subtitle
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 8)))

ggsave(paste0(
  dir_result, '/4.1_years_surv_logsize', v_suffix, '.png'), 
       plot  = g_surv_years,
       width = 4, height = 9, dpi = 150)


# Growth year specific ---------------------------------------------------------
gr_mod_yr_0 <- lmer(
  logsize_t1 ~ 1 + (logsize_t0 | year), 
  data = grow_df)
gr_mod_yr   <- lmer(
  logsize_t1 ~ logsize_t0 + (logsize_t0 | year), 
       data = grow_df)
gr_mod_yr_2 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + (logsize_t0 | year), 
       data = grow_df)
# Model failed to converge with max|grad| = 0.0297295
gr_mod_yr_3 <- lmer(
  logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + (logsize_t0 | year), 
       data = grow_df)

gr_mods <- list(gr_mod_yr_0, gr_mod_yr, gr_mod_yr_2, gr_mod_yr_3)
# Assign the best model
gr_dAIC_values <- AICtab(gr_mods, weights = T, sort = F)$dAIC
# Get the sorted indices of dAIC values
gr_sorted_indices <- order(gr_dAIC_values)

# Establish the index of model complexity
if (length(v_mod_set_gr) == 0) {
  gr_mod_ys_index_bestfit <- gr_sorted_indices[1]
  v_mod_gr_index          <- gr_mod_ys_index_bestfit - 1 
} else {
  gr_mod_ys_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index          <- v_mod_set_gr
}

# Specify the model
gr_mod_yr_bestfit <- gr_mods[[gr_mod_ys_index_bestfit]]
gr_ranef          <- data.frame(coef(gr_mod_yr_bestfit)[1])


grow_yr_plots <- function(i){
  temp_f <- function(x) {
    linear_predictor <- gr_ranef[which(rownames(gr_ranef) == i), 1]
    if (ncol(gr_ranef) >= 2) {
      linear_predictor <- linear_predictor + 
        gr_ranef[which(rownames(gr_ranef) == i), 2] * x}
    if (ncol(gr_ranef) >= 3) {
      linear_predictor <- linear_predictor + 
        gr_ranef[which(rownames(gr_ranef) == i), 3] * x^2}
    if (ncol(gr_ranef) >= 4) {
      linear_predictor <- linear_predictor + 
        gr_ranef[which(rownames(gr_ranef) == i), 4] * x^3}
    return(linear_predictor)
  }
  if (ncol(gr_ranef) == 2) {line_color <- 'red' } 
  else if (ncol(gr_ranef) == 3) {line_color <- 'green'} 
  else if (ncol(gr_ranef) == 4) {line_color <- 'blue'} 
  else {line_color <- 'black'}
  temp_plot <- grow_df %>% 
    filter(year == i) %>% 
    ggplot() +
    geom_point(aes(x = logsize_t0, y = logsize_t1), size = 0.5, alpha = 0.5) +
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

grow_yrs     <- lapply(sort(unique(grow_df$year)), grow_yr_plots)
g_grow_years <- wrap_plots(grow_yrs) + 
  plot_layout(ncol = 4) + 
  plot_annotation(
    title = "Growth - year specific",  
    subtitle = v_ggp_suffix,  
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9)))

ggsave(paste0(
  dir_result, '/4.2_years_growth_logsize', v_suffix, '.png'), 
       plot = g_grow_years,
       width = 4, height = 9, dpi = 150)


# Growth variance year specific ------------------------------------------------
x      <- fitted(gr_mod_yr_bestfit)
y      <- resid( gr_mod_yr_bestfit)^2
gr_var <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0))


# Recruitment model year specific ----------------------------------------------
recr_nona_nr_quad <- recr_df %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
rec_mod <- glmer.nb(nr_quad ~ (1 | year), data = recr_nona_nr_quad)

# predict the number of recruits per year per quad
recr_nona_nr_quad <- recr_nona_nr_quad %>% 
  mutate(pred_mod = predict(rec_mod, type = 'response')) 

# sum up the observed and predicted number of recruits per year across all quads
rec_sums_df <- recr_nona_nr_quad %>% 
  group_by(year) %>% 
  summarise(nr_quad  = sum(nr_quad),
            pred_mod = sum(pred_mod)) %>% 
  ungroup

# number of adults present in each year ----------------------------------------
indiv_yr <- surv_df %>%
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

# recruitment plot
g_recruitment <- repr_pc_yr %>% 
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

ggsave(paste0(dir_result, '/4.3_v_years_recruitment', v_suffix, '.png'), 
       plot = g_recruitment,
       width = 6, height = 4, dpi = 150)


# Exporting parameter estimates ------------------------------------------------
# Survival
# Get the coefficients matrix
su_coef_matrix <- coef(su_mod_yr_bestfit)$year

# Initialize a list to store data frames
su_data_frames <- list(
  data.frame(coefficient = paste0('year_', rownames(su_coef_matrix)), 
             value = su_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0', rownames(su_coef_matrix)), 
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
surv_out_yr <- Reduce(rbind, su_data_frames) %>%
  mutate(coefficient = as.character(coefficient))

write.csv(surv_out_yr, 
          paste0(dir_data, '/2.surv_pars', v_suffix, '.csv'), 
          row.names = F)


# Growth
# Get coefficients matrix
gr_coef_matrix <- coef(gr_mod_yr_bestfit)$year

# Initialize a list for year and logsize data frames, including var_fe
gr_data_frames <- list(
  data.frame(coefficient = names(coef(gr_var)), value = coef(gr_var)),
  data.frame(coefficient = paste0('year_', rownames(gr_coef_matrix)), 
             value = gr_coef_matrix[, '(Intercept)']),
  data.frame(coefficient = paste0('logsize_t0', rownames(gr_coef_matrix)), 
             value = gr_coef_matrix[, 'logsize_t0'])
)

# Loop to create additional logsize data frames
for (i in 2:ncol(gr_coef_matrix)) {
  column_name <- paste0('logsize_t0_', i - 1)
  if (column_name %in% colnames(gr_coef_matrix)) {
    gr_data_frames[[length(gr_data_frames) + 1]] <- data.frame(
      coefficient = paste0('logsize_t0_', i - 1, rownames(gr_coef_matrix)),
      value = gr_coef_matrix[, column_name]
    )
  }
}

# Combine all data frames using Reduce and mutate the coefficient column
grow_out_yr <- Reduce(function(...) rbind(...), gr_data_frames) %>%
  mutate(coefficient = as.character(coefficient))

write.csv(grow_out_yr, 
          paste0(dir_data, '/2.grow_pars', v_suffix, '.csv'), 
          row.names = F)


# Recruitment
rc_pc <- data.frame(coefficient = paste0('rec_pc_',repr_pc_yr$year),
                    value = repr_pc_yr$repr_percapita)

rc_sz <- data.frame(coefficient = c('rec_siz', 'rec_sd'),
                    value = c(mean(rec_size$logsize_t0),
                              sd(rec_size$logsize_t0)))

recr_out_yr <- Reduce(function(...) rbind(...), list(rc_pc, rc_sz)) %>%
  mutate(coefficient = as.character(coefficient))

write.csv(recr_out_yr, 
          paste0(dir_data, '/2.recr_pars', v_suffix, '.csv'), 
          row.names = F)


# df constant parameters, fixed effects estimates, and mean parameter estimates
constants <- data.frame(coefficient = c('recr_sz',
                                        'recr_sd',
                                        'a',
                                        'b',
                                        'L',
                                        'U',
                                        'mat_siz'),
                        value = c(mean(rec_size$logsize_t0),
                                  sd(  rec_size$logsize_t0),
                                  as.numeric(coef(gr_var)[1]),
                                  as.numeric(coef(gr_var)[2]),
                                  grow_df$logsize_t0 %>% min,
                                  grow_df$logsize_t0 %>% max,
                                  200))

# Create the data frame dynamically based on the number of fixed effects
surv_fe <- data.frame(
  coefficient = paste0('surv_b', 0:(length(fixef(su_mod_yr_bestfit)) - 1)),
  value       = fixef(su_mod_yr_bestfit))

grow_fe <- data.frame(
  coefficient = paste0('grow_b', 0:(length(fixef(gr_mod_yr_bestfit)) - 1)),
  value       = fixef(gr_mod_yr_bestfit))

rec_fe  <- data.frame(coefficient = 'fecu_b0',
                      value       = mean(repr_pc_yr$repr_percapita))

pars_cons <- Reduce(function(...) rbind(...), 
                    list(surv_fe, grow_fe, rec_fe, constants)) %>%
  mutate(coefficient = as.character(coefficient))

rownames(pars_cons) <- 1:nrow(pars_cons)

pars_cons_wide <- as.list(pivot_wider(pars_cons, names_from = 'coefficient', 
                                      values_from = 'value'))

write.csv(pars_cons_wide, 
          paste0(dir_data, '/2.pars_cons', v_suffix, '.csv'), 
          row.names = F)


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
su_data_frames <- create_coef_df(su_mod_yr_bestfit, 'surv_b')
gr_data_frames <- create_coef_df(gr_mod_yr_bestfit, 'grow_b')

# Create the fecundity data frame
fecu_b0 <- data.frame(coefficient = paste0('fecu_b0_', repr_pc_yr$year),
                      value = repr_pc_yr$repr_percapita)

# Combine all data frames into one
pars_var <- Reduce(rbind, c(su_data_frames, gr_data_frames, list(fecu_b0)))

pars_var_wide <- as.list(pivot_wider(pars_var, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value') )

write.csv(pars_var_wide, 
          paste0(dir_data, '/2.pars_var', v_suffix, '.csv'), 
          row.names = F)


# Building the year-specific IPMs from scratch ---------------------------------
# Functions
# Standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b * x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_params = gr_mod_ys_index_bestfit) {
  mean_value <- 0
  for (i in 0:num_params) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}

# Inverse logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

# Survival of x-sized individual to time t1
sx <- function(x, pars, num_params = su_mod_ys_index_bestfit) {
  survival_value <- pars$surv_b0
  for (i in 1:num_params) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^(i)
    }
  }
  return(inv_logit(survival_value))
}

# Transition of x-sized individual to y-sized individual at time t1
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Per-capita production of y-sized recruits
fy <- function(y, pars, h){
  n_recr  <- pars$fecu_b0
  recr_y  <- dnorm(y, pars$recr_sz, pars$recr_sd) * h
  recr_y  <- recr_y / sum(recr_y)
  f       <- n_recr * recr_y
  return(f)
}

# Kernel
kernel <- function(pars) {
  
  n      <- pars$mat_siz
  L      <- pars$L
  U      <- pars$U
  h      <- (U - L) / n
  b      <- L + c(0:n) * h
  y      <- 0.5 * (b[1:n] + b[2:( n + 1 )])
  
  Fmat   <- matrix(0, n, n)
  Fmat[] <- matrix(fy(y, pars, h), n, n)
  
  Smat   <- c()
  Smat   <- sx(y, pars)
  
  Gmat   <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  Tmat   <- matrix(0, n, n)
  
  for(i in 1:(n / 2)) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  for(i in (n / 2 + 1):n) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  k_yx <- Fmat + Tmat
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y))
  
}


# mean population growth rate
pars_mean <- pars_cons_wide

lambda_ipm <- function(i) {
  return(Re(eigen(kernel(i)$k_yx)$value[1]))
}

lam_mean <- lambda_ipm(pars_mean)
lam_mean


# population growth rates for each year
pars_yr <- vector(mode = 'list', length = length(years_v))
extr_value_list <- function(x, field) {
  return(as.numeric(x[paste0(field)] %>% unlist()))
}

prep_pars <- function(i, num_surv_params, num_grow_params) {
  yr_now <- years_v[i]
  
  # Initialize the parameters list with the required order
  pars_year <- list(
    surv_b0  = extr_value_list(pars_var_wide, paste('surv_b0', yr_now, sep = '_')),
    grow_b0  = extr_value_list(pars_var_wide, paste('grow_b0', yr_now, sep = '_')),
    a        = extr_value_list(pars_cons_wide, 'a'),
    b        = extr_value_list(pars_cons_wide, 'b'),
    fecu_b0  = extr_value_list(pars_var_wide, paste('fecu_b0', yr_now, sep = '_')),
    recr_sz   = extr_value_list(pars_cons_wide, 'recr_sz'),
    recr_sd   = extr_value_list(pars_cons_wide, 'recr_sd'),
    L        = extr_value_list(pars_cons_wide, 'L'),
    U        = extr_value_list(pars_cons_wide, 'U'),
    mat_siz  = 200
  )
  
  # Dynamically add survival parameters based on num_surv_params
  for (j in 1:num_surv_params) {
    param_name <- paste0('surv_b', j)
    value <- extr_value_list(pars_var_wide, paste(param_name, yr_now, sep = '_'))
    if (!is.null(value)) {
      # Insert after surv_b0
      pars_year <- append(pars_year, setNames(list(value), param_name), after = 1)
    }
  }
  
  # Dynamically add growth parameters based on num_grow_params
  for (j in num_grow_params:1) {
    param_name <- paste0('grow_b', j)
    value <- extr_value_list(pars_var_wide, paste(param_name, yr_now, sep = '_'))
    if (!is.null(value)) {
      # Insert immediately after grow_b1
      pos <- which(names(pars_year) == 'grow_b0') + 1
      pars_year <- append(pars_year, setNames(list(value), param_name), after = pos - 1)
    }
  }
  
  # Return the list with the dynamic parameters included
  return(pars_year)
}

pars_yr <- lapply(1:length(years_v), 
                  num_surv_params = v_mod_su_index, 
                  num_grow_params = v_mod_gr_index, 
                  prep_pars)


# Identify which years contain parameters with numeric(0)
contains_numeric0 <- sapply(pars_yr, function(regular_list) {
  any(sapply(regular_list, function(sublist) {
    identical(sublist, numeric(0))
  }))
})

which_contains_numeric0 <- which(contains_numeric0)
v_years_rm_sugg         <- years_v[which_contains_numeric0]
v_years_og              <- years_v
# CHECK -- Exclude these years ##
pars_yr <- pars_yr[-which_contains_numeric0]
years_v <- years_v[-which_contains_numeric0]


calc_lambda <- function(i) {
  lam <- Re(eigen(kernel(pars_yr[[i]])$k_yx)$value[1])
  return(lam)
}

lambdas_yr <- lapply(1:(length(pars_yr)), calc_lambda)
names(lambdas_yr) <- years_v


# Comparing the year-specific lambdas
year_kern <- function(i) {
  return(kernel(pars_yr[[i]])$k_yx)
}

kern_yr <- lapply(1:(length(years_v)), year_kern)

all_mat <- array(dim = c(200, 200, (length(years_v))))

for(i in 1:(length(years_v))) {
  all_mat[,,i] <- as.matrix(kern_yr[[i]])
}

mean_kern <- apply(all_mat, c(1, 2), mean)
lam_mean_kern <- Re(eigen(mean_kern)$value[1])


# Population counts at time t0
pop_counts_t0 <- df %>%
  group_by(year, quad) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  group_by(year, quad ) %>%
  summarize(n_t1 = n()) %>% 
  ungroup 

# Calculate observed population growth rates, 
#   accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
  # by dropping NAs, we remove gaps in sampling!
  mutate(year = year - 1) %>% 
  drop_na %>% 
  group_by(year) %>% 
  summarise(n_t0 = sum(n_t0),
            n_t1 = sum(n_t1)) %>% 
  ungroup %>% 
  mutate(obs_pgr = n_t1 / n_t0) %>%
  full_join(data.frame(year = as.numeric(names(lambdas_yr)),
                       lambda = unlist(lambdas_yr)), 
            by = 'year') %>% 
  drop_na

lam_mean_yr    <- mean(pop_counts$lambda, na.rm = T)
lam_mean_count <- mean(pop_counts$obs_pgr, na.rm = T)

lam_mean_geom <- exp(mean(log(pop_counts$obs_pgr), na.rm = T))
lam_mean_geom

lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)
lam_mean_overall

# changing the year according to 
years_v <- pop_counts$year

# projecting a population vector for each year using the year-specific models,
# and compare the projected population to the observed population
count_indivs_by_size <- function(size_vector,
                                 lower_size,
                                 upper_size,
                                 matrix_size){
  
  size_vector %>%
    cut(breaks = seq(lower_size - 0.00001,
                     upper_size + 0.00001,
                     length.out = matrix_size + 1)) %>%
    table %>%
    as.vector
  
}

yr_pop_vec <- function(i) {
  vec_temp <- surv_df %>% filter(year == i) %>% select(logsize_t0) %>% unlist()
  min_sz   <- pars_mean$L
  max_sz   <- pars_mean$U
  pop_vec <- count_indivs_by_size(vec_temp, min_sz, max_sz, 200)
  
  return(pop_vec)
}

year_pop <- lapply(years_v, yr_pop_vec)

proj_pop <- function(i) {
  sum(all_mat[,,i] %*% year_pop[[i]])
}

projected_pop_ns  <- sapply(1:(length(years_v)), proj_pop)

pop_counts_update <- pop_counts %>%
  mutate(proj_n_t1 = projected_pop_ns) %>%
  mutate(proj_pgr  = proj_n_t1/n_t0)

g_mod_vs_obs <- ggplot(pop_counts_update) +
  geom_point( aes(x = lambda,   y = obs_pgr), color = 'brown') +
  geom_point( aes(x = proj_pgr, y = obs_pgr), color = 'red') +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title    = 'Lambda vs pgr',
       subtitle = v_ggp_suffix,
       x        = 'Modeled lambda',
       y        = 'Observed population growth rate') +
  theme_classic()

ggsave(paste0(dir_result, '/5_years_lamda_mod_vs_obs', v_suffix, '.png'), 
       plot = g_mod_vs_obs,
       width = 4, height = 3, dpi = 150)


# Building the year-specific IPMs with ipmr ------------------------------------
# all of our varying and constant parameters into a single list
all_pars <- c(pars_cons_wide, pars_var_wide)
# add the index of the best model to the parameter
all_pars[['mod_su_index']] <- v_mod_su_index
all_pars[['mod_gr_index']] <- v_mod_gr_index

write.csv(all_pars, 
          paste0(dir_data, '/all_pars', v_suffix, '.csv'), 
          row.names = F)

# proto-IPM with '_yr' suffix
proto_ipm_yr <- init_ipm(sim_gen   = 'simple',
                         di_dd     = 'di',
                         det_stoch = 'det') %>% 
  
  define_kernel(
    name             = 'P_yr',
    family           = 'CC',
    formula          = s_yr * g_yr,
    s_yr             = plogis(
      surv_b0_yr + 
        (if (mod_su_index >= 1) surv_b1_yr * size_1   else 0) +
        (if (mod_su_index >= 2) surv_b2_yr * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3_yr * size_1^3 else 0)),
    
    g_yr             = dnorm(size_2, mu_g_yr, grow_sig),
    mu_g_yr          = grow_b0_yr + 
      (if (mod_gr_index >= 1) grow_b1_yr * size_1   else 0) +
      (if (mod_gr_index >= 2) grow_b2_yr * size_1^2 else 0) +
      (if (mod_gr_index >= 3) grow_b3_yr * size_1^3 else 0),
    
    grow_sig         = sqrt(a * exp(b * size_1)),
    data_list        = all_pars,
    states           = list(c('size')),
    
    # these next two lines are new
    # the first tells ipmr that we are using parameter sets
    uses_par_sets    = TRUE,
    # the second defines the values the yr suffix can assume
    par_set_indices  = list(yr = years_v),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions(fun    = 'norm',
                                               target = 'g_yr')
  ) %>% 
  
  define_kernel(
    name             = 'F_yr',
    family           = 'CC',
    formula          = fecu_b0_yr * r_d,
    r_d              = dnorm(size_2, recr_sz, recr_sd),
    data_list        = all_pars,
    states           = list(c('size')),
    uses_par_sets    = TRUE,
    par_set_indices  = list(yr = years_v),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm', 'r_d')
  ) %>% 
  
  define_impl(
    make_impl_args_list(
      kernel_names = c(  'P_yr', 'F_yr'),
      int_rule     = rep('midpoint', 2),
      state_start  = rep('size', 2),
      state_end    = rep('size', 2)
    )
  ) %>% 
  
  define_domains(
    size = c(all_pars$L,
             all_pars$U,
             all_pars$mat_siz
    )
  ) %>% 
  
  # We also append the suffix in define_pop_state(). This will create a deterministic
  # simulation for every 'year'
  define_pop_state(
    n_size_yr = rep(1 / 200, 200)
  )


# Make a dataframe
ipmr_yr       <- make_ipm(proto_ipm  = proto_ipm_yr,
                          iterations = 200)
lam_mean_ipmr <- lambda(ipmr_yr)
lam_out       <- data.frame(coefficient = names(lam_mean_ipmr), 
                            value       = lam_mean_ipmr,
                            years       = years_v)
rownames( lam_out) <- 1:(length(years_v))
write.csv(lam_out, 
          paste0(dir_data, '/lambdas_yr_vec', v_suffix, '.csv'), 
          row.names = F)

lam_out_wide  <- as.list(pivot_wider(lam_out, names_from = 'coefficient', 
                                     values_from         = 'value'))
write.csv(lam_out_wide, 
          paste0(dir_data, '/lambdas_yr', v_suffix, '.csv'), 
          row.names = F)
