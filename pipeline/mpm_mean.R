# MPM mean pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.16

# Original models from https://doi.org/10.1111/j.1365-2745.2009.01585.x
#  check the appendices for details on model fitting!

# A mean MPM with mean survival per age class


# Packages ---------------------------------------------------------------------
# read and check  
source('helper_functions/load_packages.R')
load_packages(tidyverse, janitor,
              ggthemes, # scale_color_colorblind
              skimr)
# load_packages(patchwork, , lme4, ggthemes, boot)

# Specification ----------------------------------------------------------------
# Species abbreviation
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))
# Suffix for the folder structure
v_folder_suffix <- paste0(v_sp_abb, '_', v_mod_type)
# Prefix for the script name
v_script_prefix <- str_c(
  str_extract(v_author_year, '^[^_]+'),
  str_sub(str_extract(v_author_year, '_\\d+$'), -2, -1))
# Define prefix for two of the same author and year
if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0('^', v_author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}
# Figure subtitle
v_fig_subt <- paste(
  paste(
    paste0(
      toupper(substr(strsplit(v_author_year, '_')[[1]][1], 1, 1)), 
      substr(strsplit(
        v_author_year, '_')[[1]][1], 
        2, 
        nchar(strsplit(v_author_year, '_')[[1]][1]))), 
    strsplit(v_author_year, '_')[[1]][2]), 
  '-', 
  v_species)


# Setting the age class, max age starting from 0 to inf
if (!exists('v_age_class') || length(v_age_class) == 0) {
  v_age_class <- 1
}


# Directories ------------------------------------------------------------------ 
dir_pub       <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R         <- file.path(dir_pub, 'R',       v_folder_suffix)
dir_data      <- file.path(dir_pub, 'data',    v_folder_suffix)
dir_result    <- file.path(dir_pub, 'results', v_folder_suffix)
dir_fun       <- file.path('helper_functions')

# Create output directory if it doesn't exist
if (!dir.exists(paste0(dir_result))) {
  dir.create(paste0(dir_result))}
# Create the data folder for the species if it does not exits already
if (!dir.exists(paste0(dir_data))) {
  dir.create(paste0(dir_data))}

# Plant tracker if it does not already exist
if (!file.exists(
  file.path(dir_data, paste0(v_script_prefix, '_', v_sp_abb, '.csv')))) {
  source(
    file.path(dir_R, paste0(v_script_prefix, '_', v_sp_abb, '_tracker.R')))}


# Functions --------------------------------------------------------------------
# Estimate survival probability with CI by age class
source(file.path(dir_fun, 'fun_sur_prob_ci.R'))


# Data -------------------------------------------------------------------------
df <- read.csv(
  file.path(dir_data, paste0(v_script_prefix, '_', v_sp_abb, '.csv'))) %>% 
  clean_names() %>% 
  filter(species == v_species) %>%
  select(-c(suspect, near_edge, site)) %>% # geometry, 
  mutate(across(c(quad), as.factor)) %>% 
  rename(survives = survives_tplus1)


# Data exploration -------------------------------------------------------------
# Quadrat inventory
#  grouping by quard, year and counting the individuals
fig_inv_per_year <- df %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  ggplot() + 
  geom_point(aes(x = year, y = quad)) +
  theme_bw() +
  labs(x        = 'Year',
       y        = 'Quadrat',
       title    = 'Sample inventory',
       subtitle = v_fig_subt) +
  theme(axis.text.y = element_text(size = 5))

ggsave(file.path(
  dir_result, '0_inventory_quadrat_per_year.png'),
  plot = fig_inv_per_year, width = 6, height = 4, dpi = 150)


# Getting the DF's ready -------------------------------------------------------
# Survival data that retains all ages
df_suv <- df %>% 
  mutate(survives = replace(
    survives, is.na(survives) & year == (df$year %>% max) - 1, 0)) %>% 
  subset(!is.na(survives)) %>%
  select(quad, track_id, year, age, survives) 


# Recruits data frame
df_rec  <- df %>%
  # Total number of recruits at time t1
  count(quad, year) %>% 
  mutate(year = as.integer(year)) %>% 
  mutate(year = year + 1) %>% 
  rename(parents = n) %>% 
  left_join(
    df %>%
      subset(recruit == 1) %>% 
      count(quad, year) %>% 
      mutate(year = as.integer(year)) %>% 
      rename(recruits = n) %>% 
      # include all available quad X year combinations
      right_join(
        # Total number of individuals time t0
        df %>% 
        dplyr::select(quad, year) %>% 
        mutate(year = as.integer(year)) %>% 
        unique) %>% 
        # if is.na(recruits), then recruits == 0
        mutate(recruits = replace(
          recruits, is.na(recruits), 0))) %>% 
  # change this, so that it is year_t0 instead of year_t1
  mutate(year = year - 1) %>% 
  mutate(year = as.character(year)) %>% 
  # calculate per capita recruitment. 
  #   Support: 0 to infinity.
  mutate(pcr = recruits / parents) 


# Removing specified years -----------------------------------------------------
v_years_all  <- sort(unique(df$year))

df      <- df     %>% filter(!is.na(year) & !(year %in% v_years_re))
df_suv  <- df_suv %>% filter(!is.na(year) & !(year %in% v_years_re))
df_rec  <- df_rec %>% filter(!is.na(year) & !(year %in% v_years_re))

v_years  <- sort(unique(df$year))


# Save all data ----------------------------------------------------------------
write.csv(df,     row.names = F, paste0(dir_data, '/mean_data_df.csv'))
write.csv(df_suv, row.names = F, paste0(dir_data, '/mean_survival_df.csv'))
write.csv(df_rec, row.names = F, paste0(dir_data, '/mean_recruitment_df.csv'))


# Age --------------------------------------------------------------------------
# Age - Counts
fig_age_hist <- ggplot(df, aes(x = age)) +
  geom_histogram(binwidth = .5, na.rm = TRUE) +
  labs(x = 'Age',
       y = 'Count', 
       title    = 'Age hist',
       subtitle = v_fig_subt) +
  theme_bw()

ggsave(file.path(dir_result, '0_age_counts.png'),
       plot = fig_age_hist, width = 6, height = 9, dpi = 150)


# Age - Survival
fig_age_surv <- ggplot(
  data = na.omit(bind_rows(lapply(
    unique(df_suv$age[!is.na(df_suv$age)]), fun_sur_prob_ci))), 
  aes(x = age, y = mean_survival)) +
  geom_point(color = 'blue', size = 2) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.2, color = 'blue') + 
  geom_text(aes(label = paste('n =', n_obs)), hjust = 1.2, 
            color = 'black', size = 3) +  
  geom_smooth(data = df_suv, aes(x = age, y = survives), method = 'glm', 
              method.args = list(family = 'binomial'), color = 'red', se = T) +  
  labs(x = 'Age', y = 'Survival Probability', 
       title = 'Survival Probability by Age with 95% CI',
       subtitle = v_fig_subt) +
  theme() +
  ylim(0, 1) +
  theme_bw()

ggsave(file.path(dir_result, '0_age_surv.png'),  
       plot = fig_age_surv, width = 6, height = 9, dpi = 150)


# Survival ---------------------------------------------------------------------
# Data
df_su_mod <- df_suv %>%
  mutate(age = replace(age, age > v_age_class, v_age_class))

# Model: GLM
mod_su   <- glm(survives ~ age, data = df_su_mod, family = 'binomial')

# Mean predicted probability by age class–specific
# Get predicted survival probabilities with standard errors
mod_su_new_ages <- data.frame(age = sort(unique(df_su_mod$age)))
mod_su_pred <- predict(mod_su, newdata = mod_su_new_ages, type = "link", se.fit = TRUE)

# Create data frame with predictions and SE's
mod_su_prob <- data.frame(
  age     = mod_su_new_ages$age,
  fit     = mod_su_pred$fit,
  se      = mod_su_pred$se.fit) %>% 
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se,
    pred_su = plogis(fit),
    ci_lower = plogis(lower),
    ci_upper = plogis(upper))

# Figure: Mean survival by age class
fig_su <- ggplot(data = df_su_mod %>% filter(!is.na(age)), 
       aes(x = as.factor(age), y = survives)) +
  geom_jitter(height = 0.1, width = 0.2, alpha = 0.3) +
  geom_errorbar(data = mod_su_prob %>% mutate(age = age),
                aes(x = as.factor(age), ymin = ci_lower, ymax = ci_upper),
                inherit.aes = FALSE,
                width = 0.2, color = 'red') +
  geom_point(data = mod_su_prob %>% mutate(age = age),
             aes(x = as.factor(age), y = pred_su),
             inherit.aes = FALSE,
             color = 'black', size = 1.5) +
  geom_text(data = mod_su_prob %>% mutate(age = age),
            aes(x = as.factor(age), y = pred_su, 
                label = round(pred_su, 2)),
            inherit.aes = FALSE,
            color = 'black',
            size = 3,
            hjust = 0,  
            nudge_x = 0.15) +
  geom_text(data = df_su_mod %>%
              filter(!is.na(age)) %>%
              group_by(age) %>%
              summarise(n = n()) %>%
              mutate(age = as.factor(age)),
            aes(x = age, y = 0, label = paste0("n=", n)),
            inherit.aes = FALSE,
            hjust = 2.5,
            vjust = 11,  # only shows in the image
            size = 3) +
  theme_bw() +
  labs(title = 'Survival by age class - mean predictions and 95CI',
       subtitle = v_fig_subt,
       y = 'Survival to t1',
       x = 'Age class')

ggsave(file.path(dir_result, 'mean_survival.png'),  
       plot = fig_su, width = 6, height = 9, dpi = 150)


# Recruitment two-stage model --------------------------------------------------
#   1. Model recruitment being there or not
#   2. Model per-capita recruitment (pcr) *conditional* on recruitment happening

# 1. yes or no recruits?
df_re_yn <- df_rec %>% 
  mutate(recruits_yesno = as.numeric(recruits > 0)) %>% 
  dplyr::select(quad, year, parents, recruits_yesno) %>% 
  drop_na

# Yes-or-no recruits? (NOTE: parent number matters)
mod_re_yn <- glm(recruits_yesno ~ 1, data = df_re_yn, family = binomial)

# Predict, which is only the intercept
mod_re_yn_pred <- plogis(coef(mod_re_yn)[[1]])


# 2.1 per capita recruitment **conditional** on recruitment happening
df_re_pc_c <- df_rec %>% subset(!(recruits == 0))

# Per-capita recruitment model
mod_re_pc_c <- glm(pcr ~ 1, data = df_re_pc_c, family = Gamma(link = 'log'))

# predict conditional per-capita recruitment
mod_re_pc_c_pred <- exp(coef(mod_re_pc_c)[[1]])


# 2.2 per capita recruitment unconditional 
df_re_pc_uc <- df_rec %>% 
  mutate(recruits_uc = replace(
    recruits, recruits == 0, 0.1)) %>% 
  mutate(pcr_uc = recruits_uc / parents)

# Per-capita recruitment model
mod_re_pc_uc <- glm(pcr_uc ~ 1, 
                    data = df_re_pc_uc, family = Gamma(link = 'log'))

# predict conditional per-capita recruitment
mod_re_pc_uc_pred <- exp(coef(mod_re_pc_c)[[1]])


# hurdle model predictions
df_re_pred <- data.frame(
  model = c('yn', 'pc_c', 'pc_uc', 'pc_pred'),
  pred  = c(mod_re_yn_pred, mod_re_pc_c_pred, mod_re_pc_uc_pred,
            # pred conditional * pred yes no
            mod_re_pc_c_pred * mod_re_yn_pred))


# Plots
fig_re <- ggplot(data = df_rec, aes(x = parents, y = recruits)) + 
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
  geom_jitter(width = .5, height = .5) +
  geom_line(aes(y = df_rec$parents + df_re_pred$pred[4]), 
            color = 'blue') + 
  theme_bw() + 
  labs(title    = 'Recruitment - prediction',
       subtitle = v_fig_subt,
       x        = 'Nr. of individuals in t0',
       y        = 'Recruits in t1')

ggsave(file.path(dir_result, 'mean_recruitment.png'),  
       plot = fig_re, width = 6, height = 9, dpi = 150)


# Matrix population model ------------------------------------------------------
# MPM parameters
df_mpm <- mod_su_prob %>% 
  select(age, pred_su) %>% 
  pivot_wider(names_from  = 'age',
              values_from = 'pred_su') %>% 
  mutate(re_pc_pred = df_re_pred[4,2])

# project based on stage distribution
stage_counts <- df %>%
  mutate(age = replace(age, age > v_age_class, v_age_class)) %>% 
  group_by(age) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  drop_na %>% 
  pivot_wider(names_from  = 'age',
              values_from = 'n_t0',
              values_fill = 0) %>% 
  drop_na 


# Functions
# lambda function
fun_lambda_mpm_mean <- function(x) x %>% eigen %>% .$value %>% Re %>% max


# create the lambda matrix 
fun_lam_mean <- function(mat_df){
  
  # mat_df <- df_mpm
  tmp_pars <- mat_df
  
  mat      <- matrix(0,v_age_class + 1,v_age_class + 1)
  
  # per capita recruitment values
  mat[1,]  <- tmp_pars$re_pc_pred
  
  # survival values
  mat[2,1] <- tmp_pars$`0`
  
  for (i in 1:v_age_class) {
    mat[i + 1, i] <- tmp_pars[[as.character(i - 1)]]
  }
  
  # Fill survival for the last age class staying in place
  mat[v_age_class + 1, v_age_class + 1] <- tmp_pars[[as.character(v_age_class )]]
  
  mat %>% fun_lambda_mpm_mean
}


# Make matrix. 
#   based on parameters in rows of a data frame
fun_make_mat_mean <- function(mat_df){
  
  # mat_df <- df_mpm
  tmp_pars <- mat_df
  
  mat      <- matrix(0,v_age_class + 1,v_age_class + 1)
  
  # per capita recruitment values
  mat[1,]  <- tmp_pars$re_pc_pred
  
  # survival values
  mat[2,1] <- tmp_pars$`0`
  
  for (i in 1:v_age_class) {
    mat[i + 1, i] <- tmp_pars[[as.character(i - 1)]]
  }
  
  # Fill survival for the last age class staying in place
  mat[v_age_class + 1, v_age_class + 1] <- tmp_pars[[as.character(v_age_class )]]
  
  # output the matrix
  mat
}


# projected lambda (takes into account observed stage distribution)
fun_proj_lam_mean <- function(mat_df){
  
  # year-specific stage counts
  tmp_stage <- stage_counts
  
  # year-specific projection matrix
  # mat_df <- df_mpm
  tmp_mat <- fun_make_mat_mean(mat_df)
  
  nt0 <- sum(tmp_stage)
  nt1 <- (tmp_mat %*% as.numeric(tmp_stage)) %>% sum
  
  nt1 / nt0
}

df_lam <- df_mpm %>% 
  mutate(lam           = fun_lam_mean(df_mpm),
         proj_lam_mean = fun_proj_lam_mean(df_mpm))


# Population counts at time t0
pop_counts_t0 <- df %>%
  # subset( !is.na(age) ) %>% 
  group_by(year, quad) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  # subset( !is.na(age) ) %>% 
  group_by(year, quad) %>%
  summarize(n_t1 = n()) %>% 
  ungroup 

# Calculate observed population growth rates, 
#   accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
  mutate(year = year - 1) %>% 
  # by dropping NAs, we remove gaps in sampling!
  drop_na %>% 
  summarise(n_t0 = sum(n_t0),
            n_t1 = sum(n_t1)) %>%
  mutate(obs_pgr = n_t1 / n_t0) %>% 
  bind_cols(df_lam) %>% 
  drop_na

write.csv(pop_counts, row.names = F, file.path(dir_data, 'mean_metrics_df.csv'))


# # Save the data
# df_mpm %>% 
#   mutate(SpeciesAuthor = paste0(v_species, ' - ', v_author_year)) %>% 
#   setNames(c('surv_age0', 'surv_age1', 
#              'per_capita_recruitment',
#              'species_author')) %>% 
#   write.csv(file.path(dir_data, paste0('mpm_', v_sp_abb, '.csv')), 
#             row.names = F)

