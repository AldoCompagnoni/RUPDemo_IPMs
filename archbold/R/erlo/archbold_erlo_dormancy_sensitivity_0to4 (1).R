# Dormancy cutoff sensitivity - Archbold - Eriogonum longifolium
# Sensitivity analysis for 0-4 years of theoretical dormancy
#
# This script changes ONLY the terminal censoring rule used for observed
# population growth. The fitted vital-rate models and IPM are not changed.
#
# Interpretation of dormancy_cutoff:
#   0 = disappearance can be treated as death immediately
#   1 = allow 1 year of unresolved dormancy
#   2 = allow 2 years of unresolved dormancy
#   3 = current/default assumption
#   4 = allow 4 years of unresolved dormancy
#
# Increasing the cutoff removes additional transitions from the end of the
# time series, where absence cannot yet be confidently distinguished from
# extended dormancy.

# Packages --------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork)


# Specification ----------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Eriogonum longifolium')

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

dir_data <- file.path(v_head, 'data', v_sp_abb)


# Data -------------------------------------------------------------------------
df <- read.csv(
  file.path(
    dir_data,
    paste0('ab_', v_sp_abb, '_df_workdata_260820.csv'))) %>%
  mutate(
    year = as.numeric(year),
    row_type = as.character(row_type))

df_ind <- df %>%
  filter(row_type == 'individual')


# Annual quadrat counts ---------------------------------------------------------
# Count known living individuals. New adults are additionally counted one year
# earlier, matching the simplified recruitment assumption used in the mean IPM.
df_counts_alive <- df_ind %>%
  filter(state %in% c('active', 'dormant', 'alive_or_dormant')) %>%
  count(site, pop, qu, year, name = 'n_alive')

df_newadult_backfill <- df_ind %>%
  filter(recruit_type == 'new_adult') %>%
  transmute(site, pop, qu, year = year - 1) %>%
  count(site, pop, qu, year, name = 'n_newadult_backfill')

# Recruitment rows define the monitored quadrat-year inventory.
df_counts_quad <- df %>%
  filter(row_type == 'recruitment') %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_counts_alive,
    by = c('site', 'pop', 'qu', 'year')) %>%
  left_join(
    df_newadult_backfill,
    by = c('site', 'pop', 'qu', 'year')) %>%
  mutate(
    n_alive = replace_na(n_alive, 0L),
    n_newadult_backfill = replace_na(n_newadult_backfill, 0L),
    n = n_alive + n_newadult_backfill) %>%
  group_by(site, pop, qu) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    year_t1 = lead(year),
    n_t1 = lead(n),
    annual_transition = year_t1 == year + 1) %>%
  ungroup() %>%
  filter(annual_transition)


# Dormancy-cutoff sensitivity --------------------------------------------------
dormancy_cutoffs <- 0:4
max_ind_year <- max(df_ind$year, na.rm = TRUE)

get_obs_lambda <- function(dormancy_cutoff) {
  # Same rule as the original script:
  # last_complete_transition = max observed individual year - cutoff - 1
  last_complete_transition <- max_ind_year - dormancy_cutoff - 1

  df_year <- df_counts_quad %>%
    filter(year <= last_complete_transition) %>%
    group_by(year) %>%
    summarise(
      n_quads = n(),
      n_t0 = sum(n),
      n_t1 = sum(n_t1), .groups = 'drop') %>%
    mutate(lambda_obs = n_t1 / n_t0)

  tibble(
    dormancy_cutoff = dormancy_cutoff,
    last_complete_transition = last_complete_transition,
    first_transition = min(df_year$year, na.rm = TRUE),
    n_transitions = nrow(df_year),
    lambda_obs_arithmetic = mean(df_year$lambda_obs, na.rm = TRUE),
    lambda_obs_geometric = exp(mean(log(df_year$lambda_obs), na.rm = TRUE)),
    lambda_obs_min = min(df_year$lambda_obs, na.rm = TRUE),
    lambda_obs_max = max(df_year$lambda_obs, na.rm = TRUE))
}

df_dorm_sensitivity <- map_dfr(dormancy_cutoffs, get_obs_lambda)

df_dorm_sensitivity


