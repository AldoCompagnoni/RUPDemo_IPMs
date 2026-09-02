# Data - Archbold - Solidago odora

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.08.31

# Study organism: Solidago odora var. chapmanii
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.246.1
# Meta data link:
# https://portal.edirepository.org/nis/metadataviewer?packageid=edi.246.1
# Citing publication: Menges & Root 2004, The American Midland Naturalist
# Time period: 1991-2000


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
set.seed(100)
options(stringsAsFactors = F)


# Packages --------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(patchwork, skimr, janitor, GGally, tidyverse)


# Specification ---------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Solidago odora')
custom_delimiter <- c()

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)


# Directory -------------------------------------------------------------------
dir_pub <- file.path(paste0(v_head))
dir_R <- file.path(dir_pub, 'R', v_sp_abb)
dir_data <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

if (!dir.exists(dir_R)) {dir.create(dir_R)}
if (!dir.exists(dir_data)) {dir.create(dir_data)}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Data ------------------------------------------------------------------------
data_file <- file.path(dir_data, 'solidago_odora_data.csv')

if (!file.exists(data_file)) {
  data_file_alt <- file.path(dir_data, 'solidago_odora_data(1).csv')
  if (file.exists(data_file_alt)) {data_file <- data_file_alt}
}

df_og <- read_csv(data_file, show_col_types = FALSE) %>%
  janitor::clean_names()

fire_vars <- c('fire90', 'fire91', 'fire92', 'fire95', 'fire98')


# Metadata --------------------------------------------------------------------
df_meta <- tibble::tribble(
  ~var, ~description,
  'quad', 'Quadrat number; permanent 1 m2 quadrat',
  'fire90', 'Did fire occur in this plot in 1990?',
  'fire91', 'Did fire occur in this plot in 1991?',
  'fire92', 'Did fire occur in this plot in 1992?',
  'fire95', 'Did fire occur in this plot in 1995?',
  'fire98', 'Did fire occur in this plot in 1998?',
  'nfires', 'Number of fires recorded for the quadrat',
  'year', 'Year of observation',
  'st', 'Number of stems in plot',
  'fl', 'Number of flowering stems in plot',
  'mht', 'Mean height of stems in plot, cm')


# Initial checks ---------------------------------------------------------------
dim(df_og)
names(df_og)
skimr::skim(df_og)

df_duplicates <- df_og %>%
  count(quad, year, name = 'n_entries') %>%
  filter(n_entries > 1)

df_expected_grid <- tidyr::expand_grid(
  quad = sort(unique(df_og$quad)),
  year = seq(min(df_og$year), max(df_og$year)))

df_missing_rows <- df_expected_grid %>%
  anti_join(df_og %>% distinct(quad, year), by = c('quad', 'year'))

df_value_checks <- tibble(
  check = c(
    'st < 0', 'fl < 0', 'fl > st', 'st non-integer',
    'fl non-integer', 'mht <= 0', 'st == 0 but fl > 0',
    'st == 0 but mht observed', 'st > 0 but mht missing'),
  n = c(
    sum(df_og$st < 0, na.rm = TRUE),
    sum(df_og$fl < 0, na.rm = TRUE),
    sum(df_og$fl > df_og$st, na.rm = TRUE),
    sum(abs(df_og$st - round(df_og$st)) > 1e-8, na.rm = TRUE),
    sum(abs(df_og$fl - round(df_og$fl)) > 1e-8, na.rm = TRUE),
    sum(df_og$mht <= 0, na.rm = TRUE),
    sum(df_og$st == 0 & df_og$fl > 0, na.rm = TRUE),
    sum(df_og$st == 0 & !is.na(df_og$mht), na.rm = TRUE),
    sum(df_og$st > 0 & is.na(df_og$mht), na.rm = TRUE)))

df_duplicates
df_missing_rows
df_value_checks


# Fire-column checks -----------------------------------------------------------
# All fire columns are binary and constant within quadrat. Despite the EDI
# description of nfires as "1991-1999", the data show that nfires equals the
# sum of fire90 + fire91 + fire92 + fire95 + fire98.

df_fire_quad <- df_og %>%
  distinct(quad, across(all_of(c(fire_vars, 'nfires')))) %>%
  mutate(
    nfires_from_cols = rowSums(across(all_of(fire_vars))),
    nfires_difference = nfires_from_cols - nfires)

df_fire_quad

df_fire_consistency <- df_og %>%
  group_by(quad) %>%
  summarise(
    across(
      all_of(c(fire_vars, 'nfires')),
      ~ n_distinct(.x, na.rm = TRUE),
      .names = 'n_values_{.col}'),
    .groups = 'drop') %>%
  filter(if_any(starts_with('n_values_'), ~ .x > 1))

df_fire_consistency


# Published quadrat information -----------------------------------------------
# Menges & Root (2004), Table 1.
df_paper_quad <- tibble::tribble(
  ~quad, ~burn_unit, ~vegetation, ~soil, ~paper_years,
  ~paper_mean_stems, ~paper_mean_height, ~formerly_grazed,
  1, '41A', 'O-H Scrub', 'Tavares', '1991-2000', 4.0, 30.1, FALSE,
  2, '49B', 'O-H Scrub', 'Astatula', '1991-2000', 19.4, 27.6, FALSE,
  3, '53A', 'O-H Scrub', 'Astatula', '1991-2000', 51.1, 22.4, FALSE,
  4, '35B', 'Scrubby F.', 'Duette', '1991-1996', 12.8, 30.2, FALSE,
  5, '35B', 'Scrubby F.', 'Duette', '1991-2000', 12.0, 15.9, FALSE,
  6, '36', 'O-H Scrub', 'Orsino', '1991-2000', 6.5, 34.5, FALSE,
  7, '36', 'O-H Scrub', 'Orsino', '1991-2000', 9.2, 30.8, FALSE,
  8, '26D', 'Scrubby F.', 'Duette', '1991-2000', 4.7, 27.8, FALSE,
  9, '26A', 'Scrubby F.', 'Orsino', '1991-2000', 15.5, 19.6, FALSE,
  10, '22A', 'Scrubby F.', 'Duette', '1991-2000', 17.1, 18.8, FALSE,
  11, '22A', 'Scrubby F.', 'Duette', '1991-2000', 10.9, 21.4, FALSE,
  12, '2D', 'SP Scrub', 'Paola', '1991-2000', 13.5, 48.6, FALSE,
  13, '2C', 'Sandhill', 'Astatula', '1991-2000', 18.0, 28.2, FALSE,
  14, '57A', 'Disturbed', 'Paola', '1992-2000', 162.9, 27.9, TRUE,
  15, '57A', 'Disturbed', 'Paola', '1992-2000', 62.1, 27.4, TRUE,
  16, '24C', 'Flatwoods', 'Duette', '1992-1998', 10.2, 46.8, FALSE,
  17, '24B', 'O-H Scrub', 'Tavares', '1992-2000', 3.3, 30.4, FALSE,
  18, '25A', 'O-H Scrub', 'Tavares', '1992-2000', 29.1, 26.5, FALSE,
  19, '25A', 'O-H Scrub', 'Tavares', '1992-2000', 4.6, 43.6, FALSE,
  20, '49B', 'O-H Scrub', 'Astatula', '1995-2000', 5.7, 10.1, FALSE,
  21, '49B', 'O-H Scrub', 'Astatula', '1995-2000', 10.0, 9.4, FALSE,
  22, '2A', 'Sandhill', 'Astatula', '1995-2000', 13.0, 27.1, FALSE,
  23, '2A', 'Sandhill', 'Astatula', '1995-2000', 25.2, 26.2, FALSE)


# Published fire history -------------------------------------------------------
# Month/year values are from Menges & Root (2004), Table 1.
df_paper_fire <- tibble::tribble(
  ~quad, ~fire_year, ~fire_month,
  1, 1989, 5,
  2, 1990, 3,
  2, 1995, 6,
  4, 1990, 11,
  5, 1990, 11,
  6, 1991, 4,
  7, 1991, 4,
  8, 1991, 7,
  9, 1991, 7,
  10, 1990, 10,
  11, 1990, 10,
  12, 1991, 5,
  12, 1998, 5,
  13, 1991, 5,
  14, 1992, 2,
  15, 1992, 2,
  16, 1992, 5,
  17, 1992, 5,
  18, 1991, 7,
  19, 1991, 7,
  20, 1995, 6,
  21, 1995, 6,
  22, 1991, 5,
  22, 1995, 5,
  22, 1998, 7,
  23, 1991, 5,
  23, 1995, 5,
  23, 1998, 7) %>%
  mutate(
    fire_date = as.Date(sprintf(
      '%d-%02d-15', fire_year, fire_month)),
    fire_season = case_when(
      fire_month %in% 2:4 ~ 'Spring',
      fire_month == 5 ~ 'Late spring',
      fire_month %in% 6:7 ~ 'Summer',
      fire_month %in% c(10, 11) ~ 'Fall',
      TRUE ~ 'Other'))


# Published values vs data -----------------------------------------------------
df_paper_mean_check <- df_og %>%
  group_by(quad) %>%
  summarise(
    data_mean_stems = mean(st, na.rm = TRUE),
    data_mean_height = mean(mht, na.rm = TRUE),
    .groups = 'drop') %>%
  left_join(df_paper_quad, by = 'quad') %>%
  mutate(
    stem_mean_difference = data_mean_stems - paper_mean_stems,
    height_mean_difference = data_mean_height - paper_mean_height)

df_paper_mean_check


# CSV fire columns vs published fire history ----------------------------------
df_paper_fire_flags <- df_paper_fire %>%
  filter(fire_year %in% c(1990, 1991, 1992, 1995, 1998)) %>%
  mutate(
    fire_variable = paste0('fire', str_sub(fire_year, 3, 4)),
    paper_fire = 1) %>%
  select(quad, fire_variable, paper_fire) %>%
  complete(
    quad = sort(unique(df_og$quad)),
    fire_variable = fire_vars,
    fill = list(paper_fire = 0))

df_fire_csv_vs_paper <- df_og %>%
  distinct(quad, across(all_of(fire_vars))) %>%
  pivot_longer(
    cols = all_of(fire_vars),
    names_to = 'fire_variable',
    values_to = 'csv_fire') %>%
  left_join(
    df_paper_fire_flags,
    by = c('quad', 'fire_variable')) %>%
  mutate(matches_paper = csv_fire == paper_fire)

df_fire_csv_vs_paper


# Census dates ----------------------------------------------------------------
# Censuses were in August except September 1991 and October 1992.
df_year <- df_og %>%
  mutate(
    census_month = case_when(
      year == 1991 ~ 9L,
      year == 1992 ~ 10L,
      TRUE ~ 8L),
    census_date = as.Date(sprintf(
      '%d-%02d-15', as.integer(year), census_month))) %>%
  left_join(
    df_paper_quad %>%
      select(
        quad, burn_unit, vegetation, soil, formerly_grazed),
    by = 'quad')


# Most recent fire before each census -----------------------------------------
# Quadrat 3 has no documented fire in Table 1 and therefore retains NA for
# time-since-fire variables.
df_fire_join <- bind_rows(
  df_paper_fire,
  tibble(
    quad = sort(unique(df_og$quad)),
    fire_year = NA_real_,
    fire_month = NA_real_,
    fire_date = as.Date(NA),
    fire_season = NA_character_))

df_year <- df_year %>%
  left_join(df_fire_join, by = 'quad', relationship = 'many-to-many') %>%
  filter(is.na(fire_date) | fire_date <= census_date) %>%
  mutate(
    fire_order = if_else(
      is.na(fire_date), -Inf, as.numeric(fire_date))) %>%
  group_by(quad, year) %>%
  slice_max(fire_order, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    months_since_fire = if_else(
      !is.na(fire_year),
      12 * (year - fire_year) + census_month - fire_month,
      NA_real_),
    years_since_fire = months_since_fire / 12,
    postfire_year = floor(years_since_fire),
    recent_fire_lt4mo = !is.na(months_since_fire) &
      months_since_fire < 4,
    analysis_usable = !is.na(st) & !recent_fire_lt4mo,
    postfire_usable = analysis_usable & !is.na(years_since_fire),
    height_usable = analysis_usable & st > 0 & !is.na(mht),
    flower_usable = analysis_usable & st > 0 & !is.na(fl),
    flower_prop = case_when(
      flower_usable ~ fl / st,
      TRUE ~ NA_real_),
    log_stems = if_else(
      analysis_usable, log1p(st), NA_real_),
    postfire_period = case_when(
      is.na(postfire_year) ~ NA_character_,
      postfire_year == 0 ~ '0',
      postfire_year == 1 ~ '1',
      postfire_year >= 2 ~ '2+')) %>%
  select(-fire_order) %>%
  mutate(
    quad = factor(quad),
    formerly_grazed = factor(
      formerly_grazed, levels = c(FALSE, TRUE),
      labels = c('ungrazed', 'formerly_grazed')),
    fire_season = factor(
      fire_season,
      levels = c('Spring', 'Late spring', 'Summer', 'Fall')),
    postfire_period = factor(
      postfire_period, levels = c('0', '1', '2+')))


# Fire-exclusion check ---------------------------------------------------------
# Raw measurements less than four months postfire are retained in df_year but
# are flagged analysis_usable = FALSE, matching the exclusion in the paper.
df_recent_fire_check <- df_year %>%
  filter(recent_fire_lt4mo) %>%
  select(
    quad, year, census_month, fire_year, fire_month,
    months_since_fire, st, fl, mht, analysis_usable)

df_recent_fire_check


# Annual quadrat transitions ---------------------------------------------------
# These are transitions in quadrat-level ramet abundance and mean height, not
# individual survival/growth transitions. Individual ramets were not followed.

df_transition <- df_year %>%
  arrange(quad, year) %>%
  group_by(quad) %>%
  mutate(
    year_t1 = lead(year),
    census_date_t1 = lead(census_date),
    annual_transition = coalesce(year_t1 == year + 1, FALSE),
    st_t1 = if_else(annual_transition, lead(st), NA_real_),
    fl_t1 = if_else(annual_transition, lead(fl), NA_real_),
    mht_t1 = if_else(annual_transition, lead(mht), NA_real_),
    analysis_usable_t1 = if_else(
      annual_transition, lead(analysis_usable), FALSE),
    recent_fire_t1 = if_else(
      annual_transition, lead(recent_fire_lt4mo), FALSE),
    months_since_fire_t1 = if_else(
      annual_transition, lead(months_since_fire), NA_real_),
    years_since_fire_t1 = if_else(
      annual_transition, lead(years_since_fire), NA_real_),
    transition_usable = annual_transition &
      !is.na(st) & !is.na(st_t1),
    transition_usable_strict = transition_usable &
      analysis_usable & analysis_usable_t1) %>%
  ungroup()


# Fires between consecutive censuses ------------------------------------------
df_fire_between <- df_transition %>%
  filter(annual_transition) %>%
  select(quad, year, census_date, census_date_t1) %>%
  mutate(quad_num = as.numeric(as.character(quad))) %>%
  inner_join(
    df_paper_fire,
    by = c('quad_num' = 'quad'),
    relationship = 'many-to-many') %>%
  filter(
    fire_date > census_date,
    fire_date < census_date_t1) %>%
  group_by(quad, year) %>%
  summarise(
    fire_between = 1,
    fire_between_date = max(fire_date),
    fire_between_year = fire_year[which.max(fire_date)],
    fire_between_month = fire_month[which.max(fire_date)],
    .groups = 'drop')

df_transition <- df_transition %>%
  left_join(df_fire_between, by = c('quad', 'year')) %>%
  mutate(
    fire_between = replace_na(fire_between, 0),
    stem_change = if_else(
      transition_usable_strict, st_t1 - st, NA_real_),
    log_stems_t0 = if_else(
      transition_usable_strict, log1p(st), NA_real_),
    log_stems_t1 = if_else(
      transition_usable_strict, log1p(st_t1), NA_real_),
    log_stem_change = log_stems_t1 - log_stems_t0)


# Critical fire-transition diagnostic -----------------------------------------
# With the <4 month postfire exclusion, no directly observed prefire-to-postfire
# annual transition is usable for estimating a discrete fire-transition effect.
df_transition_fire_check <- df_transition %>%
  filter(annual_transition) %>%
  count(
    fire_between, transition_usable,
    transition_usable_strict, name = 'n')

df_transition_fire_check

df_fire_crossing <- df_transition %>%
  filter(fire_between == 1) %>%
  select(
    quad, year, year_t1, st, st_t1,
    fire_between_year, fire_between_month,
    months_since_fire_t1, recent_fire_t1,
    transition_usable, transition_usable_strict)

df_fire_crossing


# Final workdata ---------------------------------------------------------------
# One row = one quadrat-year. Transition variables describe the next annual
# census where available. Fire effects should be modeled using time since fire,
# not as an estimable binary prefire/postfire transition effect.

df_transition_keep <- df_transition %>%
  select(
    quad, year, year_t1, annual_transition,
    st_t1, fl_t1, mht_t1,
    analysis_usable_t1, recent_fire_t1,
    months_since_fire_t1, years_since_fire_t1,
    fire_between, fire_between_date,
    fire_between_year, fire_between_month,
    transition_usable, transition_usable_strict,
    stem_change, log_stems_t0, log_stems_t1,
    log_stem_change)

df <- df_year %>%
  left_join(df_transition_keep, by = c('quad', 'year')) %>%
  arrange(quad, year)


# Final checks ----------------------------------------------------------------
stopifnot(
  nrow(df_og) == 230,
  n_distinct(df_og$quad) == 23,
  n_distinct(df_og$year) == 10,
  nrow(df_duplicates) == 0,
  nrow(df_missing_rows) == 0,
  all(df_value_checks$n == 0),
  nrow(df_fire_consistency) == 0,
  all(df_fire_quad$nfires_difference == 0),
  all(df_fire_csv_vs_paper$matches_paper),
  max(abs(
    df_paper_mean_check$stem_mean_difference),
    na.rm = TRUE) < 0.1,
  max(abs(
    df_paper_mean_check$height_mean_difference),
    na.rm = TRUE) < 0.1,
  nrow(df) == 230,
  nrow(df) == nrow(distinct(df, quad, year)),
  sum(
    df_transition$fire_between == 1 &
      df_transition$transition_usable_strict,
    na.rm = TRUE) == 0)


# Workdata summaries -----------------------------------------------------------
df_work_summary <- df %>%
  summarise(
    n_rows = n(),
    n_quads = n_distinct(quad),
    n_analysis_usable = sum(analysis_usable),
    n_postfire_usable = sum(postfire_usable),
    n_height_usable = sum(height_usable),
    n_flower_usable = sum(flower_usable),
    n_annual_transitions = sum(annual_transition),
    n_transition_usable = sum(transition_usable),
    n_transition_usable_strict =
      sum(transition_usable_strict),
    n_fire_crossing = sum(fire_between == 1, na.rm = TRUE),
    n_fire_crossing_strict = sum(
      fire_between == 1 & transition_usable_strict,
      na.rm = TRUE))

df_work_summary


# Postfire summary -------------------------------------------------------------
df_postfire_summary <- df %>%
  filter(postfire_usable) %>%
  group_by(postfire_year, formerly_grazed) %>%
  summarise(
    n = n(),
    mean_st = mean(st),
    median_st = median(st),
    flower_prop = if_else(
      sum(st) > 0, sum(fl) / sum(st), NA_real_),
    mean_mht = mean(mht, na.rm = TRUE),
    .groups = 'drop')

df_postfire_summary


# Sampling inventory -----------------------------------------------------------
fig_sampling <- df %>%
  mutate(
    quad_plot = fct_rev(quad),
    status = case_when(
      analysis_usable ~ 'usable',
      recent_fire_lt4mo ~ '<4 months postfire',
      TRUE ~ 'not observed')) %>%
  ggplot(aes(x = year, y = quad_plot, fill = status)) +
  geom_tile() +
  scale_x_continuous(breaks = 1991:2000) +
  theme_bw() +
  labs(
    title = 'Sampling inventory',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = 'Quadrat',
    fill = 'Status') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank())

fig_sampling


# Stem trajectories ------------------------------------------------------------
fig_stems_time <- df %>%
  ggplot(aes(x = year, y = st, group = quad)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  facet_wrap(~quad, scales = 'free_y') +
  geom_vline(
    data = df_paper_fire,
    aes(xintercept = fire_year),
    linetype = 2, inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(1991, 1995, 2000)) +
  theme_bw() +
  labs(
    title = 'Stem abundance through time',
    subtitle = v_ggp_suffix,
    x = 'Year',
    y = 'Number of stems per 1 m2') +
  theme(
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 6))

fig_stems_time


# Postfire stem abundance ------------------------------------------------------
fig_stems_fire <- df %>%
  filter(postfire_usable) %>%
  ggplot(
    aes(
      x = years_since_fire,
      y = log1p(st),
      shape = formerly_grazed)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_bw() +
  labs(
    title = 'Stem abundance vs time since fire',
    subtitle = v_ggp_suffix,
    x = 'Years since most recent fire',
    y = 'log(1 + stems per 1 m2)',
    shape = 'Site history')

fig_stems_fire


# Postfire flowering -----------------------------------------------------------
fig_flower_fire <- df %>%
  filter(flower_usable, !is.na(years_since_fire)) %>%
  ggplot(
    aes(
      x = years_since_fire,
      y = flower_prop,
      shape = formerly_grazed)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(
    title = 'Flowering vs time since fire',
    subtitle = v_ggp_suffix,
    x = 'Years since most recent fire',
    y = 'Proportion of stems flowering',
    shape = 'Site history')

fig_flower_fire


# Postfire mean height ---------------------------------------------------------
fig_height_fire <- df %>%
  filter(height_usable, !is.na(years_since_fire)) %>%
  ggplot(
    aes(
      x = years_since_fire,
      y = mht,
      shape = fire_season)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(
    title = 'Mean ramet height vs time since fire',
    subtitle = v_ggp_suffix,
    x = 'Years since most recent fire',
    y = 'Mean ramet height (cm)',
    shape = 'Fire season')

fig_height_fire


# Annual abundance transitions ------------------------------------------------
fig_stem_transition <- df %>%
  filter(transition_usable_strict) %>%
  ggplot(aes(x = log_stems_t0, y = log_stems_t1)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    title = 'Annual quadrat-level stem transitions',
    subtitle = v_ggp_suffix,
    x = 'log(1 + stems) at t',
    y = 'log(1 + stems) at t+1')

fig_stem_transition


# Save data -------------------------------------------------------------------
# write.csv(
#   df, row.names = FALSE,
#   file.path(
#     dir_data,
#     paste0('ab_', v_sp_abb, '_df_workdata_260831.csv')))

