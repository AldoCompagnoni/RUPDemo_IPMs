# Data - Archbold - Eriogonum longifolium

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.08.19

# Study organism: Eriogonum longifolium var. gnaphalifolium
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.226.1
# Meta data link:
# https://portal.edirepository.org/nis/metadataviewer?packageid=edi.226.1
# Citing publication: Satterthwaite et al. 2002, Ecological Applications
# Time period: 1990-2013


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
load_packages(MASS, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4,
              GGally, tidyverse)


# Specification ----------------------------------------------------------------
# Define head-directory
v_head <- c('archbold')
# Define species
v_species <- c('Eriogonum longifolium')
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c()

# Create a unique species abbreviation for file naming
v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Define script prefix
v_script_prefix <- str_c(v_head)

# Plot subtitle
v_ggp_suffix <- paste(
  tools::toTitleCase(v_head), '-', v_species)


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


# Data -------------------------------------------------------------------------
df_og <- read_csv(
  file.path(dir_data, 'eriogonum_longifolium_data.csv'),
  col_types = cols(comment = col_character())) %>%
  janitor::clean_names()

# Keep the metadata table open until the EDI code definitions are confirmed.
df_meta <- tibble(
  var = names(df_og),
  description = NA_character_)

# The plant number is not a globally unique identifier in this dataset.
# Use site + population + quadrat + plant as the provisional individual ID.
df_og <- df_og %>%
  mutate(
    year = as.numeric(str_sub(date, 1, 4)),
    month = as.numeric(str_sub(date, 6, 7)),
    id = str_c(site, pop, qu, plant, sep = '_'),
    record_type = case_when(
      !is.na(s) ~ 'demography',
      is.na(s) & !is.na(burn) ~ 'burn',
      TRUE ~ 'other'))


# Initial data structure --------------------------------------------------------
dim(df_og)
names(df_og)
skimr::skim(df_og)

df_og %>%
  count(record_type)

df_og %>%
  count(year, month)

df_og %>%
  count(site, pop)

df_og %>%
  group_by(site, pop) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    nr_quads = n_distinct(qu),
    nr_plants = n_distinct(id), .groups = 'drop')


# Separate demographic and burn records ---------------------------------------
# Burn observations are stored as separate rows rather than in the annual
# demographic record. Keep them separate until the fire codes and timing are
# confirmed from the metadata.
df_demog <- df_og %>%
  filter(record_type == 'demography')

df_burn <- df_og %>%
  filter(record_type == 'burn')

# Candidate annual census records. June occurs in every year from 1990-2013,
# but extra demographic censuses also occur in some years and are retained in
# df_demog for now.
df_demog_june <- df_demog %>%
  filter(month == 6)


# Sampling inventory -----------------------------------------------------------
fig_sampling <- df_demog_june %>%
  distinct(site, pop, qu, year) %>%
  mutate(pop_quad = interaction(pop, qu, sep = '-')) %>%
  ggplot() +
  geom_point(aes(x = year, y = pop_quad)) +
  facet_wrap(~site, scales = 'free_y') +
  theme_bw() +
  labs(title = 'Sampling inventory', x = 'Year', y = 'Population-quadrat') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))
fig_sampling

# ggsave(file.path(dir_result, '0_sampling_inventory.png'),
#        plot = fig_sampling, width = 10, height = 7, dpi = 300)


# Exploring identifiers and duplicates ----------------------------------------
# Plant tag numbers are reused among quadrats/populations.
df_og %>%
  distinct(site, pop, qu, plant) %>%
  count(plant, name = 'nr_ids') %>%
  filter(nr_ids > 1) %>%
  arrange(desc(nr_ids), plant)

# Check for duplicate demographic records within individual and date.
df_demog %>%
  group_by(id, date) %>%
  mutate(n_entries = n()) %>%
  ungroup() %>%
  filter(n_entries > 1) %>%
  arrange(id, date)

# Check for duplicate burn records within individual and date.
df_burn %>%
  group_by(id, date) %>%
  mutate(n_entries = n()) %>%
  ungroup() %>%
  filter(n_entries > 1) %>%
  arrange(id, date)


# Explore coding before building transitions ----------------------------------
# Do not convert these codes into survival, recruitment, dormancy or fire
# variables until their exact definitions have been confirmed.
df_demog %>%
  count(s, stage, sort = TRUE)

df_demog %>%
  count(stage, sort = TRUE)

df_demog %>%
  count(s, sort = TRUE)

df_burn %>%
  count(date, site, pop, burn, sort = TRUE)

df_burn %>%
  group_by(date, site, pop, qu) %>%
  summarise(
    nr_records = n(),
    nr_burn_codes = n_distinct(burn),
    burn_codes = paste(sort(unique(burn)), collapse = ', '),
    .groups = 'drop') %>%
  arrange(date, site, pop, qu)



# Clean annual demographic history --------------------------------------------
# Annual demographic records are the June censuses. Plants recorded as absent
# or previously absent are back-corrected to dormant when they occur between
# two observations where the individual is known to be alive.

df_alive_range <- df_demog_june %>%
  filter(s %in% c(1, 3, 5)) %>%
  group_by(id) %>%
  summarise(
    first_alive = min(year),
    last_alive = max(year), .groups = "drop")

df_annual <- df_demog_june %>%
  left_join(df_alive_range, by = "id") %>%
  mutate(
    dormancy_repair = s %in% c(0, 9) &
      !is.na(first_alive) & year > first_alive & year < last_alive,
    s_clean = if_else(dormancy_repair, 8, s),
    stage_clean = if_else(dormancy_repair, 5, stage),
    state_clean = case_when(
      s_clean %in% c(1, 3, 5) ~ "active",
      s_clean == 8 ~ "dormant",
      s_clean == 10 ~ "alive_or_dormant",
      s_clean == 0 ~ "absent",
      s_clean == 9 ~ "previous_absent",
      s_clean %in% c(2, 6, 7) ~ "missing",
      s_clean == 11 ~ "discontinued",
      TRUE ~ NA_character_))

# Check repaired dormancy records
df_annual %>%
  filter(dormancy_repair) %>%
  count(s, s_clean)

df_dormancy_check <- df_annual %>%
  filter(any(dormancy_repair), .by = id) %>%
  select(
    site, pop, qu, plant, id, year, s, stage, s_clean, stage_clean,
    state_clean, dia, scape, comment) %>%
  arrange(id, year)

df_dormancy_check


# Dormancy duration ------------------------------------------------------------
# Each continuous dormant period is treated as one independent dormancy spell.
# Individuals can therefore contribute more than one dormancy period.

df_dormancy_spell <- df_annual %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    new_spell = row_number() == 1 |
      coalesce(state_clean != lag(state_clean), TRUE) |
      coalesce(year != lag(year) + 1, TRUE),
    spell = cumsum(new_spell)) %>%
  group_by(id, spell) %>%
  summarise(
    state = first(state_clean),
    start_year = min(year),
    end_year = max(year),
    dormancy_years = n(), .groups = "drop") %>%
  filter(state == "dormant")


# Dormancy summary -------------------------------------------------------------
df_dormancy_summary <- df_dormancy_spell %>%
  summarise(
    n_spells = n(),
    n_individuals = n_distinct(id),
    mean_dormancy = mean(dormancy_years),
    median_dormancy = median(dormancy_years),
    sd_dormancy = sd(dormancy_years),
    min_dormancy = min(dormancy_years),
    max_dormancy = max(dormancy_years))

df_dormancy_summary


# Dormancy duration distribution -----------------------------------------------
df_dormancy_distribution <- df_dormancy_spell %>%
  count(dormancy_years, name = "n_spells") %>%
  arrange(dormancy_years) %>%
  mutate(
    prop = n_spells / sum(n_spells),
    cum_prop = cumsum(prop),
    cum_percent = 100 * cum_prop)

df_dormancy_distribution


# Three-year dormancy cutoff ---------------------------------------------------
# The 3-year cutoff is an operational assumption, 
# justified because 388/409 = 94.87% of observed dormancy spells lasted ≤3 years

dormancy_cutoff <- 3

df_dormancy_cutoff <- df_dormancy_spell %>%
  summarise(
    cutoff_years = dormancy_cutoff,
    n_spells = n(),
    n_within_cutoff = sum(dormancy_years <= dormancy_cutoff),
    prop_within_cutoff =
      mean(dormancy_years <= dormancy_cutoff),
    percent_within_cutoff =
      100 * mean(dormancy_years <= dormancy_cutoff))

df_dormancy_cutoff


# Annual demographic transitions ----------------------------------------------
df_transition <- df_annual %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    year_t1 = lead(year),
    annual_transition = coalesce(year_t1 == year + 1, FALSE),
    s_t1 = if_else(annual_transition, lead(s_clean), NA_real_),
    stage_t1 = if_else(annual_transition, lead(stage_clean), NA_real_),
    state_t1 = if_else(
      annual_transition, lead(state_clean), NA_character_),
    size_t0 = if_else(state_clean == "active", dia, NA_real_),
    size_t1 = if_else(
      annual_transition & state_clean == "active" &
        lead(state_clean) == "active",
      lead(dia), NA_real_)) %>%
  ungroup()


# Check annual state transitions -----------------------------------------------
df_transition %>%
  filter(annual_transition) %>%
  count(state_clean, state_t1) |> 
  print(n=100)

df_transition %>%
  filter(
    annual_transition,
    state_clean == "active",
    state_t1 != "active") %>%
  count(year, state_t1) |> 
  print(n=100)


# Unresolved absence periods ---------------------------------------------------
df_absence_spell <- df_annual %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    absent_state =
      state_clean %in% c("absent", "previous_absent"),
    new_spell = row_number() == 1 |
      coalesce(absent_state != lag(absent_state), TRUE) |
      coalesce(year != lag(year) + 1, TRUE),
    spell = cumsum(new_spell)) %>%
  group_by(id, spell) %>%
  summarise(
    absent_state = first(absent_state),
    absence_start = min(year),
    absence_end = max(year),
    absence_years = n(), .groups = "drop") %>%
  filter(absent_state)

df_absence_t1 <- df_absence_spell %>%
  transmute(
    id,
    year_t1 = absence_start,
    absence_years)


# Demographic outcomes ---------------------------------------------------------
df_transition <- df_transition %>%
  left_join(df_absence_t1, by = c("id", "year_t1")) %>%
  mutate(
    survives = case_when(
      !annual_transition ~ NA_real_,
      !state_clean %in% c("active", "dormant") ~ NA_real_,
      state_t1 %in%
        c("active", "dormant", "alive_or_dormant") ~ 1,
      state_t1 %in% c("absent", "previous_absent") &
        absence_years > dormancy_cutoff ~ 0,
      state_t1 %in% c("absent", "previous_absent") &
        absence_years <= dormancy_cutoff ~ NA_real_,
      TRUE ~ NA_real_))

df_transition %>%
  filter(state_clean %in% c("active", "dormant")) %>%
  count(state_clean, survives)

df_transition %>%
  filter(state_clean %in% c("active", "dormant")) %>%
  count(state_clean, survives)


# Active and dormant transitions -----------------------------------------------
df_transition <- df_transition %>%
  mutate(
    enter_dormancy = case_when(
      state_clean != "active" | survives != 1 ~ NA_real_,
      state_t1 == "dormant" ~ 1,
      state_t1 == "active" ~ 0,
      TRUE ~ NA_real_),
    
    reactivate = case_when(
      state_clean != "dormant" ~ NA_real_,
      state_t1 == "active" ~ 1,
      state_t1 == "dormant" ~ 0,
      TRUE ~ NA_real_))

df_transition %>%
  count(enter_dormancy)

df_transition %>%
  count(reactivate)


# Check flowering classification -----------------------------------------------
df_annual %>%
  filter(state_clean == "active") %>%
  count(stage_clean, scape > 0, is.na(scape))

df_annual %>%
  filter(
    state_clean == "active",
    (stage_clean == 4 & (scape == 0 | is.na(scape))) |
      (stage_clean != 4 & scape > 0)) %>%
  select(
    site, pop, qu, plant, id, year, s_clean, stage_clean, dia, scape) %>%
  arrange(id, year)

df_transition %>%
  summarise(
    n_active_active = sum(
      state_clean == "active" & state_t1 == "active", na.rm = TRUE),
    n_growth = sum(
      state_clean == "active" & state_t1 == "active" &
        !is.na(size_t0) & !is.na(size_t1), na.rm = TRUE),
    n_missing_growth = sum(
      state_clean == "active" & state_t1 == "active" &
        (is.na(size_t0) | is.na(size_t1)), na.rm = TRUE))


# Growth and flowering ---------------------------------------------------------
df_transition <- df_transition %>%
  mutate(
    flower = case_when(
      state_clean != "active" ~ NA_real_,
      !is.na(scape) & scape > 0 ~ 1,
      !is.na(scape) & scape == 0 ~ 0,
      TRUE ~ NA_real_),
    fl_nr = if_else(state_clean == "active", scape, NA_real_),
    logsize_t0 = log(size_t0),
    logsize_t1 = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3)

df_transition %>%
  filter(
    state_clean == "active",
    state_t1 == "active",
    !is.na(size_t0),
    !is.na(size_t1)) %>%
  summarise(
    n = n(),
    min_t0 = min(size_t0),
    max_t0 = max(size_t0),
    min_t1 = min(size_t1),
    max_t1 = max(size_t1),
    mean_t0 = mean(size_t0),
    mean_t1 = mean(size_t1))

df_transition %>%
  filter(state_clean == "active") %>%
  count(flower)


# Recruitment structure --------------------------------------------------------
df_quad_start <- df_demog_june %>%
  group_by(site, pop, qu) %>%
  summarise(first_quad_year = min(year), .groups = "drop")

df_recruit_first <- df_demog %>%
  filter(s %in% c(3, 5)) %>%
  arrange(id, year, month) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    recruit_year = if_else(month > 6, year + 1, year)) %>%
  left_join(
    df_quad_start,
    by = c("site", "pop", "qu")) %>%
  mutate(
    baseline = recruit_year == first_quad_year,
    recruit_type = case_when(
      s == 5 ~ "seedling",
      s == 3 ~ "new_adult"))

df_recruit_first %>%
  count(baseline, recruit_type)

df_recruit_first %>%
  filter(!baseline) %>%
  count(recruit_year, recruit_type)

df_recruit_first %>%
  filter(!baseline) %>%
  count(month, recruit_type)


# Recruitment size structure ---------------------------------------------------
df_recruit_first %>%
  filter(!baseline) %>%
  mutate(
    size_class = case_when(
      is.na(dia) ~ "missing",
      dia <= 2 ~ "<=2 cm",
      dia > 2 ~ ">2 cm")) %>%
  count(recruit_type, size_class)

df_recruit_first %>%
  filter(!baseline) %>%
  group_by(recruit_type, month) %>%
  summarise(
    n = n(),
    n_dia_na = sum(is.na(dia)),
    mean_dia = mean(dia, na.rm = TRUE),
    median_dia = median(dia, na.rm = TRUE),
    min_dia = if (all(is.na(dia))) NA_real_ else min(dia, na.rm = TRUE),
    max_dia = if (all(is.na(dia))) NA_real_ else max(dia, na.rm = TRUE),
    .groups = "drop")


# Recruitment cohorts ----------------------------------------------------------
# Seedlings enter in their first annual census. New adults are treated as
# presumed 2-year-olds and assigned to the preceding seedling cohort.
df_recruit_cohort <- df_recruit_first %>%
  filter(!baseline) %>%
  rename(annual_entry_year = recruit_year) %>%
  mutate(
    recruit_year = case_when(
      recruit_type == "seedling" ~ annual_entry_year,
      recruit_type == "new_adult" ~ annual_entry_year - 1))

df_transition %>%
  filter(stage_clean == 1) %>%
  count(survives)

df_transition %>%
  filter(stage_clean == 1, !is.na(survives)) %>%
  summarise(
    n = n(),
    seedling_survival = mean(survives))

df_recruit_first %>%
  filter(
    !baseline,
    (recruit_type == "seedling" & dia > 2) |
      (recruit_type == "new_adult" & dia <= 2)) %>%
  select(
    site, pop, qu, plant, id, date, s, stage, dia,
    recruit_type, comment)

# Recruitment classification --------------------------------------------------
df_recruit_entry <- df_recruit_first %>%
  filter(!baseline) %>%
  left_join(
    df_demog_june %>%
      select(
        id, year, s_entry = s, stage_entry = stage,
        dia_entry = dia),
    by = c("id", "recruit_year" = "year")) %>%
  mutate(
    recruit_class = case_when(
      recruit_type == "seedling" ~ "observed_seedling",
      recruit_type == "new_adult" & stage_entry == 2 ~ "presumed_2yr",
      recruit_type == "new_adult" ~ "new_adult_other"),
    fecundity_year = case_when(
      recruit_class == "observed_seedling" ~ recruit_year - 1,
      recruit_class == "presumed_2yr" ~ recruit_year - 2,
      TRUE ~ NA_real_))

df_recruit_entry %>%
  count(recruit_class)

df_recruit_entry %>%
  filter(!is.na(fecundity_year)) %>%
  anti_join(
    df_demog_june %>%
      distinct(site, pop, qu, year),
    by = c(
      "site", "pop", "qu",
      "fecundity_year" = "year")) %>%
  count(recruit_class, fecundity_year)


# Recruitment observations usable for fecundity -------------------------------
df_recruit_entry <- df_recruit_entry %>%
  left_join(
    df_demog_june %>%
      distinct(site, pop, qu, year) %>%
      mutate(fecundity_observed = TRUE),
    by = c(
      "site", "pop", "qu",
      "fecundity_year" = "year")) %>%
  mutate(
    fecundity_usable = case_when(
      recruit_class == "new_adult_other" ~ FALSE,
      !is.na(fecundity_year) &
        replace_na(fecundity_observed, FALSE) ~ TRUE,
      TRUE ~ FALSE))

df_recruit_entry %>%
  count(recruit_class, fecundity_usable)


# Quadrat monitoring history --------------------------------------------------
df_quad_monitor <- df_annual %>%
  group_by(site, pop, qu, year) %>%
  summarise(
    monitored = any(state_clean != "discontinued", na.rm = TRUE),
    .groups = "drop") %>%
  filter(monitored) %>%
  select(-monitored)


# Recruitment follow-up availability ------------------------------------------
df_recruit_followup <- df_quad_monitor %>%
  left_join(
    df_quad_monitor %>%
      transmute(
        site, pop, qu,
        year = year - 1,
        monitored_t1 = TRUE),
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_quad_monitor %>%
      transmute(
        site, pop, qu,
        year = year - 2,
        monitored_t2 = TRUE),
    by = c("site", "pop", "qu", "year")) %>%
  mutate(
    monitored_t1 = replace_na(monitored_t1, FALSE),
    monitored_t2 = replace_na(monitored_t2, FALSE))


# Recruitment counts by fecundity year ----------------------------------------
df_recruit_count <- df_recruit_entry %>%
  filter(fecundity_usable) %>%
  count(
    site, pop, qu, fecundity_year, recruit_class,
    name = "nr_recruits") %>%
  pivot_wider(
    names_from = recruit_class,
    values_from = nr_recruits,
    values_fill = 0)


# Recruitment working structure -----------------------------------------------
df_recruit_quad <- df_recruit_followup %>%
  left_join(
    df_recruit_count,
    by = c(
      "site", "pop", "qu",
      "year" = "fecundity_year")) %>%
  mutate(
    observed_seedling = case_when(
      monitored_t1 ~ replace_na(observed_seedling, 0L),
      TRUE ~ NA_integer_),
    presumed_2yr = case_when(
      monitored_t2 ~ replace_na(presumed_2yr, 0L),
      TRUE ~ NA_integer_),
    recruitment_complete = monitored_t1 & monitored_t2,
    nr_recruit = case_when(
      recruitment_complete ~ observed_seedling + presumed_2yr,
      TRUE ~ NA_integer_))

df_recruit_quad %>%
  count(year, monitored_t1, monitored_t2) %>%
  print(n = 100)

df_recruit_quad %>%
  summarise(
    n_quad_years = n(),
    n_complete = sum(recruitment_complete),
    n_incomplete = sum(!recruitment_complete),
    n_seedling_unknown = sum(is.na(observed_seedling)),
    n_presumed_unknown = sum(is.na(presumed_2yr)))

# Parental reproductive abundance ---------------------------------------------
df_repro_quad <- df_quad_monitor %>%
  left_join(
    df_annual %>%
      filter(state_clean != "discontinued") %>%
      group_by(site, pop, qu, year) %>%
      summarise(
        nr_active = sum(state_clean == "active", na.rm = TRUE),
        nr_flower_obs = sum(
          state_clean == "active" & !is.na(scape), na.rm = TRUE),
        flowering_sum = sum(
          state_clean == "active" & scape > 0, na.rm = TRUE),
        scape_sum = sum(
          if_else(state_clean == "active", scape, NA_real_),
          na.rm = TRUE),
        nr_invol_obs = sum(
          state_clean == "active" & !is.na(invol), na.rm = TRUE),
        invol_sum = sum(
          if_else(state_clean == "active", invol, NA_real_),
          na.rm = TRUE), .groups = "drop") %>%
      mutate(
        nr_flowering = case_when(
          nr_active == 0 ~ 0,
          nr_flower_obs == 0 ~ NA_real_,
          TRUE ~ as.numeric(flowering_sum)),
        nr_scapes = case_when(
          nr_active == 0 ~ 0,
          nr_flower_obs == 0 ~ NA_real_,
          TRUE ~ scape_sum),
        nr_invol = case_when(
          nr_active == 0 ~ 0,
          nr_invol_obs == 0 ~ NA_real_,
          TRUE ~ invol_sum)) %>%
      select(
        site, pop, qu, year, nr_active, nr_flower_obs,
        nr_flowering, nr_scapes, nr_invol_obs, nr_invol),
    by = c("site", "pop", "qu", "year"))


# Population-level reproductive abundance -------------------------------------
df_repro_pop <- df_repro_quad %>%
  group_by(site, pop, year) %>%
  summarise(
    nr_quads = n(),
    nr_active_pop = sum(nr_active),
    nr_quads_flower_obs = sum(!is.na(nr_flowering)),
    nr_flowering_pop_obs = sum(nr_flowering, na.rm = TRUE),
    nr_quads_scape_obs = sum(!is.na(nr_scapes)),
    nr_scapes_pop_obs = sum(nr_scapes, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    nr_flowering_pop = if_else(
      nr_quads_flower_obs == nr_quads,
      nr_flowering_pop_obs, NA_real_),
    nr_scapes_pop = if_else(
      nr_quads_scape_obs == nr_quads,
      nr_scapes_pop_obs, NA_real_),
    mean_flowering_pop = if_else(
      nr_quads_flower_obs > 0,
      nr_flowering_pop_obs / nr_quads_flower_obs, NA_real_),
    mean_scapes_pop = if_else(
      nr_quads_scape_obs > 0,
      nr_scapes_pop_obs / nr_quads_scape_obs, NA_real_))


# Add reproductive abundance to recruitment data ------------------------------
df_recruit_quad <- df_recruit_quad %>%
  left_join(
    df_repro_quad,
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_repro_pop,
    by = c("site", "pop", "year")) %>%
  mutate(
    nr_flowering_other =
      nr_flowering_pop_obs - replace_na(nr_flowering, 0),
    nr_quads_flower_other =
      nr_quads_flower_obs - as.integer(!is.na(nr_flowering)),
    mean_flowering_other = if_else(
      nr_quads_flower_other > 0,
      nr_flowering_other / nr_quads_flower_other, NA_real_),
    nr_scapes_other =
      nr_scapes_pop_obs - replace_na(nr_scapes, 0),
    nr_quads_scape_other =
      nr_quads_scape_obs - as.integer(!is.na(nr_scapes)),
    mean_scapes_other = if_else(
      nr_quads_scape_other > 0,
      nr_scapes_other / nr_quads_scape_other, NA_real_))


# Recruitment without local flowering -----------------------------------------
df_recruit_quad %>%
  mutate(nr_recruit = observed_seedling + presumed_2yr) %>%
  filter(nr_recruit > 0) %>%
  count(local_flowering = nr_flowering > 0)

# Reproductive data through time -----------------------------------------------
df_repro_quad %>%
  group_by(year) %>%
  summarise(
    nr_flowering = sum(nr_flowering),
    nr_scapes = sum(nr_scapes, na.rm = TRUE),
    nr_invol_obs = sum(nr_invol_obs),
    .groups = "drop")

# Recruitment and flowering at population level -------------------------------
df_recruit_pop <- df_recruit_quad %>%
  group_by(site, pop, year) %>%
  summarise(
    observed_seedling = sum(observed_seedling),
    presumed_2yr = sum(presumed_2yr),
    nr_flowering = sum(nr_flowering),
    nr_scapes = sum(nr_scapes, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(nr_recruit = observed_seedling + presumed_2yr)

df_recruit_pop %>%
  filter(nr_recruit > 0, nr_flowering == 0)


# Fire-event structure ---------------------------------------------------------
df_fire <- df_burn %>%
  mutate(
    fire_year = year,
    fire_month = month,
    fire_event = date,
    burned = case_when(
      burn == 0 ~ 0,
      burn %in% c(1, 2, 4) ~ 1,
      TRUE ~ NA_real_),
    fire_severity = case_when(
      burn == 0 ~ "unburned",
      burn == 1 ~ "scorched",
      burn == 2 ~ "consumed",
      burn == 4 ~ "burned_unspecified",
      TRUE ~ NA_character_),
    transition_year = case_when(
      fire_event == "1998-05" ~ 1997,
      fire_event == "1998-07" ~ 1998,
      fire_event == "2001-08" ~ 2001,
      fire_event == "2002-10" ~ 2002,
      fire_event == "2005-04" ~ 2004,
      fire_event == "2007-12" ~ 2007,
      TRUE ~ NA_real_),
    timing_known = !is.na(transition_year))

df_fire %>%
  count(fire_event, fire_severity)

df_fire %>%
  group_by(fire_event, site, pop, qu) %>%
  summarise(
    nr_scored = n(),
    nr_burned = sum(burned == 1, na.rm = TRUE),
    prop_burned = mean(burned, na.rm = TRUE),
    .groups = "drop") %>%
  group_by(fire_event) %>%
  summarise(
    nr_quads = n(),
    min_prop = min(prop_burned),
    mean_prop = mean(prop_burned),
    max_prop = max(prop_burned),
    .groups = "drop")


# Fire records relative to demographic observations ---------------------------
df_first_demog <- df_demog %>%
  group_by(id) %>%
  summarise(
    first_demog_date = min(date),
    .groups = "drop")

df_fire_history_check <- df_fire %>%
  left_join(df_first_demog, by = "id") %>%
  mutate(
    fire_time = as.Date(paste0(fire_event, "-01")),
    first_demog_time = as.Date(paste0(first_demog_date, "-01")),
    first_observed_after_fire = first_demog_time > fire_time)

df_fire_history_check %>%
  count(fire_event, first_observed_after_fire)

df_fire_history_check %>%
  filter(fire_event == "2001-08") %>%
  count(first_demog_date)

# Quadrat-level fire exposure --------------------------------------------------
df_fire_quad <- df_fire %>%
  group_by(fire_event, site, pop, qu) %>%
  summarise(
    transition_year = first(transition_year),
    timing_known = first(timing_known),
    nr_scored = n(),
    nr_burned = sum(burned == 1, na.rm = TRUE),
    nr_unburned = sum(burned == 0, na.rm = TRUE),
    prop_burned = mean(burned, na.rm = TRUE),
    any_burned = as.numeric(any(burned == 1, na.rm = TRUE)),
    .groups = "drop")


# Fire-event spatial footprint -------------------------------------------------
df_fire_quad %>%
  distinct(fire_event, site, pop, timing_known, transition_year) %>%
  arrange(fire_event, site, pop)

# Quadrat coverage of known-timing fires --------------------------------------
df_fire_pop <- df_fire_quad %>%
  distinct(
    fire_event, site, pop, transition_year, timing_known)

df_fire_quad_coverage <- df_fire_pop %>%
  filter(timing_known) %>%
  inner_join(
    df_demog_june %>%
      distinct(site, pop, qu, year),
    by = c(
      "site", "pop",
      "transition_year" = "year")) %>%
  left_join(
    df_fire_quad %>%
      select(
        fire_event, site, pop, qu,
        prop_burned, any_burned),
    by = c("fire_event", "site", "pop", "qu")) %>%
  mutate(fire_scored = !is.na(prop_burned))

df_fire_quad_coverage %>%
  group_by(fire_event, site, pop) %>%
  summarise(
    nr_quads_monitored = n(),
    nr_quads_scored = sum(fire_scored),
    prop_quads_scored = mean(fire_scored),
    .groups = "drop")

df_fire_quad_coverage %>%
  filter(!fire_scored) %>%
  count(fire_event, site, pop)


# Known fire events ------------------------------------------------------------
df_fire_pop_year <- df_fire_quad %>%
  filter(timing_known) %>%
  distinct(fire_event, site, pop, transition_year)


# Ambiguous fire timing --------------------------------------------------------
# For these fires, the annual census could have occurred either before or after
# the fire. Both possible annual transitions are therefore treated as unknown.

df_fire_ambiguous <- df_fire_quad %>%
  filter(!timing_known) %>%
  distinct(fire_event, site, pop) %>%
  mutate(fire_year = as.numeric(str_sub(fire_event, 1, 4))) %>%
  {
    bind_rows(
      transmute(
        ., fire_event, site, pop,
        candidate_year = fire_year - 1),
      transmute(
        ., fire_event, site, pop,
        candidate_year = fire_year))
  }


# Annual quadrat fire exposure -------------------------------------------------
df_fire_transition_quad <- df_demog_june %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_fire_pop_year,
    by = c("site", "pop", "year" = "transition_year")) %>%
  left_join(
    df_fire_quad %>%
      filter(timing_known) %>%
      select(
        fire_event, site, pop, qu,
        nr_scored, nr_burned, nr_unburned,
        prop_burned, any_burned),
    by = c("fire_event", "site", "pop", "qu")) %>%
  left_join(
    df_fire_ambiguous %>%
      rename(ambiguous_event = fire_event),
    by = c(
      "site", "pop",
      "year" = "candidate_year")) %>%
  mutate(
    timing_ambiguous = !is.na(ambiguous_event),
    fire_event_pop = !is.na(fire_event),
    
    fire_prop_quad = case_when(
      timing_ambiguous ~ NA_real_,
      !fire_event_pop ~ 0,
      !is.na(prop_burned) ~ prop_burned,
      TRUE ~ NA_real_),
    
    dist_transition = case_when(
      timing_ambiguous ~ NA_real_,
      !fire_event_pop ~ 0,
      !is.na(any_burned) ~ any_burned,
      TRUE ~ NA_real_))


# Previous annual fire exposure ------------------------------------------------
df_fire_previous <- df_fire_transition_quad %>%
  transmute(
    site, pop, qu,
    year = year + 1,
    disturbance_prev = dist_transition,
    fire_prop_prev = fire_prop_quad)

df_transition <- df_transition %>%
  left_join(
    df_fire_transition_quad %>%
      select(
        site, pop, qu, year,
        fire_event, ambiguous_event, timing_ambiguous,
        dist_transition, fire_prop_quad),
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_fire_previous,
    by = c("site", "pop", "qu", "year"))


# Individual fire exposure -----------------------------------------------------
df_fire_individual <- df_fire %>%
  filter(timing_known) %>%
  select(
    fire_event, id,
    burned_ind = burned,
    fire_severity_ind = fire_severity)

df_transition <- df_transition %>%
  left_join(
    df_fire_individual,
    by = c("fire_event", "id")) %>%
  mutate(
    burned_ind = case_when(
      timing_ambiguous ~ NA_real_,
      is.na(fire_event) ~ 0,
      !is.na(burned_ind) ~ burned_ind,
      TRUE ~ NA_real_))


# Previous individual fire exposure -------------------------------------------
df_fire_ind_previous <- df_transition %>%
  transmute(
    id,
    year = year + 1,
    burned_ind_prev = burned_ind)

df_transition <- df_transition %>%
  left_join(
    df_fire_ind_previous,
    by = c("id", "year"))

df_fire_transition_quad %>%
  count(timing_ambiguous, dist_transition)

df_transition %>%
  filter(!is.na(fire_event)) %>%
  count(fire_event, burned_ind)

df_fire_transition_quad %>%
  filter(!timing_ambiguous, fire_event_pop) %>%
  count(fire_event, is.na(dist_transition))

df_transition %>%
  filter(!is.na(fire_event)) %>%
  count(fire_event, burned_ind, useNA = "ifany")


# Check that fire joins did not duplicate demographic records ------------------
df_transition %>%
  count(id, year) %>%
  filter(n > 1)

nrow(df_transition)

df_transition %>%
  distinct(id, year) %>%
  nrow()


# Check that fire joins did not duplicate demographic records ------------------
df_transition %>%
  count(id, year) %>%
  filter(n > 1)

nrow(df_transition)

df_transition %>%
  distinct(id, year) %>%
  nrow()


# Recruitment correction ------------------------------------------------------
seedling_survival <- df_transition %>%
  filter(stage_clean == 1, !is.na(survives)) %>%
  summarise(seedling_survival = mean(survives)) %>%
  pull(seedling_survival)

df_recruit_quad <- df_recruit_quad %>%
  mutate(
    recruits_pess = if_else(
      recruitment_complete,
      as.numeric(observed_seedling + presumed_2yr), NA_real_),
    recruits_opt = if_else(
      recruitment_complete,
      observed_seedling + presumed_2yr / seedling_survival,
      NA_real_))

df_repro_quad %>%
  summarise(
    n_quad_years = n(),
    n_zero_active = sum(nr_active == 0),
    n_flowering_unknown = sum(is.na(nr_flowering)),
    n_scapes_unknown = sum(is.na(nr_scapes)))

df_recruit_quad %>%
  summarise(
    n_complete = sum(recruitment_complete),
    mean_pess = mean(recruits_pess, na.rm = TRUE),
    mean_opt = mean(recruits_opt, na.rm = TRUE),
    max_pess = max(recruits_pess, na.rm = TRUE),
    max_opt = max(recruits_opt, na.rm = TRUE))


# Recruitment working rows ----------------------------------------------------
df_recruit <- df_recruit_quad %>%
  left_join(
    df_fire_transition_quad %>%
      select(
        site, pop, qu, year, fire_event, ambiguous_event,
        timing_ambiguous, dist_transition, fire_prop_quad),
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_fire_previous,
    by = c("site", "pop", "qu", "year")) %>%
  transmute(
    row_type = "recruitment",
    site, pop, qu, year,
    monitored_t1, monitored_t2,
    recruitment_complete,
    observed_seedling, presumed_2yr, nr_recruit,
    recruits_pess, recruits_opt,
    nr_active, nr_flowering, nr_scapes, nr_invol,
    nr_flowering_other, mean_flowering_other,
    nr_scapes_other, mean_scapes_other,
    nr_quads, nr_active_pop,
    nr_flowering_pop, mean_flowering_pop,
    nr_scapes_pop, mean_scapes_pop,
    dist_transition, disturbance_prev,
    fire_prop_quad, fire_prop_prev,
    timing_ambiguous, fire_event, ambiguous_event)


# Individual working rows -----------------------------------------------------
df_ind <- df_transition %>%
  transmute(
    row_type = "individual",
    site, pop, qu, plant, id, year,
    state = state_clean,
    state_t1,
    stage = stage_clean,
    stage_t1,
    survives,
    enter_dormancy,
    reactivate,
    size_t0,
    size_t1,
    logsize_t0,
    logsize_t1,
    logsize_t0_2,
    logsize_t0_3,
    flower,
    fl_nr,
    invol,
    herb,
    dist_transition,
    disturbance_prev,
    fire_prop_quad,
    fire_prop_prev,
    burned_ind,
    burned_ind_prev,
    fire_severity_ind,
    timing_ambiguous,
    fire_event,
    ambiguous_event)


# Working data -----------------------------------------------------------------
# The final working data contain both individual demographic observations and
# quadrat-level recruitment observations. Recruitment is kept as separate rows
# rather than repeated for every individual within a quadrat.

df <- bind_rows(
  df_ind,
  df_recruit) %>%
  mutate(
    row_type = factor(row_type),
    year = as.numeric(year),
    site = factor(site),
    pop = factor(pop),
    state = factor(state),
    state_t1 = factor(state_t1),
    stage = factor(stage),
    stage_t1 = factor(stage_t1))


# Final working-data checks ----------------------------------------------------
df %>%
  count(row_type)

stopifnot(
  df %>%
    filter(row_type == "individual") %>%
    nrow() ==
    df %>%
    filter(row_type == "individual") %>%
    distinct(id, year) %>%
    nrow(),
  
  df %>%
    filter(row_type == "recruitment") %>%
    nrow() ==
    df %>%
    filter(row_type == "recruitment") %>%
    distinct(site, pop, qu, year) %>%
    nrow())

# Check individual vital-rate structure
df %>%
  filter(row_type == "individual") %>%
  count(state, survives)

df %>%
  filter(
    row_type == "individual",
    !is.na(size_t1)) %>%
  count(state, state_t1)

df %>%
  filter(
    row_type == "individual",
    !is.na(enter_dormancy)) %>%
  count(state, state_t1, enter_dormancy)

df %>%
  filter(
    row_type == "individual",
    !is.na(reactivate)) %>%
  count(state, state_t1, reactivate)

# Check recruitment structure
df %>%
  filter(row_type == "recruitment") %>%
  summarise(
    n = n(),
    n_complete = sum(recruitment_complete),
    n_pess = sum(!is.na(recruits_pess)),
    n_opt = sum(!is.na(recruits_opt)))

# Check ambiguous fire transitions
df %>%
  filter(timing_ambiguous) %>%
  count(row_type, dist_transition, fire_prop_quad)


# Save data --------------------------------------------------------------------
# write.csv(
#   df, row.names = FALSE,
#   file.path(
#     dir_data,
#     paste0("ab_", v_sp_abb, "_df_workdata_260820.csv")))