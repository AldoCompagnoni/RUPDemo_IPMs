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
      annual_transition & lead(state_clean) == "active",
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


# Monitoring end ---------------------------------------------------------------
df_monitor_end <- df_annual %>%
  filter(state_clean != "discontinued") %>%
  group_by(site, pop, qu) %>%
  summarise(last_monitor_year = max(year), .groups = "drop")


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

df_recruit_quad <- df_demog_june %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_recruit_cohort %>%
      count(
        site, pop, qu, recruit_year, recruit_type,
        name = "nr_recruits") %>%
      pivot_wider(
        names_from = recruit_type,
        values_from = nr_recruits,
        values_fill = 0),
    by = c("site", "pop", "qu", "year" = "recruit_year")) %>%
  mutate(
    seedling = replace_na(seedling, 0),
    new_adult = replace_na(new_adult, 0))

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


# Recruitment by fecundity year and quadrat -----------------------------------
df_recruit_quad <- df_demog_june %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_recruit_entry %>%
      filter(fecundity_usable) %>%
      count(
        site, pop, qu, fecundity_year, recruit_class,
        name = "nr_recruits") %>%
      pivot_wider(
        names_from = recruit_class,
        values_from = nr_recruits,
        values_fill = 0),
    by = c(
      "site", "pop", "qu",
      "year" = "fecundity_year")) %>%
  mutate(
    observed_seedling = replace_na(observed_seedling, 0),
    presumed_2yr = replace_na(presumed_2yr, 0))


# Parental reproductive abundance ---------------------------------------------
df_repro_quad <- df_transition %>%
  filter(state_clean == "active") %>%
  group_by(site, pop, qu, year) %>%
  summarise(
    nr_active = n(),
    nr_flowering = sum(flower == 1, na.rm = TRUE),
    nr_flower_na = sum(is.na(flower)),
    nr_scapes = if_else(
      all(is.na(fl_nr)), NA_real_, sum(fl_nr, na.rm = TRUE)),
    nr_invol_obs = sum(!is.na(invol)),
    nr_invol = if_else(
      nr_invol_obs == 0, NA_real_, sum(invol, na.rm = TRUE)),
    .groups = "drop")

df_recruit_quad <- df_recruit_quad %>%
  left_join(
    df_repro_quad,
    by = c("site", "pop", "qu", "year"))


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


# Population-level reproductive abundance -------------------------------------
df_repro_pop <- df_repro_quad %>%
  group_by(site, pop, year) %>%
  summarise(
    nr_quads = n(),
    nr_flowering_pop = sum(nr_flowering),
    nr_scapes_pop = sum(nr_scapes, na.rm = TRUE),
    .groups = "drop")

df_recruit_quad <- df_recruit_quad %>%
  left_join(df_repro_pop, by = c("site", "pop", "year")) %>%
  mutate(
    nr_recruit = observed_seedling + presumed_2yr,
    nr_flowering_other = nr_flowering_pop - nr_flowering,
    nr_scapes_other = nr_scapes_pop - replace_na(nr_scapes, 0))

df_recruit_quad %>%
  filter(nr_recruit > 0, nr_flowering == 0) %>%
  count(other_flowering = nr_flowering_other > 0)

df_recruit_quad %>%
  filter(nr_recruit > 0, nr_flowering_pop == 0) %>%
  select(
    site, pop, qu, year, observed_seedling, presumed_2yr,
    nr_flowering, nr_flowering_pop)

df_repro_pop %>%
  arrange(site, pop, year) %>%
  group_by(site, pop) %>%
  summarise(
    min_quads = min(nr_quads),
    max_quads = max(nr_quads),
    .groups = "drop")
