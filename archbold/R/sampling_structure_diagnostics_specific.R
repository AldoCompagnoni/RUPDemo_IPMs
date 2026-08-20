# Sampling structure for the four demographic datasets -------------------------

library(tidyverse)

# Species ----------------------------------------------------------------------
v_spec01 <- c("Eriogonum longifolium")
v_spec02 <- c("Hypericum cumulicola")
v_spec03 <- c("Liatris ohlingerae")
v_spec04 <- c("Solidago odora")

v_sp_abb01  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_spec01, ' ')), 1, 2), collapse = '')))
v_sp_abb02  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_spec02, ' ')), 1, 2), collapse = '')))
v_sp_abb03  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_spec03, ' ')), 1, 2), collapse = '')))
v_sp_abb04  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_spec04, ' ')), 1, 2), collapse = '')))


# Directories ------------------------------------------------------------------
dir_data <- file.path("C:/code/RUPDemo_IPMs/archbold/data")


# Data -------------------------------------------------------------------------

df_erlo <- read_csv(
  list.files(
    file.path(dir_data, v_sp_abb01),
    pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)[1],
  col_types = cols(comment = col_character()))

df_hycu <- read_csv(
  list.files(
    file.path(dir_data, v_sp_abb02),
    pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)[1],
  col_types = cols(comment = col_character()))

df_lioh <- read_csv(
  list.files(
    file.path(dir_data, v_sp_abb03),
    pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)[1],
  col_types = cols(comments     = col_character(),
                   QAQCcomments = col_character()))

df_sood <- read_csv(
  list.files(
    file.path(dir_data, v_sp_abb04),
    pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)[1])




# Spatial and temporal inspection -----------------------------------------------

# Eriogonum longifolium ---------------------------------------------------------
# Spatial hierarchy: site > population > quadrat

df_erlo_space <- df_erlo %>%
  summarise(
    sites = n_distinct(site),
    populations = n_distinct(site, pop),
    quadrats = n_distinct(site, pop, qu))

df_erlo_space_site <- df_erlo %>%
  group_by(site) %>%
  summarise(
    populations = n_distinct(pop),
    quadrats = n_distinct(pop, qu), .groups = "drop")

df_erlo_sampling <- df_erlo %>%
  transmute(
    year = as.integer(str_sub(date, 1, 4)),
    site, pop, qu) %>%
  filter(!is.na(year), !is.na(site), !is.na(pop), !is.na(qu)) %>%
  distinct() %>%
  arrange(site, pop, qu, year) %>%
  mutate(
    unit = paste0("S", site, " P", pop, " Q", qu),
    unit = fct_rev(fct_inorder(unit)))

df_erlo_sampling_unit <- df_erlo_sampling %>%
  group_by(site, pop, qu, unit) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year), .groups = "drop")

df_erlo_sampling_year <- df_erlo_sampling %>%
  count(year, name = "n_quadrats")

p_erlo_sampling <- ggplot(
  df_erlo_sampling, aes(x = year, y = unit)) +
  geom_point(shape = 15, size = 1.5) +
  scale_x_continuous(
    breaks = seq(
      min(df_erlo_sampling$year),
      max(df_erlo_sampling$year), by = 1)) +
  labs(
    title = v_spec01,
    subtitle = "site > population > quadrat",
    x = "Year", y = "Sampling unit") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank())


# Hypericum cumulicola ----------------------------------------------------------
# Spatial hierarchy: site > gap

df_hycu_space <- df_hycu %>%
  summarise(
    sites = n_distinct(Site),
    gaps = n_distinct(Site, Gap))

df_hycu_space_site <- df_hycu %>%
  group_by(Site) %>%
  summarise(gaps = n_distinct(Gap), .groups = "drop")

df_hycu_sampling <- df_hycu %>%
  transmute(year = as.integer(Year), site = Site, gap = Gap) %>%
  filter(!is.na(year), !is.na(site), !is.na(gap)) %>%
  distinct() %>%
  arrange(site, gap, year) %>%
  mutate(
    unit = paste0("S", site, " G", gap),
    unit = fct_rev(fct_inorder(unit)))

df_hycu_sampling_unit <- df_hycu_sampling %>%
  group_by(site, gap, unit) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year), .groups = "drop")

df_hycu_sampling_year <- df_hycu_sampling %>%
  count(year, name = "n_gaps")

p_hycu_sampling <- ggplot(
  df_hycu_sampling, aes(x = year, y = unit)) +
  geom_point(shape = 15, size = 1.5) +
  scale_x_continuous(
    breaks = seq(
      min(df_hycu_sampling$year),
      max(df_hycu_sampling$year), by = 1)) +
  labs(
    title = v_spec02,
    subtitle = "site > gap",
    x = "Year", y = "Sampling unit") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank())


# Liatris ohlingerae ------------------------------------------------------------
# Spatial hierarchy: site > population

df_lioh_space <- df_lioh %>%
  summarise(
    sites = n_distinct(site),
    populations = n_distinct(site, population))

df_lioh_space_site <- df_lioh %>%
  group_by(site) %>%
  summarise(
    populations = n_distinct(population), .groups = "drop")

df_lioh_sampling <- df_lioh %>%
  transmute(
    year = as.integer(year),
    site, population) %>%
  filter(!is.na(year), !is.na(site), !is.na(population)) %>%
  distinct() %>%
  arrange(site, population, year) %>%
  mutate(
    unit = paste0("S", site, " P", population),
    unit = fct_rev(fct_inorder(unit)))

df_lioh_sampling_unit <- df_lioh_sampling %>%
  group_by(site, population, unit) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year), .groups = "drop")

df_lioh_sampling_year <- df_lioh_sampling %>%
  count(year, name = "n_populations")

p_lioh_sampling <- ggplot(
  df_lioh_sampling, aes(x = year, y = unit)) +
  geom_point(shape = 15, size = 2) +
  scale_x_continuous(
    breaks = seq(
      min(df_lioh_sampling$year),
      max(df_lioh_sampling$year), by = 1)) +
  labs(
    title = v_spec03,
    subtitle = "site > population",
    x = "Year", y = "Sampling unit") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())


# Solidago odora ---------------------------------------------------------------
# Spatial hierarchy: quadrat

df_sood_space <- df_sood %>%
  summarise(quadrats = n_distinct(quad))

df_sood_sampling <- df_sood %>%
  transmute(year = as.integer(year), quad) %>%
  filter(!is.na(year), !is.na(quad)) %>%
  distinct() %>%
  arrange(quad, year) %>%
  mutate(
    unit = paste0("Q", quad),
    unit = fct_rev(fct_inorder(unit)))

df_sood_sampling_unit <- df_sood_sampling %>%
  group_by(quad, unit) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year), .groups = "drop")

df_sood_sampling_year <- df_sood_sampling %>%
  count(year, name = "n_quadrats")

p_sood_sampling <- ggplot(
  df_sood_sampling, aes(x = year, y = unit)) +
  geom_point(shape = 15, size = 2) +
  scale_x_continuous(
    breaks = seq(
      min(df_sood_sampling$year),
      max(df_sood_sampling$year), by = 1)) +
  labs(
    title = v_spec04,
    subtitle = "quadrat",
    x = "Year", y = "Sampling unit") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())


# Inspect ----------------------------------------------------------------------

df_erlo_space
df_erlo_space_site
df_erlo_sampling_unit
df_erlo_sampling_year
p_erlo_sampling

df_hycu_space
df_hycu_space_site
df_hycu_sampling_unit
df_hycu_sampling_year
p_hycu_sampling

df_lioh_space
df_lioh_space_site
df_lioh_sampling_unit
df_lioh_sampling_year
p_lioh_sampling

df_sood_space
df_sood_sampling_unit
df_sood_sampling_year
p_sood_sampling