# Data exploration (only the site with the longest time series)
# McGraw 2017 - Panax quinquefolius;

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.09.22

# Study organism: Panax quinquefolius L.
# Link: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.9.4
# Meta data link: https://portal.edirepository.org/nis/metadataviewer?packageid=edi.9.4
# Citing publication: SMcGraw et al. 2017. Long Term Research in Environmental Biology: Demographic census data for thirty natural populations of American Ginseng: 1998-2016 ver 4. Environmental Data Initiative.
# Time period: 1998-2016


# Setting the stage ------------------------------------------------------------
# rm(list = ls())
set.seed(100)
options(stringsAsFactors = F)


# Packages --------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(
  MASS,
  tidyverse,
  patchwork,
  skimr,
  ipmr,
  binom,
  bbmle,
  janitor,
  lme4,
  GGally,
  brms,
  loo)


# Specification ---------------------------------------------------------------
v_brm_suffix <- ""
v_head <- c('mcgraw')
v_species <- c('Panax quinquefolius')
custom_delimiter <- c()
v_years_re <- c()

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_ggp_suffix <- paste(tools::toTitleCase(v_head), '-', v_species)

# Keep manual model choices empty for AICc selection.
# For survival, growth, dormancy, flowering, and scape number, values 0:3
# restrict selection to the corresponding polynomial degree.
v_mod_set_su   <- c()
v_mod_set_gr   <- c()
v_mod_set_do   <- c()
v_mod_set_fl   <- c()
v_mod_set_fl_n <- c()

# Recruitment candidates are indexed 0:8 below.
v_mod_set_re <- c()


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


# Functions -------------------------------------------------------------------
source('helper_functions/plot_binned_prop.R')
source('helper_functions/plot_binned_prop_year.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data ------------------------------------------------------------------------
disturbance_case_levels <- c(
  'Undisturbed',
  'Current-year fire',
  'Previous-year fire',
  'Fire in both years')

disturbance_case_cols <- c(
  'Undisturbed' = 'black',
  'Current-year fire' = 'red',
  'Previous-year fire' = 'purple',
  'Fire in both years' = 'magenta')

df_og <- read.csv(
  file.path(
    dir_data,
    paste0('df_orig.csv')))
str(df_og)

df_meta <- read.csv(
  file.path(
    dir_data,
    paste0('df_meta.csv')))
str(df_meta)

# Preserve original data -------------------------------------------------------
df_og_raw <- df_og

# Repair encoding enough for R to handle the strings ---------------------------
# fixes the following examplaries to displayable UTF-8 text:
# fr 1060 40Î_ 11 cm
# 257Î_ 0.5m frm 262
# fr 169 20Î_ 25 cm
loc_chr <- as.character(df_og$loc)

bad <- is.na(iconv(loc_chr, from = "", to = "UTF-8", sub = NA)) &
  !is.na(loc_chr)

loc_chr[bad] <- iconv(
  loc_chr[bad], from = "WINDOWS-1252", to = "UTF-8", sub = "?")

df_og$loc <- factor(loc_chr)



# Metadata --------------------------------------------------------------------

meta_ref <- tribble(
  ~variable, ~definition, ~unit,
  "population", "Unique identification number for each population.",
  "dimensionless",
  "year", "Year in which the census data were collected.", "year",
  "id", "Individual marking-nail ID; unique within population.",
  "dimensionless",
  "cluster", "Plant cluster recorded in the phototrail.", "dimensionless",
  "age", "Plant age; age 0 denotes the year of germination.", "dimensionless",
  "obs_leaf_num", "Observed leaf number at the spring census.",
  "dimensionless",
  "inf_leaf_num", "Inferred leaf number when leaves could not be observed.",
  "dimensionless",
  "lflt_arr", "Leaflet arrangement among leaves, e.g. 553.", "text",
  "lflt_tot", "Total number of leaflets.", "dimensionless",
  "stalk_ht", "Stalk height from ground to leaf branching point.", "cm",
  "lll1", "Length of longest leaflet on leaf 1.", "cm",
  "wll1", "Width of longest leaflet on leaf 1.", "cm",
  "lll2", "Length of longest leaflet on leaf 2.", "cm",
  "wll2", "Width of longest leaflet on leaf 2.", "cm",
  "lll3", "Length of longest leaflet on leaf 3.", "cm",
  "wll3", "Width of longest leaflet on leaf 3.", "cm",
  "lll4", "Length of longest leaflet on leaf 4.", "cm",
  "wll4", "Width of longest leaflet on leaf 4.", "cm",
  "obs_la", "Observed total leaf area calculated from leaflet dimensions.",
  "cm2",
  "inf_la", "Observed or inferred leaf area when the plant was not measurable.",
  "cm2",
  "f_buds", "Reproductive status based on presence of buds or flowers.",
  "nominal",
  "seeds", "Estimated actual seed count.", "dimensionless",
  "red_frts", "Number of red or partially ripened fruits.", "dimensionless",
  "grn_frts", "Number of green fruits.", "dimensionless",
  "tot_frts", "Total number of red and green fruits.", "dimensionless",
  "loc", "Location description used to relocate plants within clusters.",
  "text",
  "harvest", "Evidence that a plant had been harvested.", "nominal",
  "harv_time", "Timing of observed harvest.", "nominal",
  "browse", "Evidence of browsing based on damaged stalks or petioles.",
  "nominal",
  "browse_prcnt", "Estimated percentage of leaves removed by browsing.",
  "percent",
  "browse_time", "Timing of browsing relative to censuses.", "nominal",
  "persistence", "Persistence/status class used to describe plant fate.",
  "nominal",
  "insect", "Presence of insect damage.", "nominal",
  "thrips", "Presence of damage attributed to thrips.", "nominal",
  "fungal", "Presence of fungal disease.", "nominal",
  "fungal_cat", "Category of fungal disease.", "nominal")


# Properties of downloaded data -----------------------------------------------

meta_data <- tibble(variable = names(df_og)) %>%
  mutate(
    class = map_chr(
      variable, ~ paste(class(df_og[[.x]]), collapse = ", ")),
    n = nrow(df_og),
    n_missing = map_int(variable, ~ sum(is.na(df_og[[.x]]))),
    pct_missing = round(100 * n_missing / n, 2),
    n_unique = map_int(
      variable, ~ n_distinct(df_og[[.x]], na.rm = TRUE)),
    minimum = map_chr(variable, function(x) {
      z <- df_og[[x]]
      if (!is.numeric(z) || all(is.na(z))) return(NA_character_)
      as.character(min(z, na.rm = TRUE))
    }),
    maximum = map_chr(variable, function(x) {
      z <- df_og[[x]]
      if (!is.numeric(z) || all(is.na(z))) return(NA_character_)
      as.character(max(z, na.rm = TRUE))
    }),
    encoding_issues = map_int(variable, function(x) {
      z <- df_og[[x]]
      if (!is.character(z) && !is.factor(z)) return(0L)
      
      z <- as.character(z)
      sum(
        is.na(iconv(z, from = "", to = "UTF-8", sub = NA)) &
          !is.na(z))
    })) %>%
  left_join(meta_ref, by = "variable")

meta_data



# Sampling structure over time ------------------------------------------------

df_sampling <- df_og %>%
  left_join(
    df_meta %>%
      select(population, state),
    by = "population") %>%
  mutate(
    population_id = paste0(state, " - Pop ", population),
    cluster_id = paste0("Pop ", population, " - Cl ", cluster))


# State -----------------------------------------------------------------------

p_state <- df_sampling %>%
  distinct(year, state) %>%
  ggplot(aes(x = year, y = state)) +
  geom_tile() +
  scale_x_continuous(breaks = 1998:2016) +
  labs(x = NULL, y = "State", title = "State") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


# Population ------------------------------------------------------------------

p_population <- df_sampling %>%
  distinct(year, population_id, population) %>%
  mutate(
    population_id = forcats::fct_reorder(
      population_id, population, .fun = min)) %>%
  ggplot(aes(x = year, y = population_id)) +
  geom_tile() +
  scale_x_continuous(breaks = 1998:2016) +
  labs(x = NULL, y = "Population", title = "Population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


# Cluster ---------------------------------------------------------------------

p_cluster <- df_sampling %>%
  filter(!is.na(cluster)) %>%
  distinct(year, population, cluster, cluster_id) %>%
  ggplot(aes(x = year, y = factor(cluster))) +
  geom_tile() +
  facet_wrap(~ population, scales = "free_y") +
  scale_x_continuous(breaks = seq(1998, 2016, 3)) +
  labs(
    x = "Year",
    y = "Cluster",
    title = "Clusters within populations") +
  theme_bw()


# Combined overview -----------------------------------------------------------

p_state / p_population / p_cluster +
  plot_layout(heights = c(1, 3, 7))


# Variables and where they go in the IPM --------------------------------------

df_ipm_variables <- meta_data %>%
  mutate(
    ipm_role = case_when(
      variable %in% c("population", "year", "id", "cluster") ~
        "Structure",
      variable %in% c(
        "obs_la", "inf_la", "obs_leaf_num", "inf_leaf_num",
        "lflt_tot", "stalk_ht") ~
        "Size / growth",
      variable %in% c(
        "lll1", "wll1", "lll2", "wll2",
        "lll3", "wll3", "lll4", "wll4") ~
        "Raw size measurement",
      variable == "persistence" ~
        "Survival / dormancy / recruitment",
      variable == "f_buds" ~
        "Flowering",
      variable %in% c(
        "seeds", "red_frts", "grn_frts", "tot_frts") ~
        "Fecundity",
      variable %in% c(
        "browse", "browse_prcnt", "browse_time",
        "harvest", "harv_time", "insect", "thrips",
        "fungal", "fungal_cat", "age") ~
        "Potential covariate",
      variable %in% c("loc", "lflt_arr") ~
        "Not core IPM",
      TRUE ~ "Unclassified"))

df_ipm_variables


# Create working data for simple mean IPM -------------------------------------

df_work <- df_og %>%
  dplyr::left_join(df_meta, by = "population") %>%
  dplyr::arrange(population, id, year) %>%
  dplyr::mutate(
    size = inf_la,
    recruit = as.integer(persistence == "NS"),
    alive = dplyr::case_when(
      persistence == "DEAD" ~ 0,
      !is.na(persistence) ~ 1,
      TRUE ~ NA_real_)) %>%
  dplyr::group_by(population, id) %>%
  dplyr::mutate(
    year_t1 = dplyr::lead(year),
    size_t1 = dplyr::lead(size),
    persistence_t1 = dplyr::lead(persistence),
    survives = dplyr::lead(alive)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    consecutive = year_t1 == year + 1,
    survives = dplyr::if_else(consecutive, survives, NA_real_),
    size_t1 = dplyr::if_else(consecutive, size_t1, NA_real_),
    reliable_survival = consecutive & year_t1 <= 2014,
    logsize_t0 = log(size),
    logsize_t1 = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  dplyr::rename(
    size_t0 = size,
    persistence_t0 = persistence) %>%
  dplyr::select(
    state, population,
    cluster, id, year, year_t1,
    persistence_t0, persistence_t1,
    size_t0, size_t1,
    logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
    survives, recruit, consecutive, reliable_survival)



# Save

write.csv(
  df_work,
  file.path(dir_data, "mcgraw_paqu_df_workdata.csv"),
  row.names = FALSE)
