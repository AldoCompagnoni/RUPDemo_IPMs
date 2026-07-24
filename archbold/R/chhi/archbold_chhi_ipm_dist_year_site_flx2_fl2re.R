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
v_mod_set_gr <- c()
v_mod_set_fl <- c()
v_mod_set_fr <- c()
v_mod_set_re <- c()


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


# Recruitment classification --------------------------------------------------
# Recruits include:
# 1. plants classified as seedlings;
# 2. plants first observed after site monitoring began and recorded as
#    "new plant" in at least one quarterly census.
#
# New plants from a site's first census are excluded because they may have
# existed before monitoring started.
# Recruitment was initially identified solely from the annual stage classification, 
# with individuals recorded as seedlings (`astg = 1`) classified as recruits. 
# Diagnostic checks showed that no individual was classified 
# as a seedling recruit in more than one year and that 
# all seedling records occurred in the individual's first observed year. 
# However, examination of the quarterly census status variables revealed 
# additional first-observed individuals coded as “new plant” but already 
# classified as vegetative or bolting at the annual census. Among the records 
# that could be assigned relative to the start of site monitoring, 828 occurred 
# during a site's initial census and were not considered recruits because their 
# previous presence was unknown. A further 280 new plants appeared after 
# monitoring had begun, comprising 264 vegetative and 16 bolting individuals. 
# These individuals were classified as recruits because they were first observed 
# during an established monitoring period and were explicitly recorded as new 
# plants. Recruitment was therefore defined as either seedling stage or first 
# observation after site establishment with at least one quarterly census status 
# indicating a new plant.


df_recruit_flag <- df_og %>%
  rename(
    plant_id = identifier,
    year = year0) %>%
  filter(
    !is.na(site),
    !is.na(plant_id),
    !is.na(year)) %>%
  mutate(
    new_plant_status = if_any(
      c(s_03, s_06, s_09, s_12), ~ .x == 3)) %>%
  group_by(site) %>%
  mutate(site_first_year = min(year, na.rm = TRUE)) %>%
  group_by(site, plant_id) %>%
  mutate(first_plant_year = min(year, na.rm = TRUE)) %>%
  group_by(site, plant_id, year) %>%
  summarise(
    seedling = any(astg == 1, na.rm = TRUE),
    new_plant_status = any(new_plant_status, na.rm = TRUE),
    site_first_year = first(site_first_year),
    first_plant_year = first(first_plant_year),
    .groups = "drop") %>%
  mutate(
    recruit = as.integer(
      seedling |
        (new_plant_status &
           year == first_plant_year &
           year > site_first_year)),
    recruit_source = case_when(
      seedling ~ "Seedling stage",
      recruit == 1 ~ "New plant after site establishment",
      TRUE ~ "Not recruit"))


# Mean data frame --------------------------------------------------------------
df <- df_og %>%
  rename(
    plant_id = identifier,
    year = year0,
    survival = survival_1) %>%
  arrange(site, plant_id, year, survival) %>%
  group_by(site, plant_id, year) %>%
  summarise(
    survives = if (all(is.na(survival))) {
      NA_real_
    } else {
      min(survival, na.rm = TRUE)
    },
    size_t0 = if (all(is.na(dia))) {
      NA_real_
    } else {
      max(dia, na.rm = TRUE)
    },
    size_t1 = if (all(is.na(dia_1))) {
      NA_real_
    } else {
      max(dia_1, na.rm = TRUE)
    },
    flower = if (all(is.na(hd))) {
      NA_real_
    } else {
      max(hd, na.rm = TRUE)
    },
    .groups = "drop") %>%
  left_join(
    df_recruit_flag %>%
      select(site, plant_id, year, recruit, recruit_source),
    by = c("site", "plant_id", "year")) %>%
  mutate(
    plant_id = factor(plant_id),
    recruit = replace_na(recruit, 0L),
    recruit_source = replace_na(recruit_source, "Not recruit"),
    fl_nr = flower,
    flower = if_else(flower > 0, 1, flower),
    logsize_t0 = log(size_t0),
    logsize_t1 = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  full_join(df_fire, by = c("site", "year")) %>%
  mutate(
    fire = case_when(
      is.na(fire) ~ "No fire",
      fire == "burn" ~ "Fire",
      TRUE ~ NA_character_),
    fire = factor(fire, levels = c("No fire", "Fire")))


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
      color = "black", size = 1) +
    # fire points
    geom_point(
      data = grow_pts %>% filter(fire == "Fire"),
      aes(logsize_t0, logsize_t1),
      color = "red", size = 1) +
    # no fire line
    geom_line(aes(logsize_t0, pred),
              data = pred_no_fire, color = "black", lwd = 1) +
    # fire line
    geom_line(aes(logsize_t0, pred), 
              data = pred_fire, color = "red", lwd = 1) +
    geom_abline(intercept = 0, slope = 1, color = "blue", lty = 2) +
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


# Flowering data ---------------------------------------------------------------
df_fl <- df %>%
  mutate(
    flower = if_else(!is.na(fl_nr) & fl_nr > 0, 1, 0),
    fire   = factor(fire, levels = c("No fire", "Fire")),
    year   = factor(year)) %>%
  filter(
    !is.na(logsize_t0),
    !is.na(logsize_t0_2),
    !is.na(logsize_t0_3),
    !is.na(fire),
    !is.na(year))


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
  filter(!is.na(fl_nr), fl_nr > 0) %>%
  filter(fl_nr == round(fl_nr)) %>%
  mutate(year = factor(year))
  

ctrl_fl_n <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

mod_fl_n_0 <- glmer.nb(
  fl_nr ~ fire + (1 | year),
  data = df_fl_cond, control = ctrl_fl_n)

mod_fl_n_1 <- glmer.nb(
  fl_nr ~ logsize_t0 + fire + (1 | year),
  data = df_fl_cond, control = ctrl_fl_n)

mod_fl_n_2 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + fire + (1 | year),
  data = df_fl_cond, control = ctrl_fl_n)

mod_fl_n_3 <- glmer.nb(
  fl_nr ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire + (1 | year),
  data = df_fl_cond, control = ctrl_fl_n)

mods_fl_n <- list(mod_fl_n_0, mod_fl_n_1, mod_fl_n_2, mod_fl_n_3)

mods_fl_n_dAIC <- bbmle::AICctab(mods_fl_n, weights = TRUE, sort = FALSE)$dAIC

mods_fl_n_sorted <- order(mods_fl_n_dAIC)

mod_fl_n_best <- mods_fl_n[[mods_fl_n_sorted[1]]]
v_mod_fl_n_index <- mods_fl_n_sorted[1] - 1


# Flowerhead production in t0 affects recruits in t1 --------------------------

df_ind <- df %>%
  filter(!is.na(site), !is.na(plant_id), !is.na(year)) %>%
  transmute(site, plant_id, year, fh_nr = replace_na(fl_nr, 0),
            recruit = replace_na(recruit, 0L))

df_site_year <- df_ind %>% 
  group_by(site, year) %>% 
  summarise(
    fh = sum(fh_nr[recruit == 0], na.rm = TRUE),
    re = sum(recruit == 1, na.rm = TRUE),
    .groups = "drop")

df_re <- df_site_year %>%
  transmute(site, year_t1 = year, re_t1 = re) %>%
  left_join(
    df_site_year %>%
      transmute(site, year_t0 = year, year_t1 = year + 1, fh_t0 = fh),
    by = c("site", "year_t1")) %>%
  mutate(
    site = factor(site),
    year_t0 = factor(year_t0),
    year_t1 = factor(year_t1),
    re_t1_log = log1p(re_t1),
    fh_t0_log = log1p(fh_t0)) %>% 
  filter(
    !is.na(re_t1_log),
    !is.na(fh_t0_log),
    !is.na(site),
    !is.na(year_t0),
    !is.na(year_t1))

mod_re_00 <- glm(
  re_t1_log ~ fh_t0_log,
  data = df_re, family = gaussian())

mod_re_10 <- lmer(
  re_t1_log ~ fh_t0_log + (1 | site),
  data = df_re, REML = FALSE, 
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_11 <- lmer(
  re_t1_log ~ fh_t0_log + (1 | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_12 <- lmer(
  re_t1_log ~ fh_t0_log + (1 | site) + (1 | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_20 <- lmer(
  re_t1_log ~ fh_t0_log + (fh_t0_log | site),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_21 <- lmer(
  re_t1_log ~ fh_t0_log + (fh_t0_log | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_22 <- lmer(
  re_t1_log ~ fh_t0_log + (fh_t0_log | site) + (fh_t0_log | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_30 <- lmer(
  re_t1_log ~ fh_t0_log + (0 + fh_t0_log | site),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_31 <- lmer(
  re_t1_log ~ fh_t0_log + (0 + fh_t0_log | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_re_32 <- lmer(
  re_t1_log ~ fh_t0_log + (0 + fh_t0_log | site) + (0 + fh_t0_log | year_t1),
  data = df_re, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


mods_re <- list(mod_re_00, 
                mod_re_10, mod_re_11, mod_re_12, 
                mod_re_20, mod_re_21, mod_re_22, 
                mod_re_30, mod_re_31, mod_re_32)
mods_re_dAIC <- bbmle::AICctab(mods_re, weights = TRUE, sort = FALSE)$dAIC
mods_re_sorted <- order(mods_re_dAIC)
if (length(v_mod_set_re) == 0) {
  mod_re_index_bestfit <- mods_re_sorted[1]
  v_mod_re_index       <- mod_re_index_bestfit - 1 
} else {
  mod_re_index_bestfit <- v_mod_set_re +1
  v_mod_re_index       <- v_mod_set_re
}

mod_re_best <- mods_re[[mod_re_index_bestfit]]
mod_re_ranef_site <- if (
  inherits(mod_re_best, "merMod") &&
  "site" %in% names(coef(mod_re_best))) {
  coef(mod_re_best)$site
} else {
  NULL
}

mod_re_ranef_year_t1 <- if (
  inherits(mod_re_best, "merMod") &&
  "year_t1" %in% names(coef(mod_re_best))) {
  coef(mod_re_best)$year_t1
} else {
  NULL
}

# Recruitment plots: flowerheads at t0 -> recruits at t1

predict_re_safe <- function(model, newdata, re.form = NULL) {
  
  if (inherits(model, "merMod")) {
    out <- predict(
      model,
      newdata = newdata,
      re.form = re.form,
      allow.new.levels = TRUE
    )
  } else {
    out <- predict(
      model,
      newdata = newdata
    )
  }
  
  return(out)
}


# site-specific point shapes
shape_values_re <- c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1)

site_shapes_re <- setNames(
  shape_values_re[seq_along(levels(df_re$site))],
  levels(df_re$site)
)


re_yr_plots <- function(i) {
  
  re_temp <- df_re %>%
    filter(year_t1 == i)
  
  x_temp_re <- seq(
    min(re_temp$fh_t0_log, na.rm = TRUE),
    max(re_temp$fh_t0_log, na.rm = TRUE),
    length.out = 100
  )
  
  # population-level fixed-effect trend
  pred_fixed_re <- data.frame(
    fh_t0_log = x_temp_re,
    site      = factor(levels(df_re$site)[1], levels = levels(df_re$site)),
    year_t1   = factor(i, levels = levels(df_re$year_t1))
  )
  
  pred_fixed_re$pred_re <- predict_re_safe(
    model   = mod_re_best,
    newdata = pred_fixed_re,
    re.form = NA
  )
  
  # full conditional trends:
  # includes site random effects, year_t1 random effects, or both if present
  pred_site_re <- tidyr::expand_grid(
    fh_t0_log = x_temp_re,
    site      = levels(df_re$site)
  ) %>%
    mutate(
      site    = factor(site, levels = levels(df_re$site)),
      year_t1 = factor(i, levels = levels(df_re$year_t1))
    )
  
  pred_site_re$pred_re <- predict_re_safe(
    model   = mod_re_best,
    newdata = pred_site_re,
    re.form = NULL
  )
  
  ggplot() +
    geom_point(
      data = re_temp,
      aes(
        x = fh_t0_log,
        y = re_t1_log,
        color = site,
        shape = site
      ),
      alpha = 0.75,
      size = 2
    ) +
    geom_line(
      data = pred_site_re,
      aes(
        x = fh_t0_log,
        y = pred_re,
        color = site,
        group = site
      ),
      linewidth = 0.7,
      alpha = 0.65
    ) +
    geom_line(
      data = pred_fixed_re,
      aes(
        x = fh_t0_log,
        y = pred_re
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.2
    ) +
    scale_color_viridis_d(option = "turbo", end = 0.95) +
    scale_shape_manual(values = site_shapes_re) +
    theme_bw() +
    labs(
      title = i,
      x = expression("log(1 + flowerheads)"[t0]),
      y = expression("log(1 + recruits)"[t1]),
      color = "Site",
      shape = "Site"
    ) +
    theme(
      text = element_text(size = 5),
      legend.position = "none"
    )
}


re_yrs <- lapply(
  sort(unique(df_re$year_t1)),
  re_yr_plots
)

g_re_years <- wrap_plots(re_yrs) +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = "Recruitment - year specific",
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  )

g_re_years

# Recruit size distribution -----------------------------------------------

df_recr <- df %>%
  filter(recruit == 1) %>%
  filter(size_t0 > 0) %>%
  mutate(logsize_t0 = log(size_t0))

rc_sz <- data.frame(
  mean = mean(df_recr$logsize_t0, na.rm = TRUE),
  sd   = sd(df_recr$logsize_t0, na.rm = TRUE))


# Exporting parameter estimates ------------------------------------------------

# Expected model objects:
# mod_su_best
# mod_gr_best
# mod_fl_best
# mod_fl_n_best
# mod_re_best


# Helper functions -------------------------------------------------------------

empty_coef_df <- function() {
  return(data.frame(
    coefficient = character(),
    value       = numeric()
  ))
}


is_mixed_model <- function(model) {
  return(inherits(model, "merMod"))
}


get_fixef_safe <- function(model) {
  
  if(is_mixed_model(model)) {
    out <- fixef(model)
  } else {
    out <- coef(model)
  }
  
  return(out)
}


get_group_coef_safe <- function(model, group_var) {
  
  if(!is_mixed_model(model)) {
    return(NULL)
  }
  
  coef_list <- coef(model)
  
  if(!(group_var %in% names(coef_list))) {
    return(NULL)
  }
  
  out <- coef_list[[group_var]]
  
  return(out)
}


find_model_term <- function(model_terms, aliases) {
  
  aliases <- unlist(aliases)
  
  for(alias in aliases) {
    
    if(grepl("^regex:", alias)) {
      pattern <- sub("^regex:", "", alias)
      hit <- grep(pattern, model_terms, value = TRUE)
    } else {
      hit <- model_terms[model_terms == alias]
    }
    
    if(length(hit) > 0) {
      return(hit[1])
    }
  }
  
  return(NA_character_)
}


bind_coef_rows <- function(x) {
  
  x <- x[!vapply(x, is.null, logical(1))]
  
  if(length(x) == 0) {
    return(empty_coef_df())
  }
  
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  
  return(out)
}


extract_fixed_pars <- function(model, term_map, fill_missing = TRUE) {
  
  model_coefs <- get_fixef_safe(model)
  model_terms <- names(model_coefs)
  
  out <- lapply(names(term_map), function(coef_name) {
    
    model_term <- find_model_term(
      model_terms = model_terms,
      aliases     = term_map[[coef_name]]
    )
    
    if(is.na(model_term)) {
      
      if(fill_missing) {
        value <- 0
      } else {
        return(NULL)
      }
      
    } else {
      value <- unname(model_coefs[model_term])
    }
    
    data.frame(
      coefficient = coef_name,
      value       = value
    )
  })
  
  out <- bind_coef_rows(out)
  
  return(out)
}


extract_group_pars <- function(model, group_var, term_map, fill_missing = FALSE) {
  
  coef_matrix <- get_group_coef_safe(
    model     = model,
    group_var = group_var
  )
  
  if(is.null(coef_matrix)) {
    return(empty_coef_df())
  }
  
  model_terms <- colnames(coef_matrix)
  
  out <- lapply(names(term_map), function(coef_prefix) {
    
    model_term <- find_model_term(
      model_terms = model_terms,
      aliases     = term_map[[coef_prefix]]
    )
    
    if(is.na(model_term)) {
      
      if(fill_missing) {
        value <- rep(0, nrow(coef_matrix))
      } else {
        return(NULL)
      }
      
    } else {
      value <- coef_matrix[, model_term]
    }
    
    data.frame(
      coefficient = paste0(coef_prefix, "_", rownames(coef_matrix)),
      value       = value
    )
  })
  
  out <- bind_coef_rows(out)
  
  return(out)
}


check_duplicate_coefficients <- function(x, object_name = "parameter table") {
  
  dup <- x$coefficient[duplicated(x$coefficient)]
  
  if(length(dup) > 0) {
    stop(
      object_name,
      " contains duplicated coefficient names: ",
      paste(unique(dup), collapse = ", ")
    )
  }
}


coef_df_to_list <- function(x) {
  
  if(nrow(x) == 0) {
    return(list())
  }
  
  check_duplicate_coefficients(x)
  
  out <- x %>%
    mutate(coefficient = as.character(coefficient)) %>%
    pivot_wider(
      names_from  = coefficient,
      values_from = value
    ) %>%
    as.list()
  
  return(out)
}


extract_re_ranef_devs <- function(model, group_var, group_label) {
  
  if (!inherits(model, "merMod")) {
    return(empty_coef_df())
  }
  
  ranef_list <- ranef(model)
  
  if (!(group_var %in% names(ranef_list))) {
    return(empty_coef_df())
  }
  
  ranef_matrix <- ranef_list[[group_var]]
  
  out <- list()
  
  if ("(Intercept)" %in% colnames(ranef_matrix)) {
    
    out[["intercept"]] <- data.frame(
      coefficient = paste0(
        "re_u0_", group_label, "_",
        rownames(ranef_matrix)
      ),
      value = ranef_matrix[, "(Intercept)"]
    )
  }
  
  if ("fh_t0_log" %in% colnames(ranef_matrix)) {
    
    out[["slope"]] <- data.frame(
      coefficient = paste0(
        "re_u1_", group_label, "_",
        rownames(ranef_matrix)
      ),
      value = ranef_matrix[, "fh_t0_log"]
    )
  }
  
  out <- bind_coef_rows(out)
  
  return(out)
}

# Expected term maps -----------------------------------------------------------

size_term_map <- list(
  b0 = c("(Intercept)"),
  b1 = c("logsize_t0"),
  b2 = c("logsize_t0_2", "I(logsize_t0^2)"),
  b3 = c("logsize_t0_3", "I(logsize_t0^3)"),
  bf = c("fireFire", "regex:^fire"))


make_size_term_map <- function(prefix) {
  
  out <- size_term_map
  names(out) <- paste0(prefix, names(out))
  
  return(out)
}


su_term_map  <- make_size_term_map("surv_")
gr_term_map  <- make_size_term_map("grow_")
fl_term_map  <- make_size_term_map("fl_")
fln_term_map <- make_size_term_map("fln_")
re_term_map  <- list(re_b0 = c("(Intercept)"), re_b1 = c("fh_t0_log"))


# Fixed parameters -------------------------------------------------------------

# Survival fixed effects
su_fe <- extract_fixed_pars(
  model    = mod_su_best,
  term_map = su_term_map,
  fill_missing = TRUE)


# Growth fixed effects
gr_fe <- extract_fixed_pars(
  model    = mod_gr_best,
  term_map = gr_term_map,
  fill_missing = TRUE)


# Flowering fixed effects
fl_fe <- extract_fixed_pars(
  model    = mod_fl_best,
  term_map = fl_term_map,
  fill_missing = TRUE)


# Flower number fixed effects
fln_fe <- extract_fixed_pars(
  model    = mod_fl_n_best,
  term_map = fln_term_map,
  fill_missing = TRUE)


# Recruitment fixed effects
re_fe <- extract_fixed_pars(
  model    = mod_re_best,
  term_map = re_term_map,
  fill_missing = TRUE)


# Constants -------------------------------------------------------------------
gr_var_coef <- coef(gr_var)

constants <- tibble::tribble(
  ~coefficient,      ~value,
  "recr_sz",         rc_sz$mean,
  "recr_sd",         rc_sz$sd,
  "a",               as.numeric(gr_var_coef[1]),
  "b",               as.numeric(gr_var_coef[2]),
  "L",               min(df_gr$logsize_t0, na.rm = TRUE),
  "U",               max(df_gr$logsize_t0, na.rm = TRUE),
  "mat_siz",         200,
  "fh_t0_ref",       mean(df_re$fh_t0, na.rm = TRUE),
  "fh_t0_ref_log",   log1p(mean(df_re$fh_t0, na.rm = TRUE)),
  "re_sigma",        sigma(mod_re_best)) %>%
  mutate(
    coefficient = as.character(coefficient),
    value       = as.numeric(value))


# Constant parameter table -----------------------------------------------------

pars_cons <- Reduce(
  function(...) rbind(...), list(
    su_fe, gr_fe, fl_fe, fln_fe, re_fe, constants)) %>%
  mutate(coefficient = as.character(coefficient))

check_duplicate_coefficients(pars_cons, object_name = "pars_cons")

pars_cons_wide <- coef_df_to_list(pars_cons)

pars_mean <- pars_cons_wide


# Year-varying parameters ------------------------------------------------------

# These use conditional coefficients from coef(model)$year.
# If a model has no year-level random effects, it returns an empty table.

su_out_yr <- extract_group_pars(
  model     = mod_su_best,
  group_var = "year",
  term_map  = su_term_map,
  fill_missing = FALSE)


gr_out_yr <- extract_group_pars(
  model     = mod_gr_best,
  group_var = "year",
  term_map  = gr_term_map,
  fill_missing = FALSE)


fl_out_yr <- extract_group_pars(
  model     = mod_fl_best,
  group_var = "year",
  term_map  = fl_term_map,
  fill_missing = FALSE)


fln_out_yr <- extract_group_pars(
  model     = mod_fl_n_best,
  group_var = "year",
  term_map  = fln_term_map,
  fill_missing = FALSE)


pars_var <- Reduce(
  function(...) rbind(...), list(
    su_out_yr, gr_out_yr, fl_out_yr, fln_out_yr)) %>%
  mutate(coefficient = as.character(coefficient))

check_duplicate_coefficients(pars_var, object_name = "pars_var")

pars_var_wide <- coef_df_to_list(pars_var)


# Site-varying recruitment random-effect deviations ---------------------------

re_out_site <- extract_re_ranef_devs(
  model       = mod_re_best,
  group_var   = "site",
  group_label = "site")

pars_re_site <- re_out_site %>%
  mutate(coefficient = as.character(coefficient))

check_duplicate_coefficients(
  pars_re_site,
  object_name = "pars_re_site")

pars_re_site_wide <- coef_df_to_list(pars_re_site)

# Recruitment-year random-effect deviations
re_out_year_t1 <- extract_re_ranef_devs(
  model       = mod_re_best,
  group_var   = "year_t1",
  group_label = "year_t1")

pars_re_year_t1 <- re_out_year_t1 %>%
  mutate(coefficient = as.character(coefficient))

check_duplicate_coefficients(
  pars_re_year_t1, object_name = "pars_re_year_t1")

pars_re_year_t1_wide <- coef_df_to_list(pars_re_year_t1)


# Combined parameter objects -----------------------------------------
pars_all_mean       <- pars_mean
pars_all_year       <- pars_var_wide
pars_all_site       <- pars_re_site_wide
pars_all_re_year_t1 <- pars_re_year_t1_wide


# Building the year-specific IPMs from scratch ---------------------------------
# Add missing recruitment reference constants

# The recruitment model was fitted at the site-year level.
# Therefore, inside the IPM, we use a reference site-year flowerhead production
# to convert predicted total recruits into recruits per individual.

if(is.null(pars_all_mean$fh_t0_ref)) {
  pars_all_mean$fh_t0_ref <- mean(df_re$fh_t0, na.rm = TRUE)
}

if(is.null(pars_all_mean$fh_t0_ref_log)) {
  pars_all_mean$fh_t0_ref_log <- log1p(pars_all_mean$fh_t0_ref)
}

if(is.null(pars_all_mean$re_sigma)) {
  pars_all_mean$re_sigma <- sigma(mod_re_best)
}


# Functions --------------------------------------------------------------------

# Inverse logit
inv_logit <- function(x) {
  plogis(x)
}


# Safe parameter extraction
get_par <- function(pars, par_name, default = 0) {
  
  if(!is.null(pars[[par_name]])) {
    out <- pars[[par_name]]
  } else {
    out <- default
  }
  
  return(out)
}


# Get highest polynomial order available for a vital rate
get_poly_order <- function(pars, prefix) {
  
  par_names <- names(pars)
  
  hits <- grep(paste0("^", prefix, "b[0-9]+$"), par_names, value = TRUE)
  
  if(length(hits) == 0) {
    return(0)
  }
  
  powers <- as.numeric(
    gsub(paste0("^", prefix, "b"), "", hits))
  
  powers <- powers[!is.na(powers)]
  
  return(max(powers))
}


# Generic linear predictor for size-based vital rates
vital_lp <- function(x, pars, prefix, fire = 0) {
  
  value <- get_par(
    pars = pars, par_name = paste0(prefix, "b0"), default  = 0)
  
  num_params <- get_poly_order(pars = pars, prefix = prefix)
  
  if(num_params > 0) {
    
    for(i in 1:num_params) {
      
      param_name <- paste0(prefix, "b", i)
      
      value <- value +
        get_par(pars = pars, par_name = param_name, default  = 0) * x^i
    }
  }
  
  value <- value + get_par(
    pars = pars, par_name = paste0(prefix, "bf"), default  = 0) * fire
  
  return(value)
}


# Build parameter list for mean, year-specific, or site-specific IPMs
make_ipm_pars <- function(
    pars_mean, pars_year = NULL, pars_site = NULL, pars_re_year_t1 = NULL, 
    year = NULL, site = NULL, re_year_t1 = NULL) {
  
  pars <- pars_mean
  
  # Add year-specific survival, growth, flowering, and flower-number parameters
  if(!is.null(year) && !is.null(pars_year)) {
    
    year <- as.character(year)
    
    year_hits <- grep(paste0("_", year, "$"), names(pars_year), value = TRUE)
    
    if(length(year_hits) > 0) {
      
      for(i in year_hits) {
        
        new_name <- sub(paste0("_", year, "$"), "", i)
        
        pars[[new_name]] <- pars_year[[i]]
      }
    }
  }
  
  # Add site-specific recruitment random-effect deviations
  if(!is.null(site) && !is.null(pars_site)) {
    
    site <- as.character(site)
    
    pars$re_b0 <- get_par(pars, "re_b0", 0) +
      get_par(pars_site, paste0("re_u0_site_", site), 0)
    
    pars$re_b1 <- get_par(pars, "re_b1", 0) +
      get_par(pars_site, paste0("re_u1_site_", site), 0)
  }
  
  # If not supplied directly, recruitment year is t1 = t0 + 1
  if(is.null(re_year_t1) && !is.null(year)) {
    re_year_t1 <- as.numeric(as.character(year)) + 1
  }
  
  # Add recruitment-year random-effect deviations
  # These belong to the recruit census year t1.
  if(!is.null(re_year_t1) && !is.null(pars_re_year_t1)) {
    
    re_year_t1 <- as.character(re_year_t1)
    
    pars$re_b0 <- get_par(pars, "re_b0", 0) +
      get_par(pars_re_year_t1, paste0("re_u0_year_t1_", re_year_t1), 0)
    
    pars$re_b1 <- get_par(pars, "re_b1", 0) +
      get_par(pars_re_year_t1, paste0("re_u1_year_t1_", re_year_t1), 0)
  }
  
  
  return(pars)
}


# Survival of x-sized individual to time t1
sx <- function(x, pars, fire = 0) {
  
  survival_value <- vital_lp(
    x = x, pars = pars, prefix = "surv_", fire   = fire)
  
  return(inv_logit(survival_value))
}


# Standard deviation of growth model
grow_sd <- function(x, pars) {
  
  sd_value <- pars$a * exp(pars$b * x)
  sd_value <- sqrt(sd_value)
  
  return(sd_value)
}


# Growth from size x to size y
gxy <- function(x, y, pars, fire = 0) {
  
  mean_value <- vital_lp(
    x = x, pars = pars, prefix = "grow_", fire   = fire)
  
  sd_value <- grow_sd(x = x, pars = pars)
  
  return(dnorm(y, mean = mean_value, sd = sd_value))
}


# Transition of x-sized individual to y-sized individual at time t1
pxy <- function(x, y, pars, fire = 0) {
  
  return(sx(x, pars, fire = fire) * gxy(x, y, pars, fire = fire))
}


# Probability that an x-sized individual flowers
fl_x <- function(x, pars, fire = 0) {
  
  fl_value <- vital_lp(x = x, pars = pars, prefix = "fl_", fire   = fire)
  
  return(inv_logit(fl_value))
}


# Expected flowerhead number of a flowering x-sized individual
fln_x <- function(x, pars, fire = 0) {
  
  fln_value <- vital_lp(
    x      = x, pars   = pars, prefix = "fln_", fire   = fire)
  
  fln_value <- exp(fln_value)
  
  return(fln_value)
}


# Expected flowerhead production of an x-sized individual
fhx <- function(x, pars, fire = 0) {
  
  fh_value <- fl_x(x = x, pars = pars, fire = fire) *
    fln_x(x = x, pars = pars, fire = fire)
  
  return(fh_value)
}


# Predicted total recruitment at reference site-year flowerhead production
re_total_ref <- function(pars) {
  
  re_log <- get_par(pars, "re_b0", 0) +
    get_par(pars, "re_b1", 0) *
    get_par(pars, "fh_t0_ref_log", 0)
  
  re_value <- expm1(
    re_log + 0.5 * get_par(pars, "re_sigma", 0)^2
  )
  
  re_value <- pmax(re_value, 0)
  
  return(re_value)
}


# Recruitment function
rx <- function(x, pars, fire = 0) {
  # Flowerheads produced by an x-sized individual
  fh_x <- fhx(x = x, pars = pars, fire = fire)
  
  # Predicted total recruits for a reference site-year
  re_ref <- re_total_ref(pars = pars)
  
  # Convert total recruitment into recruitment per flowerhead
  re_per_fh <- re_ref /
    pmax(get_par(pars, "fh_t0_ref", 0), .Machine$double.eps)
  
  # Expected recruits contributed by an individual of size x
  re_value <- fh_x * re_per_fh
  
  return(re_value)
}


# Recruit size distribution
recr_y <- function(y, pars) {
  
  p_y <- dnorm(y, mean = pars$recr_sz, sd= pars$recr_sd)
  
  return(p_y)
}


# Fertility function
fy <- function(y, x, pars, h, fire = 0) {
  # Recruit size distribution
  Pvec <- recr_y(y = y, pars = pars)
  
  # Correct for eviction of recruits outside the mesh
  Pvec <- Pvec / sum(Pvec * h)
  
  # Expected recruits produced by individuals of size x
  Rvec <- rx(x = x, pars = pars, fire = fire)
  
  # Fertility matrix
  Fmat <- outer(Pvec, Rvec) * h
  
  return(Fmat)
}


# Kernel -----------------------------------------------------------------------

kernel <- function(pars,
                   fire = 0) {
  
  # number of bins over which to integrate
  n   <- pars$mat_siz
  # lower limit of integration
  L   <- pars$L
  # upper limit of integration
  U   <- pars$U
  # bin size
  h   <- (U - L) / n
  # lower boundaries of bins
  b   <- L + c(0:n) * h
  # midpoints of bins
  y   <- 0.5 * (b[1:n] + b[2:(n + 1)])
  
  # Fertility matrix
  Fmat <- matrix(0, n, n)
  Fmat <- fy(y = y, x = y, pars = pars, h = h, fire = fire)
  
  # Survival vector
  Smat <- c()
  Smat <- sx(x = y, pars = pars, fire = fire)
  
  # Growth matrix
  Gmat <- matrix(0, n, n)
  Gmat[] <- t(
    outer(y, y, Vectorize(function(x, y) {
      gxy(x = x, y = y, pars = pars, fire = fire)
        }))) * h
  
  # Growth/survival transition matrix
  Tmat <- matrix(0, n, n)
  
  # Correct for eviction of small individuals
  for(i in 1:(n / 2)) {
    Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    Tmat[, i]  <- Gmat[, i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for(i in (n / 2 + 1):n) {
    Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    Tmat[, i]  <- Gmat[, i] * Smat[i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrices
  k_yx <- Fmat + Tmat
  
  return(list(
    k_yx    = k_yx,
    Fmat    = Fmat,
    Tmat    = Tmat,
    Gmat    = Gmat,
    meshpts = y))
}


# Mean population growth rate --------------------------------------------------

lambda_ipm <- function(i, fire = 0) {
  return(Re(eigen(kernel(i, fire = fire)$k_yx)$value[1]))
}


# Mean IPM ---------------------------------------------------------------------

pars_mean <- pars_all_mean

lambda_ipm(pars_mean, fire = 0)
lambda_ipm(pars_mean, fire = 1)


# Year-specific IPMs -----------------------------------------------------------

lambda_ipm_year <- function(year,
                            fire = 0,
                            site = NULL,
                            re_year_t1 = NULL) {
  pars_i <- make_ipm_pars(
    pars_mean       = pars_all_mean,
    pars_year       = pars_all_year,
    pars_site       = pars_all_site,
    pars_re_year_t1 = pars_all_re_year_t1,
    year            = year,
    site            = site,
    re_year_t1      = re_year_t1)
  return(lambda_ipm(pars_i, fire = fire))
}


# Example: year-specific lambda
lambda_ipm_year(year = 1999, fire = 0)
lambda_ipm_year(year = 1999, fire = 1)


# Example: year-specific lambda with site-specific recruitment
lambda_ipm_year(year = 1999, fire = 0, site = 1)
lambda_ipm_year(year = 1999, fire = 1, site = 1)


# Build all year-specific lambdas ---------------------------------------------

ipm_years <- sort(unique(
  as.numeric(
    gsub("surv_b0_", "",
      grep("^surv_b0_", names(pars_all_year), value = TRUE)))))

lambda_year <- data.frame(
  year           = ipm_years,
  lambda_no_fire = sapply(
    ipm_years, function(i) lambda_ipm_year(year = i, fire = 0, site = NULL)),
  lambda_fire    = sapply(
    ipm_years, function(i) lambda_ipm_year(year = i, fire = 1, site = NULL)))

lambda_year


# Observed PGR vs modeled projected lambda, site-paired ------------------------

# Keep NULL for population-level recruitment.
# Set TRUE only if you explicitly want site-specific recruitment coefficients.
use_site_specific_recruitment <- FALSE


# Site-year fire lookup --------------------------------------------------------

fire_lookup_site <- df %>%
  group_by(site, year) %>%
  summarise(
    fire = if_else(any(fire == "Fire", na.rm = TRUE), "Fire", "No fire"),
    .groups = "drop") %>%
  mutate(
    fire = factor(fire, levels = c("No fire", "Fire")),
    fire_num = if_else(fire == "Fire", 1, 0))


get_fire_site_y <- function(site_i, yr) {
  
  fire_y <- fire_lookup_site %>%
    filter(
      site == site_i,
      year == yr) %>%
    pull(fire_num)
  
  if(length(fire_y) == 0) {
    fire_y <- 0
  }
  return(fire_y)
}


# Site-year observed counts ----------------------------------------------------

df_counts_site <- df %>%
  filter(
    !is.na(logsize_t0),
    is.finite(logsize_t0)) %>%
  group_by(site, year) %>%
  summarise(
    n = n(),
    .groups = "drop")


# Keep only site-years where both t and t+1 are observed -----------------------

df_obs_pgr_site <- df_counts_site %>%
  rename(n_t0 = n) %>%
  left_join(
    df_counts_site %>%
      transmute(
        site,
        year = year - 1,
        n_t1 = n),
    by = c("site", "year")) %>%
  filter(!is.na(n_t1)) %>%
  left_join(
    fire_lookup_site,
    by = c("site", "year"))


# Initial site-specific size distribution --------------------------------------

make_initial_n_site <- function(year0, site_i, pars) {
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  
  df0 <- df %>%
    filter(
      year == year0,
      site == site_i,
      !is.na(logsize_t0),
      is.finite(logsize_t0))
  
  adult_counts <- hist(
    pmin(pmax(df0$logsize_t0, L), U),
    breaks = seq(L, U, length.out = n + 1),
    plot = FALSE,
    include.lowest = TRUE)$counts
  
  adult_density <- adult_counts / h
  
  return(adult_density)}


# Project one site for one year ------------------------------------------------

project_one_site_year <- function(yr, site_i) {
  
  site_re <- if(use_site_specific_recruitment) site_i else NULL
  
  pars_y <- make_ipm_pars(
    pars_mean       = pars_all_mean,
    pars_year       = pars_all_year,
    pars_site       = pars_all_site,
    pars_re_year_t1 = pars_all_re_year_t1,
    year            = yr,
    site            = site_re,
    re_year_t1      = yr + 1)
  
  fire_y <- get_fire_site_y(
    site_i = site_i,
    yr     = yr)
  
  n_obs <- make_initial_n_site(
    year0  = yr,
    site_i = site_i,
    pars   = pars_y)
  
  K <- kernel(
    pars = pars_y,
    fire = fire_y)$k_yx
  
  h <- (pars_y$U - pars_y$L) / pars_y$mat_siz
  
  n_proj <- K %*% n_obs
  
  out <- data.frame(
    year = yr,
    site = site_i,
    fire_num = fire_y,
    n_obs = sum(n_obs) * h,
    n_proj = sum(n_proj) * h,
    asym_lambda = Re(eigen(K)$values[1]),
    proj_lambda = as.numeric((sum(n_proj) * h) / (sum(n_obs) * h)))
  
  return(out)
}


# Project all comparable site-years -------------------------------------------

df_proj_site <- bind_rows(
  lapply(seq_len(nrow(df_obs_pgr_site)), function(i) {
    
    project_one_site_year(
      yr     = df_obs_pgr_site$year[i],
      site_i = df_obs_pgr_site$site[i])
  }))


# Combine site-level observed and projected values -----------------------------

df_compare_site <- df_obs_pgr_site %>%
  left_join(
    df_proj_site,
    by = c("year", "site", "fire_num")) %>%
  mutate(
    obs_pgr = n_t1 / n_t0)


# Whole-population annual comparison ------------------------------------------

df_compare <- df_compare_site %>%
  group_by(year) %>%
  summarise(
    n_t0 = sum(n_t0, na.rm = TRUE),
    n_t1 = sum(n_t1, na.rm = TRUE),
    obs_pgr = n_t1 / n_t0,
    n_obs_model = sum(n_obs, na.rm = TRUE),
    n_proj_model = sum(n_proj, na.rm = TRUE),
    proj_lambda = n_proj_model / n_obs_model,
    asym_lambda = weighted.mean(
      asym_lambda,
      w = n_obs,
      na.rm = TRUE),
    fire = if_else(any(fire == "Fire", na.rm = TRUE), "Fire", "No fire"),
    .groups = "drop") %>%
  mutate(
    fire = factor(fire, levels = c("No fire", "Fire")))


# Long format for plotting -----------------------------------------------------

df_plot <- df_compare %>%
  select(
    year,
    obs_pgr,
    asym_lambda,
    proj_lambda,
    fire) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = "lambda_type",
    values_to = "lambda") %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = "Abundance-weighted asymptotic lambda",
      proj_lambda = "Projected lambda from observed site size distributions"))


# Plot -------------------------------------------------------------------------

g_mod_vs_obs <- ggplot(
  df_plot,
  aes(x = lambda, y = obs_pgr, color = fire)) +
  scale_color_manual(values = c("Fire" = "#FF5733", "No fire" = "black")) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = "free_x") +
  labs(
    title = "Observed population growth vs modeled lambda",
    x     = expression("Modeled " * lambda),
    y     = "Observed population growth rate",
    color = "Fire") +
  theme_classic()

g_mod_vs_obs


g_mod_vs_obs2 <- ggplot(
  df_plot,
  aes(x = lambda, y = obs_pgr, color = fire)) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Fire" = "#FF5733", "No fire" = "black")) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = "free_x") +
  labs(
    title = "Observed population growth vs modeled lambda",
    x     = expression("Modeled " * lambda * " (log scale)"),
    y     = "Observed population growth rate (log scale)",
    color = "Fire") +
  theme_classic()

g_mod_vs_obs2


# Summary statistics -----------------------------------------------------------

df_compare %>%
  summarise(
    arithmetic_mean_obs_pgr = mean(obs_pgr, na.rm = TRUE),
    geometric_mean_obs_pgr  = exp(mean(log(obs_pgr), na.rm = TRUE)),
    arithmetic_mean_asym_lambda = mean(asym_lambda, na.rm = TRUE),
    geometric_mean_asym_lambda  = exp(mean(log(asym_lambda), na.rm = TRUE)),
    arithmetic_mean_proj_lambda = mean(proj_lambda, na.rm = TRUE),
    geometric_mean_proj_lambda  = exp(mean(log(proj_lambda), na.rm = TRUE)))


# Inspect comparison table -----------------------------------------------------

df_compare %>%
  print(n = 100)



# Site-year observed versus modeled lambda ------------------------------------
use_site_specific_recruitment <- TRUE

df_plot_site <- df_compare_site %>%
  select(
    site, year, fire, obs_pgr,
    asym_lambda, proj_lambda) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = "lambda_type",
    values_to = "lambda") %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = "Asymptotic lambda",
      proj_lambda = "Projected lambda"))

g_mod_vs_obs_site <- ggplot(
  df_plot_site,
  aes(x = lambda, y = obs_pgr, color = factor(year))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed") +
  facet_wrap(~ lambda_type, scales = "free") +
  scale_color_viridis_d() +
  labs(
    title = "Observed versus modeled site-year population growth",
    subtitle = "Site and year-specific recruitment effects",
    x = expression("Modeled " * lambda),
    y = expression("Observed " * lambda),
    color = "Year") +
  theme_bw()

g_mod_vs_obs_site


# Investigation into observed population growth rates above 2.5 
#  for site and year combinations
extreme_pairs <- df_compare_site %>%
  filter(obs_pgr > 2.5) %>%
  select(site, year)

plant_presence <- df %>%
  filter(!is.na(logsize_t0), is.finite(logsize_t0)) %>%
  distinct(site, year, plant_id)

df_id_transition <- full_join(
  plant_presence %>%
    mutate(in_t0 = TRUE),
  plant_presence %>%
    transmute(
      site,
      year = year - 1,
      plant_id,
      in_t1 = TRUE),
  by = c("site", "year", "plant_id")) %>%
  mutate(
    in_t0 = replace_na(in_t0, FALSE),
    in_t1 = replace_na(in_t1, FALSE)) %>%
  semi_join(extreme_pairs, by = c("site", "year")) %>%
  group_by(site, year) %>%
  summarise(
    n_t0 = sum(in_t0),
    n_t1 = sum(in_t1),
    persisted = sum(in_t0 & in_t1),
    disappeared = sum(in_t0 & !in_t1),
    new_ids = sum(!in_t0 & in_t1),
    .groups = "drop") %>%
  mutate(obs_pgr = n_t1 / n_t0) %>%
  arrange(desc(obs_pgr))

df_id_transition


df_new_ids <- df %>%
  filter(!is.na(logsize_t0), is.finite(logsize_t0)) %>%
  transmute(
    site,
    year = year - 1,
    plant_id,
    recruit) %>%
  anti_join(
    df %>%
      filter(!is.na(logsize_t0), is.finite(logsize_t0)) %>%
      distinct(site, year, plant_id),
    by = c("site", "year", "plant_id")) %>%
  semi_join(extreme_pairs, by = c("site", "year")) %>%
  group_by(site, year) %>%
  summarise(
    new_ids = n(),
    recruits = sum(recruit == 1, na.rm = TRUE),
    non_recruits = sum(recruit == 0, na.rm = TRUE),
    recruit_unknown = sum(is.na(recruit)),
    .groups = "drop") %>%
  arrange(desc(new_ids))

df_new_ids
