# =============================================================================
# Short IPM for Chrysopsis highlandsensis
# Archbold demographic data
#
# Workflow:
#   1. Load and clean demographic data.
#   2. Fit vital-rate models:
#        survival
#        growth
#        flowering probability
#        flower-head number conditional on flowering
#   3. Build an IPM with:
#        one seed-bank state
#        continuous aboveground size states
#   4. Validate using matched-site observed PGR vs projected lambda.
#   5. Make all model-effect plots at the end.
#
# Important:
#   The seed bank here is an EFFECTIVE recruit-producing seed bank.
#   It is not a direct estimate of literal soil seed density.
# =============================================================================


# Setup ------------------------------------------------------------------------

set.seed(100)
options(stringsAsFactors = FALSE)

# Load MASS before tidyverse so dplyr functions keep priority.
library(MASS)
library(tidyverse)
library(janitor)
library(lme4)
library(patchwork)

v_head <- "archbold"
v_species <- "Chrysopsis highlandsensis"

v_sp_abb <- tolower(
  gsub(
    " ",
    "",
    paste(substr(unlist(strsplit(v_species, " ")), 1, 2), collapse = "")
  )
)

dir_data <- file.path(v_head, "data", v_sp_abb)

# Years used for fitting and validation.
# compare_years needs year + 1 to exist for observed population growth.
fit_years <- 1999:2023
compare_years <- 1999:2023

# Number of mesh cells in the continuous size distribution.
# 200 is fine enough for smooth IPM integration but still fast.
mat_siz <- 200


# Biological constants ---------------------------------------------------------

# Annual emergence / field germination from seed bank.
# Paper reports field germination experiments from 1.0% to 8.6%.
# We use the conservative lower bound as the default.
germ_p_default <- 0.01

# Dormant seed survival / persistence from one year to the next.
# This is NOT 1 - germ_p.
# The paper observed germination up to 5 years after sowing, but late
# germination was infrequent. So this value is an assumption for the
# effective seed bank, not a directly measured parameter.
seed_surv_default <- 0.60

# Intact seeds per flowering head.
# Paper estimate: mean intact seeds per head = 50.9.
# Do NOT use 3181 here, because 3181 is already approximately:
#   mean heads per reproductive plant * intact seeds per head.
# In this IPM, head number is predicted separately by mod_heads.
intact_seeds_per_head <- 50.9

# Required for a valid seed-bank transition.
stopifnot(seed_surv_default + germ_p_default <= 1)


# Read data --------------------------------------------------------------------

df_og <- read_csv(
  file.path(dir_data, "chrysopsis_highlandsensis_data.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

df_fire <- read_csv(
  file.path(dir_data, "chrysopsis_highlandsensis_fire.csv"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  rename(
    year = burn_yr,
    fire = treatment
  ) %>%
  dplyr::select(site, year, fire) %>%
  mutate(
    fire = case_when(
      fire == "burn" ~ "Fire",
      TRUE ~ "No fire"
    )
  )


# Clean demographic data -------------------------------------------------------

df_all <- df_og %>%
  rename(
    plant_id = identifier,
    year = year0,
    survival = survival_1
  ) %>%
  mutate(
    plant_id = factor(plant_id),
    
    # Recruits are individuals with astg == 1.
    recruit = as.integer(!is.na(astg) & astg == 1)
  ) %>%
  group_by(site, plant_id, year) %>%
  summarise(
    survives = if (all(is.na(survival))) NA_real_ else min(survival, na.rm = TRUE),
    size_t0  = if (all(is.na(dia)))      NA_real_ else max(dia, na.rm = TRUE),
    size_t1  = if (all(is.na(dia_1)))    NA_real_ else max(dia_1, na.rm = TRUE),
    
    # Crucial choice:
    # If all hd values are NA, treat this as 0 heads / non-flowering.
    # This keeps non-flowering plants in the data instead of dropping them.
    flower_raw = if (all(is.na(hd))) 0 else max(hd, na.rm = TRUE),
    
    recruit = max(recruit, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    fl_nr = replace_na(flower_raw, 0),
    flower = as.integer(fl_nr > 0),
    
    logsize_t0 = log(size_t0),
    logsize_t1 = log(size_t1)
  ) %>%
  left_join(df_fire, by = c("site", "year")) %>%
  mutate(
    fire = replace_na(fire, "No fire"),
    fire = factor(fire, levels = c("No fire", "Fire"))
  )

df <- df_all %>%
  filter(year %in% fit_years)


# Fire proportion by year ------------------------------------------------------

# Do not treat the whole population as burned when only one site burned.
# For IPM validation, fire is represented as the proportion of observed
# plants in burned sites.
fire_prop_by_year <- df %>%
  filter(is.finite(logsize_t0)) %>%
  group_by(year) %>%
  summarise(
    fire_prop = mean(fire == "Fire", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    fire_prop = replace_na(fire_prop, 0)
  )


# Vital-rate models ------------------------------------------------------------

# Survival ---------------------------------------------------------------------

# Survival is a binomial response:
#   0 = dead
#   1 = alive
#
# Predictors:
#   logsize_t0: larger plants generally survive better
#   flower: captures post-flowering survival cost
#   fire: direct fire effect
#   (1 | year): year-specific survival variation
#
# The flower term is important because the species is largely monocarpic.
df_su <- df %>%
  filter(
    !is.na(survives),
    is.finite(logsize_t0),
    !is.na(fire)
  ) %>%
  mutate(year = factor(year))

mod_su <- glmer(
  survives ~ logsize_t0 + flower + fire + (1 | year),
  data = df_su,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)


# Growth -----------------------------------------------------------------------

# Growth is conditional on being alive and observed at t + 1.
# Response is log size at t + 1.
# Random year-specific intercept and slope allow annual variation in growth.
df_gr <- df %>%
  filter(
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1),
    !is.na(fire)
  ) %>%
  mutate(year = factor(year))

mod_gr <- lmer(
  logsize_t1 ~ logsize_t0 + fire + (logsize_t0 | year),
  data = df_gr
)


# Growth variance --------------------------------------------------------------

# The IPM needs a full growth distribution, not only mean growth.
# We model residual variance as a function of fitted growth.
# If this fails, we use constant residual variance.
gr_var_dat <- tibble(
  mu = fitted(mod_gr),
  r2 = resid(mod_gr)^2
)

gr_var_fit <- try(
  nls(
    r2 ~ a * exp(b * mu),
    data = gr_var_dat,
    start = list(a = mean(gr_var_dat$r2, na.rm = TRUE), b = 0)
  ),
  silent = TRUE
)

if (inherits(gr_var_fit, "try-error")) {
  gr_a <- mean(gr_var_dat$r2, na.rm = TRUE)
  gr_b <- 0
} else {
  gr_a <- as.numeric(coef(gr_var_fit)["a"])
  gr_b <- as.numeric(coef(gr_var_fit)["b"])
}


# Flowering probability --------------------------------------------------------

# Flowering is binary:
#   0 = did not flower
#   1 = flowered
#
# We use glm rather than glmer here because year-random flowering fits were
# unstable. Year-specific variation still enters the IPM through survival
# and growth.
df_fl <- df %>%
  filter(
    is.finite(logsize_t0),
    !is.na(fire)
  ) %>%
  mutate(year = factor(year))

mod_fl <- glm(
  flower ~ logsize_t0 + fire,
  data = df_fl,
  family = binomial
)


# Flower heads conditional on flowering ---------------------------------------

# Flower-head number is a count, conditional on flowering.
# A negative-binomial GLM is appropriate for overdispersed count data.
# One non-integer head value is treated as a data-entry issue and rounded.
df_heads <- df_fl %>%
  filter(flower == 1, !is.na(fl_nr), fl_nr > 0) %>%
  mutate(
    fl_nr = round(fl_nr),
    year = factor(year),
    fire = factor(fire, levels = levels(df_fl$fire))
  )

mod_heads <- MASS::glm.nb(
  fl_nr ~ logsize_t0 + fire,
  data = df_heads
)


# Recruit size distribution ----------------------------------------------------

# Recruits enter the aboveground population with a log-size distribution.
# This distribution is estimated directly from observed recruits.
df_recr <- df %>%
  filter(recruit == 1, size_t0 > 0, is.finite(logsize_t0))

recr_sz <- mean(df_recr$logsize_t0, na.rm = TRUE)
recr_sd <- sd(df_recr$logsize_t0, na.rm = TRUE)


# Seed-bank calibration --------------------------------------------------------

# We do not know actual soil seed density.
# Instead we estimate the effective standing seed bank needed to produce
# the observed mean number of recruits:
#
#   seedbank_size = mean annual recruits / germ_p
#
# This is an effective recruit-producing seed bank, not literal seed density.

annual_seed_obs <- df %>%
  group_by(year) %>%
  summarise(
    recruits = sum(recruit == 1, na.rm = TRUE),
    
    # Raw intact seed production per year:
    # observed flower heads * paper estimate of intact seeds per head.
    raw_intact_seeds = sum(fl_nr, na.rm = TRUE) * intact_seeds_per_head,
    
    .groups = "drop"
  )

mean_rec <- mean(annual_seed_obs$recruits, na.rm = TRUE)
mean_raw_intact_seed_input <- mean(annual_seed_obs$raw_intact_seeds, na.rm = TRUE)

calibrate_seedbank <- function(germ_p_i, seed_surv_i) {
  
  stopifnot(seed_surv_i + germ_p_i <= 1)
  
  # Effective seed-bank size implied by observed recruitment.
  seedbank_size_i <- mean_rec / germ_p_i
  
  # Equilibrium effective input:
  #
  #   B[t + 1] = seed_surv * B[t] + effective_input
  #   B = seed_surv * B + effective_input
  #   effective_input = B * (1 - seed_surv)
  #
  # This replaces seeds lost to emergence and death.
  effective_seed_input_mean_i <- seedbank_size_i * (1 - seed_surv_i)
  
  # Convert raw intact seed production into effective seed-bank entrants.
  # This absorbs losses from predation, nonviability, microsite failure,
  # dispersal outside plots, etc.
  seedbank_entry_prob_i <- effective_seed_input_mean_i /
    mean_raw_intact_seed_input
  
  seedbank_entry_prob_i <- min(max(seedbank_entry_prob_i, 0), 1)
  
  # Reconstruct year-specific effective seed bank:
  #
  #   B[t + 1] = seed_surv * B[t] + raw_intact_seeds[t] * entry_prob
  seedbank_by_year_i <- tibble(
    year = fit_years,
    seedbank = NA_real_
  )
  
  B <- seedbank_size_i
  
  for (i in seq_along(fit_years)) {
    
    yr <- fit_years[i]
    seedbank_by_year_i$seedbank[i] <- B
    
    raw_input_yr <- annual_seed_obs %>%
      filter(year == yr) %>%
      pull(raw_intact_seeds)
    
    if (length(raw_input_yr) == 0 || is.na(raw_input_yr)) {
      raw_input_yr <- mean_raw_intact_seed_input
    }
    
    effective_input_yr <- raw_input_yr * seedbank_entry_prob_i
    
    B <- seed_surv_i * B + effective_input_yr
  }
  
  list(
    germ_p = germ_p_i,
    seed_surv = seed_surv_i,
    seedbank_size = seedbank_size_i,
    seedbank_entry_prob = seedbank_entry_prob_i,
    seedbank_by_year = seedbank_by_year_i
  )
}

seed_par_default <- calibrate_seedbank(
  germ_p_i = germ_p_default,
  seed_surv_i = seed_surv_default
)

cat("\nSeed-bank default assumptions\n")
cat("mean recruits/year:       ", mean_rec, "\n")
cat("germ_p:                   ", seed_par_default$germ_p, "\n")
cat("seed_surv:                ", seed_par_default$seed_surv, "\n")
cat("seedbank_size:            ", seed_par_default$seedbank_size, "\n")
cat("intact_seeds_per_head:    ", intact_seeds_per_head, "\n")
cat("seedbank_entry_prob:      ", seed_par_default$seedbank_entry_prob, "\n\n")


# IPM helper functions ---------------------------------------------------------

fire_fac <- function(fire_num) {
  factor(
    if_else(fire_num == 1, "Fire", "No fire"),
    levels = levels(df$fire)
  )
}

make_newdata <- function(x, fire_num, year = NULL, flower = NULL) {
  
  out <- tibble(
    logsize_t0 = x,
    fire = fire_fac(fire_num)
  )
  
  if (!is.null(year)) {
    out$year <- factor(year, levels = levels(df_su$year))
  }
  
  if (!is.null(flower)) {
    out$flower <- flower
  }
  
  out
}

pred_su <- function(x, flower_state, fire_num, year = NULL,
                    year_specific = TRUE) {
  
  nd <- make_newdata(
    x = x,
    fire_num = fire_num,
    year = year,
    flower = flower_state
  )
  
  predict(
    mod_su,
    newdata = nd,
    type = "response",
    re.form = if (year_specific) NULL else NA,
    allow.new.levels = TRUE
  )
}

pred_growth_mu <- function(x, fire_num, year = NULL,
                           year_specific = TRUE) {
  
  nd <- make_newdata(
    x = x,
    fire_num = fire_num,
    year = year
  )
  
  predict(
    mod_gr,
    newdata = nd,
    re.form = if (year_specific) NULL else NA,
    allow.new.levels = TRUE
  )
}

pred_growth_sd <- function(mu) {
  sqrt(gr_a * exp(gr_b * mu))
}

pred_flower <- function(x, fire_num) {
  
  nd <- tibble(
    logsize_t0 = x,
    fire = fire_fac(fire_num)
  )
  
  predict(
    mod_fl,
    newdata = nd,
    type = "response"
  )
}

pred_heads <- function(x, fire_num) {
  
  nd <- tibble(
    logsize_t0 = x,
    fire = fire_fac(fire_num)
  )
  
  predict(
    mod_heads,
    newdata = nd,
    type = "response"
  )
}

recruit_size_density <- function(y) {
  dnorm(y, mean = recr_sz, sd = recr_sd)
}

lambda_from_K <- function(K) {
  Re(eigen(K)$values[1])
}


# Mesh -------------------------------------------------------------------------

L <- min(
  df_gr$logsize_t0,
  df_gr$logsize_t1,
  df_recr$logsize_t0,
  na.rm = TRUE
) - 0.1

U <- max(
  df_gr$logsize_t0,
  df_gr$logsize_t1,
  df_recr$logsize_t0,
  na.rm = TRUE
) + 0.1

h <- (U - L) / mat_siz

breaks <- seq(L, U, length.out = mat_siz + 1)
meshpts <- 0.5 * (breaks[-1] + breaks[-length(breaks)])


# IPM kernel -------------------------------------------------------------------

# State vector:
#   [1]      seed bank
#   [2:n+1]  aboveground plants distributed by log size
#
# Blocks:
#   S11: seed bank -> seed bank
#   S12: aboveground reproductive plants -> seed bank
#   S21: seed bank -> recruit size distribution
#   S22: aboveground plant survival/growth

kernel_one_fire <- function(fire_num,
                            year = NULL,
                            year_specific = TRUE,
                            germ_p_i = germ_p_default,
                            seed_surv_i = seed_surv_default,
                            seedbank_entry_prob_i = seed_par_default$seedbank_entry_prob) {
  
  y <- meshpts
  n <- length(y)
  
  # Flowering probability by size.
  fl <- pred_flower(y, fire_num)
  
  # Survival if not flowering.
  s_nf <- pred_su(
    x = y,
    flower_state = 0,
    fire_num = fire_num,
    year = year,
    year_specific = year_specific
  )
  
  # Survival if flowering.
  # This captures the post-flowering survival cost.
  s_f <- pred_su(
    x = y,
    flower_state = 1,
    fire_num = fire_num,
    year = year,
    year_specific = year_specific
  )
  
  # Total survival-growth probability.
  s_total <- (1 - fl) * s_nf + fl * s_f
  
  # Growth distribution.
  mu_g <- pred_growth_mu(
    x = y,
    fire_num = fire_num,
    year = year,
    year_specific = year_specific
  )
  
  sd_g <- pred_growth_sd(mu_g)
  
  Gmat <- matrix(0, n, n)
  
  for (j in seq_len(n)) {
    
    Gmat[, j] <- dnorm(y, mean = mu_g[j], sd = sd_g[j]) * h
    
    # Eviction correction.
    Gmat[1, j] <- Gmat[1, j] + pnorm(L, mean = mu_g[j], sd = sd_g[j])
    Gmat[n, j] <- Gmat[n, j] + 1 - pnorm(U, mean = mu_g[j], sd = sd_g[j])
  }
  
  Tmat <- sweep(Gmat, 2, s_total, "*")
  
  # Seed bank remains seed bank.
  S11 <- matrix(seed_surv_i, 1, 1)
  
  # Aboveground plants produce effective seed-bank entrants.
  S12 <- matrix(
    fl *
      pred_heads(y, fire_num) *
      intact_seeds_per_head *
      seedbank_entry_prob_i *
      h,
    nrow = 1
  )
  
  # Seed bank produces recruits.
  # No *h here because this is a density over recruit size.
  S21 <- matrix(
    germ_p_i * recruit_size_density(y),
    ncol = 1
  )
  
  K <- rbind(
    cbind(S11, S12),
    cbind(S21, Tmat)
  )
  
  list(
    k_yx = K,
    Tmat = Tmat,
    meshpts = y,
    h = h
  )
}

kernel_mixed_fire <- function(year,
                              fire_prop,
                              year_specific = TRUE,
                              germ_p_i = germ_p_default,
                              seed_surv_i = seed_surv_default,
                              seedbank_entry_prob_i = seed_par_default$seedbank_entry_prob) {
  
  K0 <- kernel_one_fire(
    fire_num = 0,
    year = year,
    year_specific = year_specific,
    germ_p_i = germ_p_i,
    seed_surv_i = seed_surv_i,
    seedbank_entry_prob_i = seedbank_entry_prob_i
  )$k_yx
  
  K1 <- kernel_one_fire(
    fire_num = 1,
    year = year,
    year_specific = year_specific,
    germ_p_i = germ_p_i,
    seed_surv_i = seed_surv_i,
    seedbank_entry_prob_i = seedbank_entry_prob_i
  )$k_yx
  
  (1 - fire_prop) * K0 + fire_prop * K1
}


# Matched-site validation ------------------------------------------------------

# Matched-site comparison avoids artificial PGR jumps caused by sites entering
# or leaving the data. This is the validation target for the IPM.

counts_site <- df_all %>%
  filter(
    year %in% c(compare_years, compare_years + 1),
    is.finite(logsize_t0)
  ) %>%
  group_by(site, year) %>%
  summarise(n = n(), .groups = "drop")

matched_site_years <- counts_site %>%
  rename(year_t = year, n_t0 = n) %>%
  mutate(year_t1 = year_t + 1) %>%
  inner_join(
    counts_site %>%
      rename(year_t1 = year, n_t1 = n),
    by = c("site", "year_t1")
  ) %>%
  filter(year_t %in% compare_years)

df_obs_pgr_matched <- matched_site_years %>%
  group_by(year = year_t) %>%
  summarise(
    n_obs_t0 = sum(n_t0),
    n_obs_t1 = sum(n_t1),
    obs_pgr = n_obs_t1 / n_obs_t0,
    n_sites = n_distinct(site),
    sites = paste(sort(unique(site)), collapse = ", "),
    .groups = "drop"
  )

fire_prop_matched <- df %>%
  semi_join(
    matched_site_years %>%
      transmute(year = year_t, site),
    by = c("year", "site")
  ) %>%
  filter(is.finite(logsize_t0)) %>%
  group_by(year) %>%
  summarise(
    fire_prop = mean(fire == "Fire", na.rm = TRUE),
    .groups = "drop"
  )

obs_survivors <- df_all %>%
  semi_join(
    matched_site_years %>%
      transmute(year = year_t, site),
    by = c("year", "site")
  ) %>%
  filter(is.finite(logsize_t0), !is.na(survives)) %>%
  group_by(year) %>%
  summarise(
    obs_survivors = sum(survives == 1, na.rm = TRUE),
    .groups = "drop"
  )

obs_recruits_next <- df_all %>%
  mutate(year_t = year - 1) %>%
  semi_join(
    matched_site_years %>%
      dplyr::select(year_t, site),
    by = c("year_t", "site")
  ) %>%
  filter(recruit == 1, is.finite(logsize_t0)) %>%
  group_by(year = year_t) %>%
  summarise(
    obs_recruits_next = n(),
    .groups = "drop"
  )

make_initial_n_matched <- function(year0, seed_par) {
  
  sites_y <- matched_site_years %>%
    filter(year_t == year0) %>%
    pull(site) %>%
    unique()
  
  adult_df <- df_all %>%
    filter(
      year == year0,
      site %in% sites_y,
      is.finite(logsize_t0)
    )
  
  adult_counts <- hist(
    pmin(pmax(adult_df$logsize_t0, L), U),
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE
  )$counts
  
  adult_density <- adult_counts / h
  
  B0_total <- seed_par$seedbank_by_year %>%
    filter(year == year0) %>%
    pull(seedbank)
  
  if (length(B0_total) == 0 || is.na(B0_total)) {
    B0_total <- seed_par$seedbank_size
  }
  
  total_adults_year <- df_all %>%
    filter(year == year0, is.finite(logsize_t0)) %>%
    nrow()
  
  frac_matched <- ifelse(
    total_adults_year > 0,
    nrow(adult_df) / total_adults_year,
    0
  )
  
  # Approximate matched-site seed bank by matched-site aboveground fraction.
  B0_matched <- B0_total * frac_matched
  
  c(seedbank = B0_matched, adult_density)
}

project_one_year_matched <- function(yr, seed_par) {
  
  fire_prop_y <- fire_prop_matched %>%
    filter(year == yr) %>%
    pull(fire_prop)
  
  if (length(fire_prop_y) == 0 || is.na(fire_prop_y)) {
    fire_prop_y <- 0
  }
  
  K <- kernel_mixed_fire(
    year = yr,
    fire_prop = fire_prop_y,
    year_specific = TRUE,
    germ_p_i = seed_par$germ_p,
    seed_surv_i = seed_par$seed_surv,
    seedbank_entry_prob_i = seed_par$seedbank_entry_prob
  )
  
  n0 <- make_initial_n_matched(yr, seed_par)
  n1 <- K %*% n0
  
  adult_from_seedbank <- sum(K[-1, 1] * n0[1]) * h
  adult_from_plants <- sum((K[-1, -1, drop = FALSE] %*% n0[-1])) * h
  
  adult_n0 <- sum(n0[-1]) * h
  adult_n1 <- sum(n1[-1]) * h
  
  tibble(
    year = yr,
    asym_lambda = lambda_from_K(K),
    proj_lambda = adult_n1 / adult_n0,
    proj_from_seedbank = adult_from_seedbank,
    proj_from_plants = adult_from_plants,
    fire_prop = fire_prop_y,
    seedbank_start = n0[1],
    adult_start = adult_n0
  )
}

run_matched_projection <- function(germ_p_i, seed_surv_i) {
  
  seed_par_i <- calibrate_seedbank(
    germ_p_i = germ_p_i,
    seed_surv_i = seed_surv_i
  )
  
  df_lams_i <- purrr::map_dfr(
    df_obs_pgr_matched$year,
    ~ project_one_year_matched(.x, seed_par_i)
  )
  
  df_obs_pgr_matched %>%
    left_join(df_lams_i, by = "year") %>%
    left_join(obs_survivors, by = "year") %>%
    left_join(obs_recruits_next, by = "year") %>%
    mutate(
      obs_recruits_next = replace_na(obs_recruits_next, 0),
      obs_n1_decomp = obs_survivors + obs_recruits_next,
      obs_pgr_decomp = obs_n1_decomp / n_obs_t0,
      
      germ_p = germ_p_i,
      seed_surv = seed_surv_i,
      seedbank_size = seed_par_i$seedbank_size,
      seedbank_entry_prob = seed_par_i$seedbank_entry_prob
    ) %>%
    arrange(year)
}

df_compare_final <- run_matched_projection(
  germ_p_i = germ_p_default,
  seed_surv_i = seed_surv_default
)

cat("\nDefault matched-site summary\n")
print(
  df_compare_final %>%
    summarise(
      geom_mean_obs = exp(mean(log(obs_pgr), na.rm = TRUE)),
      geom_mean_obs_decomp = exp(mean(log(obs_pgr_decomp), na.rm = TRUE)),
      geom_mean_asym = exp(mean(log(asym_lambda), na.rm = TRUE)),
      geom_mean_proj = exp(mean(log(proj_lambda), na.rm = TRUE)),
      mean_abs_log_error = mean(abs(log(obs_pgr) - log(proj_lambda)), na.rm = TRUE),
      cor_obs_proj = cor(obs_pgr, proj_lambda, use = "complete.obs")
    )
)


# Recruitment-pulse diagnostic -------------------------------------------------

df_emerg_diag <- df_compare_final %>%
  mutate(
    obs_emerg_rate = obs_recruits_next / seedbank_start,
    pred_emerg_rate = proj_from_seedbank / seedbank_start
  ) %>%
  dplyr::select(
    year,
    obs_recruits_next,
    proj_from_seedbank,
    seedbank_start,
    obs_emerg_rate,
    pred_emerg_rate,
    obs_pgr,
    proj_lambda
  )


# =============================================================================
# PLOTS
# =============================================================================

# Shared plot settings ---------------------------------------------------------

fire_cols <- c("No fire" = "black", "Fire" = "red")
fire_shapes <- c("No fire" = 16, "Fire" = 17)
fire_lines <- c("No fire" = "solid", "Fire" = "dashed")


# Survival plot ----------------------------------------------------------------
# Raw points are not binned.
# Each year has two panels: non-flowering and flowering.
# Fire points are red; fire predictions are dashed.

surv_bin_yrs <- split(df_su, df_su$year)
v_surv <- names(surv_bin_yrs)

surv_yr_plots <- function(i) {
  
  yr_i <- v_surv[i]
  
  surv_temp <- as.data.frame(surv_bin_yrs[[i]]) %>%
    mutate(
      flower_status = if_else(flower == 1, "Flowering", "Non-flowering"),
      fire = factor(fire, levels = levels(df_su$fire)),
      year = factor(year, levels = levels(df_su$year))
    )
  
  x_temp <- seq(
    min(surv_temp$logsize_t0, na.rm = TRUE),
    max(surv_temp$logsize_t0, na.rm = TRUE),
    length.out = 100
  )
  
  pred_dat <- expand_grid(
    logsize_t0 = x_temp,
    flower = c(0, 1),
    fire = factor(c("No fire", "Fire"), levels = levels(df_su$fire))
  ) %>%
    mutate(
      flower_status = if_else(flower == 1, "Flowering", "Non-flowering")
    )
  
  pred_dat <- as.data.frame(pred_dat)
  pred_dat$year <- factor(yr_i, levels = levels(df_su$year))
  
  pred_dat$surv <- predict(
    mod_su,
    newdata = pred_dat,
    type = "response",
    re.form = NULL,
    allow.new.levels = TRUE
  )
  
  ggplot() +
    geom_jitter(
      data = surv_temp,
      aes(
        x = logsize_t0,
        y = survives,
        color = fire,
        shape = fire
      ),
      width = 0,
      height = 0.055,
      size = 0.9,
      alpha = 0.35
    ) +
    geom_line(
      data = pred_dat,
      aes(
        x = logsize_t0,
        y = surv,
        color = fire,
        linetype = fire
      ),
      linewidth = 1
    ) +
    facet_wrap(~ flower_status, nrow = 1) +
    scale_color_manual(values = fire_cols) +
    scale_shape_manual(values = fire_shapes) +
    scale_linetype_manual(values = fire_lines) +
    coord_cartesian(ylim = c(-0.05, 1.05)) +
    labs(
      title = yr_i,
      x = NULL,
      y = NULL,
      color = "Fire treatment",
      shape = "Fire treatment",
      linetype = "Fire treatment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

surv_yrs <- lapply(seq_along(surv_bin_yrs), surv_yr_plots)

g_surv_years <- wrap_plots(surv_yrs, ncol = 5, guides = "collect") +
  plot_annotation(
    title = "Year-specific survival model with flowering cost",
    subtitle = "Each year has two panels: non-flowering and flowering. Points are raw 0/1 survival observations; lines are fitted GLMM predictions.",
    caption = "X-axis: log(size) at t. Y-axis: survival probability at t + 1."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10)
  )

print(g_surv_years)


# Growth plot ------------------------------------------------------------------

grow_yr_plots <- function(i) {
  
  yr_i <- sort(unique(df_gr$year))[i]
  
  grow_temp <- df_gr %>%
    filter(year == yr_i) %>%
    mutate(
      fire = factor(fire, levels = levels(df_gr$fire)),
      year = factor(year, levels = levels(df_gr$year))
    )
  
  x_temp <- seq(
    min(grow_temp$logsize_t0, na.rm = TRUE),
    max(grow_temp$logsize_t0, na.rm = TRUE),
    length.out = 100
  )
  
  pred_dat <- expand_grid(
    logsize_t0 = x_temp,
    fire = factor(c("No fire", "Fire"), levels = levels(df_gr$fire))
  )
  
  pred_dat <- as.data.frame(pred_dat)
  pred_dat$year <- factor(yr_i, levels = levels(df_gr$year))
  
  pred_dat$pred_logsize_t1 <- predict(
    mod_gr,
    newdata = pred_dat,
    re.form = NULL,
    allow.new.levels = TRUE
  )
  
  ggplot() +
    geom_point(
      data = grow_temp,
      aes(
        x = logsize_t0,
        y = logsize_t1,
        color = fire,
        shape = fire
      ),
      size = 0.9,
      alpha = 0.35
    ) +
    geom_line(
      data = pred_dat,
      aes(
        x = logsize_t0,
        y = pred_logsize_t1,
        color = fire,
        linetype = fire
      ),
      linewidth = 1
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "blue",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    scale_color_manual(values = fire_cols) +
    scale_shape_manual(values = fire_shapes) +
    scale_linetype_manual(values = fire_lines) +
    labs(
      title = as.character(yr_i),
      x = NULL,
      y = NULL,
      color = "Fire treatment",
      shape = "Fire treatment",
      linetype = "Fire treatment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

grow_yrs <- lapply(seq_along(sort(unique(df_gr$year))), grow_yr_plots)

g_grow_years <- wrap_plots(grow_yrs, ncol = 5, guides = "collect") +
  plot_annotation(
    title = "Year-specific growth model",
    subtitle = "Points are raw observed transitions. Lines are fitted mixed-model predictions. Blue dashed line is no growth.",
    caption = "X-axis: log(size) at t. Y-axis: log(size) at t + 1."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10)
  )

print(g_grow_years)


# Flowering probability plot ---------------------------------------------------

fl_bin_yrs <- split(df_fl, df_fl$year)
v_fl <- names(fl_bin_yrs)

flower_yr_plots <- function(i) {
  
  yr_i <- v_fl[i]
  
  fl_temp <- as.data.frame(fl_bin_yrs[[i]]) %>%
    mutate(
      fire = factor(fire, levels = levels(df_fl$fire)),
      year = factor(year, levels = levels(df_fl$year))
    )
  
  x_temp <- seq(
    min(fl_temp$logsize_t0, na.rm = TRUE),
    max(fl_temp$logsize_t0, na.rm = TRUE),
    length.out = 100
  )
  
  pred_dat <- expand_grid(
    logsize_t0 = x_temp,
    fire = factor(c("No fire", "Fire"), levels = levels(df_fl$fire))
  )
  
  pred_dat <- as.data.frame(pred_dat)
  
  pred_dat$pred_flower <- predict(
    mod_fl,
    newdata = pred_dat,
    type = "response"
  )
  
  ggplot() +
    geom_jitter(
      data = fl_temp,
      aes(
        x = logsize_t0,
        y = flower,
        color = fire,
        shape = fire
      ),
      width = 0,
      height = 0.055,
      size = 0.9,
      alpha = 0.35
    ) +
    geom_line(
      data = pred_dat,
      aes(
        x = logsize_t0,
        y = pred_flower,
        color = fire,
        linetype = fire
      ),
      linewidth = 1
    ) +
    scale_color_manual(values = fire_cols) +
    scale_shape_manual(values = fire_shapes) +
    scale_linetype_manual(values = fire_lines) +
    coord_cartesian(ylim = c(-0.05, 1.05)) +
    labs(
      title = as.character(yr_i),
      x = NULL,
      y = NULL,
      color = "Fire treatment",
      shape = "Fire treatment",
      linetype = "Fire treatment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

flower_yrs <- lapply(seq_along(fl_bin_yrs), flower_yr_plots)

g_flower_years <- wrap_plots(flower_yrs, ncol = 5, guides = "collect") +
  plot_annotation(
    title = "Year-specific flowering probability model",
    subtitle = "Points are raw 0/1 flowering observations. Lines are fitted binomial-model predictions.",
    caption = "X-axis: log(size) at t. Y-axis: flowering probability."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10)
  )

print(g_flower_years)


# Flower-head number plot ------------------------------------------------------

heads_bin_yrs <- split(df_heads, df_heads$year)
v_heads <- names(heads_bin_yrs)

heads_yr_plots <- function(i) {
  
  yr_i <- v_heads[i]
  
  heads_temp <- as.data.frame(heads_bin_yrs[[i]]) %>%
    mutate(
      fire = factor(fire, levels = levels(df_heads$fire)),
      year = factor(year, levels = levels(df_heads$year))
    )
  
  x_temp <- seq(
    min(heads_temp$logsize_t0, na.rm = TRUE),
    max(heads_temp$logsize_t0, na.rm = TRUE),
    length.out = 100
  )
  
  pred_dat <- expand_grid(
    logsize_t0 = x_temp,
    fire = factor(c("No fire", "Fire"), levels = levels(df_heads$fire))
  )
  
  pred_dat <- as.data.frame(pred_dat)
  
  pred_dat$pred_heads <- predict(
    mod_heads,
    newdata = pred_dat,
    type = "response"
  )
  
  ggplot() +
    geom_point(
      data = heads_temp,
      aes(
        x = logsize_t0,
        y = fl_nr,
        color = fire,
        shape = fire
      ),
      size = 1,
      alpha = 0.45
    ) +
    geom_line(
      data = pred_dat,
      aes(
        x = logsize_t0,
        y = pred_heads,
        color = fire,
        linetype = fire
      ),
      linewidth = 1
    ) +
    scale_color_manual(values = fire_cols) +
    scale_shape_manual(values = fire_shapes) +
    scale_linetype_manual(values = fire_lines) +
    labs(
      title = as.character(yr_i),
      x = NULL,
      y = NULL,
      color = "Fire treatment",
      shape = "Fire treatment",
      linetype = "Fire treatment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

heads_yrs <- lapply(seq_along(heads_bin_yrs), heads_yr_plots)

g_heads_years <- wrap_plots(heads_yrs, ncol = 5, guides = "collect") +
  plot_annotation(
    title = "Year-specific flower-head number model",
    subtitle = "Conditional on flowering plants only. Points are raw head counts. Lines are fitted negative-binomial predictions.",
    caption = "X-axis: log(size) at t. Y-axis: flower-head number."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10)
  )

print(g_heads_years)


# Final plot: observed PGR vs modeled lambda ----------------------------------
# This resembles the old final plot, but uses the matched-site validation table.
# The projected lambda starts from the observed aboveground size distribution
# and reconstructed effective seed bank.

df_plot_final <- df_compare_final %>%
  mutate(
    fire = if_else(fire_prop > 0, "Fire", "No fire"),
    fire = factor(fire, levels = c("No fire", "Fire"))
  ) %>%
  dplyr::select(year, obs_pgr, asym_lambda, proj_lambda, fire) %>%
  pivot_longer(
    cols = c(asym_lambda, proj_lambda),
    names_to = "lambda_type",
    values_to = "lambda"
  ) %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      asym_lambda = "Asymptotic lambda",
      proj_lambda = "Projected lambda from observed size distribution"
    )
  )

g_mod_vs_obs <- ggplot(
  df_plot_final,
  aes(x = lambda, y = obs_pgr, color = fire)
) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(~ lambda_type, scales = "free_x") +
  scale_color_manual(values = fire_cols) +
  labs(
    title = "Observed population growth vs modeled lambda",
    subtitle = "Matched-site comparison; projected lambda uses observed size distribution and reconstructed effective seed bank",
    x = expression("Modeled " * lambda),
    y = "Observed matched-site population growth rate",
    color = "Fire"
  ) +
  theme_classic()

print(g_mod_vs_obs)