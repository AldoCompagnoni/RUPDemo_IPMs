# Dormancy sensitivity - observed population growth and mean IPM
# Archbold - Eriogonum longifolium
#
# Sensitivity analysis for 0-4 years of unresolved terminal dormancy.
#
# The dormancy assumption is applied only where it matters:
#   1. Observed population growth: terminal annual transitions are removed when
#      there is not enough follow-up to distinguish death from dormancy.
#   2. Mean IPM: terminal deaths are censored from the survival model under the
#      same rule, and the mean IPM is rebuilt for every cutoff.
#
# Growth, flowering, scape production, recruitment, dormancy entry and known
# reactivation observations are not discarded simply because they occur near
# the end of the study. Their observed values do not depend on deciding whether
# a disappearing individual died or remained dormant.
#
# Cutoff interpretation:
#   0 = no unresolved dormancy allowance
#   1 = allow 1 terminal year of unresolved dormancy
#   2 = allow 2 terminal years of unresolved dormancy
#   3 = current/default assumption
#   4 = allow 4 terminal years of unresolved dormancy


# Setting the stage ------------------------------------------------------------
set.seed(100)
options(stringsAsFactors = FALSE)


# Packages ---------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(MASS, tidyverse, patchwork, bbmle, lme4)


# Specification ----------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Eriogonum longifolium')

v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_ggp_suffix <- paste(
  tools::toTitleCase(v_head), '-', v_species)

# Same model-selection settings as the mean IPM.
v_mod_set_su <- c(2)
v_mod_set_gr <- c()
v_mod_set_do <- c()
v_mod_set_fl <- c()
v_mod_set_fl_n <- c()

dormancy_cutoffs <- 0:20
reference_cutoff <- 3


# Directory --------------------------------------------------------------------
dir_pub <- file.path(paste0(v_head))
dir_data <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)


# Raw data ---------------------------------------------------------------------
df_og <- read_csv(
  file.path(dir_data, "eriogonum_longifolium_data.csv"),
  col_types = cols(comment = col_character())) %>%
  janitor::clean_names() %>%
  mutate(
    year = as.numeric(str_sub(date, 1, 4)),
    month = as.numeric(str_sub(date, 6, 7)),
    id = str_c(site, pop, qu, plant, sep = "_"),
    record_type = case_when(
      !is.na(s) ~ "demography",
      is.na(s) & !is.na(burn) ~ "burn",
      TRUE ~ "other"))

df_demog <- df_og %>%
  filter(record_type == "demography")

df_demog_june <- df_demog %>%
  filter(month == 6)


# Clean annual demographic history --------------------------------------------
# Absences bracketed by known living observations are confirmed dormancy.
# These remain dormant under every sensitivity scenario.

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


# Unresolved absence spells ----------------------------------------------------
df_absence_spell <- df_annual %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    absent_state = state_clean %in% c("absent", "previous_absent"),
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


# Base annual transitions ------------------------------------------------------
df_transition_base <- df_annual %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    year_t1 = lead(year),
    annual_transition = coalesce(year_t1 == year + 1, FALSE),
    stage_t1 = if_else(
      annual_transition, lead(stage_clean), NA_real_),
    state_t1 = if_else(
      annual_transition, lead(state_clean), NA_character_),
    size_t0 = if_else(state_clean == "active", dia, NA_real_),
    size_t1 = if_else(
      annual_transition & state_clean == "active" &
        lead(state_clean) == "active",
      lead(dia), NA_real_),
    size_reactivate_t1 = if_else(
      annual_transition & state_clean == "dormant" &
        lead(state_clean) == "active",
      lead(dia), NA_real_)) %>%
  ungroup() %>%
  left_join(df_absence_t1, by = c("id", "year_t1"))


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
  left_join(df_quad_start, by = c("site", "pop", "qu")) %>%
  mutate(
    baseline = recruit_year == first_quad_year,
    recruit_type = case_when(
      s == 5 ~ "seedling",
      s == 3 ~ "new_adult"))

df_quad_monitor <- df_annual %>%
  group_by(site, pop, qu, year) %>%
  summarise(
    monitored = any(state_clean != "discontinued", na.rm = TRUE),
    .groups = "drop") %>%
  filter(monitored) %>%
  select(-monitored)

df_recruit_followup <- df_quad_monitor %>%
  left_join(
    df_quad_monitor %>%
      transmute(
        site, pop, qu, year = year - 1,
        monitored_t1 = TRUE),
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_quad_monitor %>%
      transmute(
        site, pop, qu, year = year - 2,
        monitored_t2 = TRUE),
    by = c("site", "pop", "qu", "year")) %>%
  mutate(
    monitored_t1 = replace_na(monitored_t1, FALSE),
    monitored_t2 = replace_na(monitored_t2, FALSE))

df_recruit_count_simple <- df_recruit_first %>%
  filter(!baseline) %>%
  mutate(
    fecundity_year_simple = case_when(
      recruit_type == "seedling" ~ recruit_year - 1,
      recruit_type == "new_adult" ~ recruit_year - 2)) %>%
  semi_join(
    df_quad_monitor,
    by = c(
      "site", "pop", "qu",
      "fecundity_year_simple" = "year")) %>%
  count(
    site, pop, qu, fecundity_year_simple, recruit_type,
    name = "nr_recruits") %>%
  pivot_wider(
    names_from = recruit_type,
    values_from = nr_recruits,
    values_fill = 0) %>%
  rename(
    observed_seedling_simple = seedling,
    new_adult_simple = new_adult)

df_repro_quad <- df_quad_monitor %>%
  left_join(
    df_annual %>%
      filter(state_clean != "discontinued") %>%
      group_by(site, pop, qu, year) %>%
      summarise(
        nr_active = sum(state_clean == "active", na.rm = TRUE),
        nr_flower_obs = sum(
          state_clean == "active" & !is.na(scape), na.rm = TRUE),
        scape_sum = sum(
          if_else(state_clean == "active", scape, NA_real_),
          na.rm = TRUE), .groups = "drop") %>%
      mutate(
        nr_scapes = case_when(
          nr_active == 0 ~ 0,
          nr_flower_obs == 0 ~ NA_real_,
          TRUE ~ scape_sum)) %>%
      select(site, pop, qu, year, nr_active, nr_scapes),
    by = c("site", "pop", "qu", "year"))

df_recruit_quad <- df_recruit_followup %>%
  left_join(
    df_recruit_count_simple,
    by = c(
      "site", "pop", "qu",
      "year" = "fecundity_year_simple")) %>%
  left_join(
    df_repro_quad,
    by = c("site", "pop", "qu", "year")) %>%
  mutate(
    observed_seedling_simple = case_when(
      monitored_t1 ~ replace_na(observed_seedling_simple, 0L),
      TRUE ~ NA_integer_),
    new_adult_simple = case_when(
      monitored_t2 ~ replace_na(new_adult_simple, 0L),
      TRUE ~ NA_integer_),
    recruitment_complete = monitored_t1 & monitored_t2,
    recruits_simple = case_when(
      recruitment_complete ~
        observed_seedling_simple + new_adult_simple,
      TRUE ~ NA_integer_))

df_recruit <- df_recruit_quad %>%
  transmute(
    row_type = "recruitment",
    site, pop, qu, year,
    recruitment_complete,
    recruits_simple,
    nr_scapes)

df_recruit_ind <- df_recruit_first %>%
  filter(!baseline) %>%
  transmute(
    id,
    year = recruit_year,
    recruit_type)


# Scenario-specific workdata ---------------------------------------------------
make_scenario_data <- function(dormancy_cutoff) {
  
  df_transition <- df_transition_base %>%
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
        TRUE ~ NA_real_),
      
      enter_dormancy = case_when(
        state_clean != "active" | survives != 1 ~ NA_real_,
        state_t1 == "dormant" ~ 1,
        state_t1 == "active" ~ 0,
        TRUE ~ NA_real_),
      
      reactivate = case_when(
        state_clean != "dormant" ~ NA_real_,
        state_t1 == "active" ~ 1,
        state_t1 == "dormant" ~ 0,
        TRUE ~ NA_real_),
      
      flower = case_when(
        state_clean != "active" ~ NA_real_,
        !is.na(scape) & scape > 0 ~ 1,
        !is.na(scape) & scape == 0 ~ 0,
        TRUE ~ NA_real_),
      
      fl_nr = if_else(
        state_clean == "active", scape, NA_real_),
      logsize_t0 = log(size_t0),
      logsize_t1 = log(size_t1),
      logsize_t0_2 = logsize_t0^2,
      logsize_t0_3 = logsize_t0^3)
  
  df_ind <- df_transition %>%
    left_join(df_recruit_ind, by = c("id", "year")) %>%
    transmute(
      row_type = "individual",
      site, pop, qu, id, year,
      state = state_clean,
      state_t1,
      stage = stage_clean,
      stage_t1,
      survives,
      enter_dormancy,
      recruit_type,
      size_reactivate_t1,
      reactivate,
      size_t0,
      size_t1,
      logsize_t0,
      logsize_t1,
      logsize_t0_2,
      logsize_t0_3,
      flower,
      fl_nr)
  
  bind_rows(df_ind, df_recruit) %>%
    mutate(
      row_type = as.character(row_type),
      year = as.numeric(year))
}


# Model helpers ----------------------------------------------------------------
poly_formulas <- function(response) {
  list(
    as.formula(paste(response, '~ 1')),
    as.formula(paste(response, '~ logsize_t0')),
    as.formula(paste(
      response, '~ logsize_t0 + logsize_t0_2')),
    as.formula(paste(
      response,
      '~ logsize_t0 + logsize_t0_2 + logsize_t0_3')))
}

get_daicc <- function(mods) {
  aicc <- sapply(mods, function(mod) {
    k <- attr(logLik(mod), "df")
    n <- nobs(mod)
    AIC(mod) + 2 * k * (k + 1) / (n - k - 1)
  })
  
  aicc - min(aicc)
}

pick_model <- function(mods, manual = c()) {
  d_aicc <- get_daicc(mods)
  index <- if (length(manual) == 0) {
    which.min(d_aicc)
  } else {
    manual[1] + 1
  }

  list(
    model = mods[[index]],
    index = index - 1,
    d_aicc = d_aicc)
}

fit_poly_glm <- function(data, response, manual = c()) {
  forms <- poly_formulas(response)
  mods <- lapply(
    forms,
    function(form) glm(form, data = data, family = 'binomial'))
  pick_model(mods, manual)
}

fit_poly_lm <- function(data, response, manual = c()) {
  forms <- poly_formulas(response)
  mods <- lapply(forms, function(form) lm(form, data = data))
  pick_model(mods, manual)
}

fit_poly_nb <- function(data, response, manual = c()) {
  forms <- poly_formulas(response)
  mods <- lapply(forms, function(form) MASS::glm.nb(form, data = data))
  pick_model(mods, manual)
}

get_coef <- function(model, term) {
  values <- coef(model)
  if (term %in% names(values)) {
    unname(values[[term]])
  } else {
    NULL
  }
}


# IPM functions ----------------------------------------------------------------
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

sx <- function(x, pars, num_pars = pars$mod_su_index) {
  val <- pars$surv_b0

  for (i in seq_len(num_pars)) {
    param <- paste0('surv_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }

  inv_logit(val)
}

grow_mu <- function(x, pars, num_pars = pars$mod_gr_index) {
  val <- pars$grow_b0

  for (i in seq_len(num_pars)) {
    param <- paste0('grow_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }

  val
}

grow_sd <- function(mu, pars) {
  sqrt(pars$a * exp(pars$b * mu))
}

gxy <- function(x, y, pars) {
  mu <- grow_mu(x, pars)
  dnorm(y, mean = mu, sd = grow_sd(mu, pars))
}

dorm_x <- function(x, pars, num_pars = pars$mod_do_index) {
  val <- pars$dorm_b0

  for (i in seq_len(num_pars)) {
    param <- paste0('dorm_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }

  inv_logit(val)
}

fl_x <- function(x, pars, num_pars = pars$mod_fl_index) {
  val <- pars$fl_b0

  for (i in seq_len(num_pars)) {
    param <- paste0('fl_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }

  inv_logit(val)
}

fl_n_x <- function(x, pars, num_pars = pars$mod_fl_n_index) {
  val <- pars$fln_b0

  for (i in seq_len(num_pars)) {
    param <- paste0('fln_b', i)
    if (!is.null(pars[[param]])) {
      val <- val + pars[[param]] * x^i
    }
  }

  exp(val)
}

re_y_dist <- function(y, pars) {
  norm <- pnorm(
    pars$U, mean = pars$recr_sz, sd = pars$recr_sd) -
    pnorm(
      pars$L, mean = pars$recr_sz, sd = pars$recr_sd)

  dnorm(
    y, mean = pars$recr_sz, sd = pars$recr_sd) / norm
}

react_y_dist <- function(y, pars) {
  norm <- pnorm(
    pars$U, mean = pars$react_sz, sd = pars$react_sd) -
    pnorm(
      pars$L, mean = pars$react_sz, sd = pars$react_sd)

  dnorm(
    y, mean = pars$react_sz, sd = pars$react_sd) / norm
}

fyx <- function(y, x, pars) {
  fl_x(x, pars) *
    fl_n_x(x, pars) *
    pars$fecu_b0 *
    re_y_dist(y, pars)
}

kernel <- function(pars) {
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n

  b <- L + c(0:n) * h
  y <- 0.5 * (b[1:n] + b[2:(n + 1)])

  Smat <- sx(y, pars)
  Dmat <- dorm_x(y, pars)

  Gmat <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h

  for (i in seq_len(n)) {
    if (i <= n / 2) {
      Gmat[1, i] <- Gmat[1, i] + 1 - sum(Gmat[, i])
    } else {
      Gmat[n, i] <- Gmat[n, i] + 1 - sum(Gmat[, i])
    }
  }

  Tmat <- sweep(
    Gmat, 2, Smat * (1 - Dmat), '*')

  Dorm_row <- matrix(
    Smat * Dmat, nrow = 1)

  React_vec <- react_y_dist(y, pars) * h
  React_vec <- React_vec / sum(React_vec)

  React_col <- matrix(
    pars$p_reactivate * React_vec, ncol = 1)

  Dorm_stasis <- 1 - pars$p_reactivate

  Fmat <- outer(
    y, y,
    Vectorize(function(y, x) {
      fyx(y, x, pars)
    })) * h

  K <- rbind(
    cbind(Tmat + Fmat, React_col),
    cbind(Dorm_row, Dorm_stasis))

  K
}

lambda_ipm <- function(pars) {
  Re(eigen(
    kernel(pars), only.values = TRUE)$values[1])
}


# Reference data for scenario-independent population counts -------------------
df_reference <- make_scenario_data(reference_cutoff)

df_ind_reference <- df_reference %>%
  filter(row_type == "individual")

max_ind_year <- max(df_annual$year, na.rm = TRUE)


# Observed population counts ---------------------------------------------------
df_counts_alive <- df_ind_reference %>%
  filter(state %in% c("active", "dormant", "alive_or_dormant")) %>%
  count(site, pop, qu, year, name = "n_alive")

df_newadult_backfill <- df_ind_reference %>%
  filter(recruit_type == "new_adult") %>%
  transmute(site, pop, qu, year = year - 1) %>%
  count(site, pop, qu, year, name = "n_newadult_backfill")

df_counts_quad <- df_reference %>%
  filter(row_type == "recruitment") %>%
  distinct(site, pop, qu, year) %>%
  left_join(
    df_counts_alive,
    by = c("site", "pop", "qu", "year")) %>%
  left_join(
    df_newadult_backfill,
    by = c("site", "pop", "qu", "year")) %>%
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


# One sensitivity scenario -----------------------------------------------------
run_dormancy_scenario <- function(dormancy_cutoff) {
  
  df_i <- make_scenario_data(dormancy_cutoff)
  
  df_ind_i <- df_i %>%
    filter(row_type == "individual")
  
  df_re_i <- df_i %>%
    filter(
      row_type == "recruitment",
      recruitment_complete,
      !is.na(recruits_simple))
  
  
  # Survival -------------------------------------------------------------------
  df_su_i <- df_ind_i %>%
    filter(
      state == "active",
      !is.na(survives),
      size_t0 > 0,
      is.finite(logsize_t0))
  
  fit_su_i <- fit_poly_glm(
    df_su_i, "survives", v_mod_set_su)
  
  mod_su_i <- fit_su_i$model
  
  
  # Growth ---------------------------------------------------------------------
  df_gr_i <- df_ind_i %>%
    filter(
      state == "active",
      state_t1 == "active",
      size_t0 > 0,
      size_t1 > 0,
      is.finite(logsize_t0),
      is.finite(logsize_t1))
  
  fit_gr_i <- fit_poly_lm(
    df_gr_i, "logsize_t1", v_mod_set_gr)
  
  mod_gr_i <- fit_gr_i$model
  
  mod_gr_x_i <- fitted(mod_gr_i)
  mod_gr_y_i <- resid(mod_gr_i)^2
  
  mod_gr_var_i <- nls(
    mod_gr_y_i ~ a * exp(b * mod_gr_x_i),
    start = list(a = 1, b = 0),
    control = nls.control(
      maxiter = 1000, tol = 1e-6, warnOnly = TRUE))
  
  
  # Dormancy entry -------------------------------------------------------------
  df_do_i <- df_ind_i %>%
    filter(
      !is.na(enter_dormancy),
      size_t0 > 0,
      is.finite(logsize_t0))
  
  fit_do_i <- fit_poly_glm(
    df_do_i, "enter_dormancy", v_mod_set_do)
  
  mod_do_i <- fit_do_i$model
  
  
  # Reactivation ---------------------------------------------------------------
  df_ra_i <- df_ind_i %>%
    filter(!is.na(reactivate))
  
  mod_ra_i <- glm(
    reactivate ~ 1,
    data = df_ra_i, family = "binomial")
  
  p_reactivate_i <- predict(
    mod_ra_i,
    newdata = data.frame(x = 1),
    type = "response")[1]
  
  df_ra_size_i <- df_ind_i %>%
    filter(
      reactivate == 1,
      size_reactivate_t1 > 0) %>%
    mutate(
      logsize_reactivate = log(size_reactivate_t1))
  
  react_sz_i <- mean(
    df_ra_size_i$logsize_reactivate, na.rm = TRUE)
  
  react_sd_i <- sd(
    df_ra_size_i$logsize_reactivate, na.rm = TRUE)
  
  
  # Flowering ------------------------------------------------------------------
  df_fl_i <- df_ind_i %>%
    filter(
      state == "active",
      !is.na(flower),
      size_t0 > 0,
      is.finite(logsize_t0))
  
  fit_fl_i <- fit_poly_glm(
    df_fl_i, "flower", v_mod_set_fl)
  
  mod_fl_i <- fit_fl_i$model
  
  
  # Flower number --------------------------------------------------------------
  df_fl_cond_i <- df_fl_i %>%
    filter(
      flower == 1,
      !is.na(fl_nr),
      fl_nr > 0,
      fl_nr %% 1 == 0)
  
  fit_fln_i <- fit_poly_nb(
    df_fl_cond_i, "fl_nr", v_mod_set_fl_n)
  
  mod_fln_i <- fit_fln_i$model
  
  
  # Recruitment ----------------------------------------------------------------
  df_sc2re_i <- df_re_i %>%
    filter(
      !is.na(nr_scapes),
      !is.na(recruits_simple)) %>%
    group_by(site, pop, year) %>%
    summarise(
      total_scapes = sum(nr_scapes),
      recruit_count = sum(recruits_simple), .groups = "drop")
  
  mods_re_i <- list(
    MASS::glm.nb(recruit_count ~ 1, data = df_sc2re_i),
    MASS::glm.nb(recruit_count ~ total_scapes, data = df_sc2re_i),
    MASS::glm.nb(
      recruit_count ~ log1p(total_scapes), data = df_sc2re_i))
  
  mod_re_i <- mods_re_i[[which.min(get_daicc(mods_re_i))]]
  
  df_sc2re_i <- df_sc2re_i %>%
    mutate(
      recruits_pred = predict(
        mod_re_i, type = "response"))
  
  recr_per_scape_i <- sum(df_sc2re_i$recruits_pred) /
    sum(df_sc2re_i$total_scapes)
  
  
  # Recruit size ---------------------------------------------------------------
  df_re_size_i <- df_ind_i %>%
    filter(
      recruit_type == "seedling",
      size_t0 > 0,
      is.finite(logsize_t0))
  
  recr_sz_i <- mean(
    df_re_size_i$logsize_t0, na.rm = TRUE)
  
  recr_sd_i <- sd(
    df_re_size_i$logsize_t0, na.rm = TRUE)
  
  
  # IPM parameters -------------------------------------------------------------
  pars_i <- Filter(function(x) length(x) > 0, list(
    surv_b0 = get_coef(mod_su_i, "(Intercept)"),
    surv_b1 = get_coef(mod_su_i, "logsize_t0"),
    surv_b2 = get_coef(mod_su_i, "logsize_t0_2"),
    surv_b3 = get_coef(mod_su_i, "logsize_t0_3"),
    
    grow_b0 = get_coef(mod_gr_i, "(Intercept)"),
    grow_b1 = get_coef(mod_gr_i, "logsize_t0"),
    grow_b2 = get_coef(mod_gr_i, "logsize_t0_2"),
    grow_b3 = get_coef(mod_gr_i, "logsize_t0_3"),
    a = unname(coef(mod_gr_var_i)[["a"]]),
    b = unname(coef(mod_gr_var_i)[["b"]]),
    
    dorm_b0 = get_coef(mod_do_i, "(Intercept)"),
    dorm_b1 = get_coef(mod_do_i, "logsize_t0"),
    dorm_b2 = get_coef(mod_do_i, "logsize_t0_2"),
    dorm_b3 = get_coef(mod_do_i, "logsize_t0_3"),
    
    fl_b0 = get_coef(mod_fl_i, "(Intercept)"),
    fl_b1 = get_coef(mod_fl_i, "logsize_t0"),
    fl_b2 = get_coef(mod_fl_i, "logsize_t0_2"),
    fl_b3 = get_coef(mod_fl_i, "logsize_t0_3"),
    
    fln_b0 = get_coef(mod_fln_i, "(Intercept)"),
    fln_b1 = get_coef(mod_fln_i, "logsize_t0"),
    fln_b2 = get_coef(mod_fln_i, "logsize_t0_2"),
    fln_b3 = get_coef(mod_fln_i, "logsize_t0_3"),
    
    fecu_b0 = recr_per_scape_i,
    recr_sz = recr_sz_i,
    recr_sd = recr_sd_i,
    react_sz = react_sz_i,
    react_sd = react_sd_i,
    p_reactivate = p_reactivate_i,
    
    L = min(df_gr_i$logsize_t0),
    U = max(df_gr_i$logsize_t0),
    mat_siz = 200,
    
    mod_su_index = fit_su_i$index,
    mod_gr_index = fit_gr_i$index,
    mod_do_index = fit_do_i$index,
    mod_fl_index = fit_fl_i$index,
    mod_fl_n_index = fit_fln_i$index))
  
  lambda_ipm_i <- lambda_ipm(pars_i)
  
  
  # Observed population growth -------------------------------------------------
  last_complete_transition <- max_ind_year - dormancy_cutoff - 1
  
  df_counts_year_i <- df_counts_quad %>%
    filter(year <= last_complete_transition) %>%
    group_by(year) %>%
    summarise(
      n_quads = n(),
      n_t0 = sum(n),
      n_t1 = sum(n_t1), .groups = "drop") %>%
    filter(n_t0 > 0) %>%
    mutate(lambda_obs = n_t1 / n_t0)
  
  lambda_obs_arithmetic <- mean(
    df_counts_year_i$lambda_obs, na.rm = TRUE)
  
  lambda_obs_geometric <- exp(mean(
    log(df_counts_year_i$lambda_obs), na.rm = TRUE))
  
  
  # Scenario summary -----------------------------------------------------------
  n_absence_dead <- df_ind_i %>%
    filter(
      state == "active",
      state_t1 %in% c("absent", "previous_absent"),
      survives == 0) %>%
    nrow()
  
  n_absence_unresolved <- df_ind_i %>%
    filter(
      state == "active",
      state_t1 %in% c("absent", "previous_absent"),
      is.na(survives)) %>%
    nrow()
  
  tibble(
    dormancy_cutoff = dormancy_cutoff,
    last_complete_transition = last_complete_transition,
    n_observed_transitions = nrow(df_counts_year_i),
    n_absence_dead = n_absence_dead,
    n_absence_unresolved = n_absence_unresolved,
    n_survival = nrow(df_su_i),
    n_survival_dead = sum(df_su_i$survives == 0),
    n_survival_alive = sum(df_su_i$survives == 1),
    survival_model_degree = fit_su_i$index,
    lambda_ipm = lambda_ipm_i,
    lambda_obs_arithmetic = lambda_obs_arithmetic,
    lambda_obs_geometric = lambda_obs_geometric)
}


# Run 0-4 years ----------------------------------------------------------------
df_sensitivity <- map_dfr(
  dormancy_cutoffs, run_dormancy_scenario)


# Summary relative to the current 3-year assumption ----------------------------
ref <- df_sensitivity %>%
  filter(dormancy_cutoff == reference_cutoff) %>%
  transmute(
    ref_ipm = lambda_ipm,
    ref_obs_arithmetic = lambda_obs_arithmetic,
    ref_obs_geometric = lambda_obs_geometric)

df_sensitivity_summary <- df_sensitivity %>%
  crossing(ref) %>%
  mutate(
    delta_ipm = lambda_ipm - ref_ipm,
    delta_obs_arithmetic =
      lambda_obs_arithmetic - ref_obs_arithmetic,
    delta_obs_geometric =
      lambda_obs_geometric - ref_obs_geometric,
    pct_delta_ipm = 100 * delta_ipm / ref_ipm,
    pct_delta_obs_arithmetic =
      100 * delta_obs_arithmetic / ref_obs_arithmetic,
    pct_delta_obs_geometric =
      100 * delta_obs_geometric / ref_obs_geometric) %>%
  select(-starts_with('ref_'))

print(df_sensitivity_summary, n = Inf)


# Final sensitivity plot -------------------------------------------------------
df_sensitivity_plot <- df_sensitivity_summary %>%
  select(
    dormancy_cutoff,
    `Mean IPM` = lambda_ipm,
    `Observed arithmetic` = lambda_obs_arithmetic,
    `Observed geometric` = lambda_obs_geometric) %>%
  pivot_longer(
    -dormancy_cutoff,
    names_to = 'estimate',
    values_to = 'lambda')

fig_dormancy_sensitivity <- ggplot(
  df_sensitivity_plot,
  aes(
    x = dormancy_cutoff,
    y = lambda,
    linetype = estimate,
    shape = estimate)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_vline(
    xintercept = reference_cutoff,
    linetype = 'dotted') +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_x_continuous(breaks = dormancy_cutoffs) +
  theme_bw() +
  labs(
    title = 'Sensitivity to terminal dormancy assumption',
    subtitle = v_ggp_suffix,
    x = 'Maximum unresolved dormancy (years)',
    y = expression(lambda),
    linetype = NULL,
    shape = NULL)

fig_dormancy_sensitivity


# Optional saves ---------------------------------------------------------------
# write.csv(
#   df_sensitivity_summary,
#   file.path(
#     dir_result,
#     paste0('ab_', v_sp_abb,
#            '_dormancy_sensitivity_obs_ipm_0to4.csv')),
#   row.names = FALSE)
#
# ggsave(
#   file.path(
#     dir_result,
#     paste0('ab_', v_sp_abb,
#            '_dormancy_sensitivity_obs_ipm_0to4.png')),
#   plot = fig_dormancy_sensitivity,
#   width = 8, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# Size-data diagnostics and main findings --------------------------------------
# Diagnostic only. No diameter observation is corrected, replaced or removed.
# The purpose is to identify and summarize the unusual size observations that
# matter most for the growth model.


# Low active non-seedling measurements -----------------------------------------
df_annual <- df_annual %>%
  mutate(
    flag_low_nonseedling =
      state_clean == "active" &
      stage_clean != 1 &
      !is.na(dia) &
      dia < 1)

df_low_summary <- df_annual %>%
  count(flag_low_nonseedling)

print(df_low_summary, n = Inf)


df_low_nonseedling <- df_annual %>%
  filter(flag_low_nonseedling) %>%
  select(
    site, pop, qu, plant, id, year,
    s, stage, stage_clean,
    dia, scape, invol, nb, nr, herb, burn, comment) %>%
  arrange(year, site, pop, qu, plant)

print(df_low_nonseedling, n = Inf, width = Inf)


# Diameter heaping and measurement precision ----------------------------------
df_heaping_summary <- df_annual %>%
  filter(
    state_clean == "active",
    !is.na(dia),
    dia > 0) %>%
  summarise(
    n = n(),
    exact_integer = sum(dia %% 1 == 0),
    prop_integer = mean(dia %% 1 == 0),
    exact_half = sum((dia * 2) %% 1 == 0),
    prop_half = mean((dia * 2) %% 1 == 0))

print(df_heaping_summary)


df_heap_year <- df_annual %>%
  filter(
    state_clean == "active",
    !is.na(dia),
    dia > 0) %>%
  group_by(year) %>%
  summarise(
    n = n(),
    prop_integer = mean(dia %% 1 == 0),
    prop_half = mean((dia * 2) %% 1 == 0),
    n_unique = n_distinct(dia),
    .groups = "drop")

print(df_heap_year, n = Inf)


fig_heaping_year <- df_heap_year %>%
  select(
    year,
    `Exact integer` = prop_integer,
    `Multiple of 0.5` = prop_half) %>%
  pivot_longer(
    -year,
    names_to = "precision",
    values_to = "proportion") %>%
  ggplot(
    aes(
      x = year,
      y = proportion,
      linetype = precision,
      shape = precision)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "Diameter measurement precision through time",
    subtitle = v_ggp_suffix,
    x = "Year",
    y = "Proportion of active diameter records",
    linetype = NULL,
    shape = NULL)

fig_heaping_year


# Three-year size valleys -------------------------------------------------------
# A positive valley_depth means that the focal year's diameter is below both
# neighboring years on the log scale. Larger values indicate a deeper valley.

df_three_year <- df_annual %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    year_prev = lag(year),
    year_next = lead(year),
    state_prev = lag(state_clean),
    state_next = lead(state_clean),
    dia_prev = lag(dia),
    dia_next = lead(dia),
    stage_prev = lag(stage_clean),
    stage_next = lead(stage_clean),
    scape_prev = lag(scape),
    scape_next = lead(scape),
    herb_prev = lag(herb),
    herb_next = lead(herb),
    comment_prev = lag(comment),
    comment_next = lead(comment)) %>%
  ungroup() %>%
  filter(
    year_prev == year - 1,
    year_next == year + 1,
    state_prev == "active",
    state_clean == "active",
    state_next == "active",
    dia_prev > 0,
    dia > 0,
    dia_next > 0) %>%
  mutate(
    logdia_prev = log(dia_prev),
    logdia = log(dia),
    logdia_next = log(dia_next),
    valley_depth =
      pmin(logdia_prev, logdia_next) - logdia,
    local_valley =
      dia < dia_prev & dia < dia_next,
    below_log0 = logdia < 0,
    below0_between_above0 =
      dia < 1 & dia_prev >= 1 & dia_next >= 1) %>%
  left_join(
    df_heap_year %>%
      select(
        year,
        year_prop_integer = prop_integer,
        year_prop_half = prop_half),
    by = "year")


# Strict big-small-big observations -------------------------------------------
df_suspicious <- df_three_year %>%
  filter(below0_between_above0) %>%
  arrange(desc(valley_depth))

print(
  df_suspicious %>%
    select(
      site, pop, qu, plant, id, year,
      stage_prev, stage_clean, stage_next,
      scape_prev, scape, scape_next,
      herb_prev, herb, herb_next,
      dia_prev, dia, dia_next,
      valley_depth,
      comment_prev, comment, comment_next),
  n = Inf,
  width = Inf)


# Strongest valleys in the complete three-year dataset ------------------------
df_top_valleys <- df_three_year %>%
  filter(local_valley) %>%
  arrange(desc(valley_depth)) %>%
  select(
    site, pop, qu, plant, id, year,
    stage_prev, stage_clean, stage_next,
    scape_prev, scape, scape_next,
    herb,
    dia_prev, dia, dia_next,
    valley_depth) %>%
  slice_head(n = 20)

print(df_top_valleys, n = Inf, width = Inf)


# Raw variables for the strict suspicious observations ------------------------
df_suspicious_full <- df_annual %>%
  inner_join(
    df_suspicious %>%
      select(id, year, valley_depth),
    by = c("id", "year")) %>%
  arrange(desc(valley_depth))

print(df_suspicious_full, n = Inf, width = Inf)


# Original records in the year before, during and after each suspicious event --
df_suspicious_year <- df_suspicious %>%
  select(id, focal_year = year, valley_depth)


df_suspicious_raw <- df_og %>%
  inner_join(
    df_suspicious_year,
    by = "id") %>%
  filter(
    year >= focal_year - 1,
    year <= focal_year + 1) %>%
  arrange(id, focal_year, year, month)

print(df_suspicious_raw, n = Inf, width = Inf)


# Plot the seven strict big-small-big histories --------------------------------
df_suspicious_long <- df_suspicious %>%
  select(id, year, dia_prev, dia, dia_next) %>%
  pivot_longer(
    cols = c(dia_prev, dia, dia_next),
    names_to = "position",
    values_to = "diameter") %>%
  mutate(
    census_year = case_when(
      position == "dia_prev" ~ year - 1,
      position == "dia" ~ year,
      position == "dia_next" ~ year + 1))

fig_suspicious_trajectories <- ggplot(
  df_suspicious_long,
  aes(
    x = census_year,
    y = diameter,
    group = id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ id, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Flagged big-small-big diameter histories",
    subtitle = v_ggp_suffix,
    x = "Year",
    y = "Diameter")

fig_suspicious_trajectories


# Broad association screen for valley depth -----------------------------------
# This is descriptive only. It is not used to infer causation.
# Herbivory is handled separately below because the codes should not be treated
# as a continuous severity scale without explicit justification.

v_valley_cor <- intersect(
  c(
    "scape",
    "invol",
    "nb",
    "nr",
    "year_prop_integer",
    "year_prop_half"),
  names(df_three_year))


df_valley_cor <- map_dfr(
  v_valley_cor,
  function(v) {
    x <- df_three_year[[v]]
    y <- df_three_year$valley_depth

    keep <- is.finite(x) & is.finite(y)

    if (sum(keep) < 10 ||
        n_distinct(x[keep]) < 2) {
      return(tibble(
        variable = v,
        n = sum(keep),
        rho = NA_real_))
    }

    tibble(
      variable = v,
      n = sum(keep),
      rho = cor(
        x[keep],
        y[keep],
        method = "spearman"))
  })

print(
  df_valley_cor %>%
    arrange(desc(abs(rho))),
  n = Inf)


# Stage sequences containing strict suspicious valleys -------------------------
df_stage_valleys <- df_three_year %>%
  mutate(suspicious = below0_between_above0) %>%
  group_by(
    stage_prev,
    stage_clean,
    stage_next) %>%
  summarise(
    n = n(),
    n_suspicious = sum(suspicious),
    prop_suspicious = mean(suspicious),
    median_valley = median(valley_depth),
    max_valley = max(valley_depth),
    .groups = "drop") %>%
  filter(n_suspicious > 0) %>%
  arrange(desc(n_suspicious), desc(max_valley))

print(df_stage_valleys, n = Inf)


# Herbivory --------------------------------------------------------------------
# herb == 1 is recorded herbivory. Codes 2 and 3 are retained and displayed,
# but are not assumed here to represent increasing severity.

df_herb_codes <- df_three_year %>%
  count(herb, .drop = FALSE)

print(df_herb_codes, n = Inf)


df_herb_test <- df_three_year %>%
  filter(herb %in% c(0, 1)) %>%
  mutate(
    shrink_in = log(dia) - log(dia_prev),
    rebound_out = log(dia_next) - log(dia))


df_herb_summary <- df_herb_test %>%
  group_by(herb) %>%
  summarise(
    n = n(),
    median_dia_prev = median(dia_prev),
    median_dia = median(dia),
    median_dia_next = median(dia_next),
    mean_shrink_in = mean(shrink_in),
    median_shrink_in = median(shrink_in),
    mean_rebound_out = mean(rebound_out),
    median_rebound_out = median(rebound_out),
    mean_valley_depth = mean(valley_depth),
    median_valley_depth = median(valley_depth),
    prop_local_valley = mean(local_valley),
    .groups = "drop")

print(df_herb_summary, width = Inf)


mod_herb_shrink <- lme4::lmer(
  log(dia) ~ log(dia_prev) + herb +
    (1 | year) +
    (1 | id),
  data = df_herb_test)

mod_herb_rebound <- lme4::lmer(
  log(dia_next) ~ log(dia) + herb +
    (1 | year) +
    (1 | id),
  data = df_herb_test)

summary(mod_herb_shrink)
summary(mod_herb_rebound)


# Suspicious observations with herbivory information ---------------------------
df_suspicious_herb <- df_suspicious %>%
  select(
    site, pop, qu, plant, id, year,
    dia_prev, dia, dia_next,
    stage_clean,
    scape,
    invol,
    nb,
    nr,
    herb,
    valley_depth)

print(df_suspicious_herb, n = Inf, width = Inf)


# Compact summary of the main findings ----------------------------------------
herb_shrink_coef <- unname(lme4::fixef(mod_herb_shrink)[["herb"]])
herb_rebound_coef <- unname(lme4::fixef(mod_herb_rebound)[["herb"]])


df_main_findings <- tibble(
  finding = c(
    "Positive active diameter observations",
    "Exact integer diameter observations",
    "Exact half-centimeter diameter observations",
    "Flagged active non-seedling diameters below 1",
    "Strict >=1 -> <1 -> >=1 three-year valleys",
    "Strict valleys with herbivory recorded as 1",
    "Strict valleys with herbivory missing",
    "Herb effect in shrinkage mixed model",
    "Herb effect in rebound mixed model"),
  value = c(
    df_heaping_summary$n,
    df_heaping_summary$prop_integer,
    df_heaping_summary$prop_half,
    sum(df_annual$flag_low_nonseedling, na.rm = TRUE),
    nrow(df_suspicious),
    sum(df_suspicious$herb == 1, na.rm = TRUE),
    sum(is.na(df_suspicious$herb)),
    herb_shrink_coef,
    herb_rebound_coef))

print(df_main_findings, n = Inf)


# Key interpretation -----------------------------------------------------------
# 1. Diameter heaping is substantial, but strongly changes through time and is
#    therefore best treated as a measurement-resolution feature rather than a
#    unique problem at diameter = 1.
# 2. The main growth-data concern is the small set of active non-seedling
#    observations below 1, especially the strict big-small-big histories.
# 3. These observations are flagged only; no value is corrected or discarded.
# 4. The strongest strict valley includes herb == 1 in the focal year, providing
#    a biologically relevant observation for that individual. However, herbivory
#    is missing for most strict valleys and the 0-vs-1 mixed models do not show a
#    strong general herbivory effect on shrinkage or rebound.
# 5. The flagged observations should therefore remain identifiable for later
#    sensitivity analysis of the growth model rather than being automatically
#    altered.




# Flowering and unusually small individuals ------------------------------------
# Investigation summary:
#
# - Flowering at very small size is rare.
# - Only 1 flowering observation occurs below 1 cm.
# - Only 9 flowering observations occur below 3 cm.
# - Flowering in general is NOT associated with a temporary reduction in
#   diameter. Flowering plants are, on average, slightly larger than expected
#   from their surrounding-year sizes.
# - However, several of the strongest flowering-year size valleys are extremely
#   unusual, including:
#
#     1_7_8_538: 24.0 -> 0.5 -> 19.9, 3 scapes, herbivory recorded
#     1_5_4_912: 24.0 -> 1.0 -> 23.0, 1 scape, no herbivory recorded
#     1_4_1_947: 16.8 -> 2.0 -> 13.5, 1 scape
#     1_5_1_80:  12.0 -> 2.0 -> 13.0, 1 scape
#
# These observations are retained unchanged and flagged for reference.


# Three-year flowering context -------------------------------------------------
df_fl_valley <- df_annual %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    year_prev = lag(year),
    year_next = lead(year),
    dia_prev = lag(dia),
    dia_next = lead(dia),
    stage_prev = lag(stage_clean),
    stage_next = lead(stage_clean),
    scape_prev = lag(scape),
    scape_next = lead(scape)) %>%
  ungroup() %>%
  filter(
    year_prev == year - 1,
    year_next == year + 1,
    state_clean == "active",
    dia > 0,
    dia_prev > 0,
    dia_next > 0,
    !is.na(scape)) %>%
  mutate(
    flower = scape > 0,
    log_change_in = log(dia / dia_prev),
    log_change_out = log(dia_next / dia),
    expected_dia = sqrt(dia_prev * dia_next),
    deviation_middle = log(dia / expected_dia),
    local_valley = dia < dia_prev & dia < dia_next)


# Main flowering-size summary --------------------------------------------------
df_flowering_size_summary <- df_annual %>%
  filter(
    state_clean == "active",
    !is.na(dia),
    dia > 0,
    !is.na(scape),
    scape > 0) %>%
  summarise(
    n_flowering = n(),
    n_flowering_lt1 = sum(dia < 1),
    n_flowering_lt3 = sum(dia < 3),
    n_flowering_le5 = sum(dia <= 5),
    min_flowering_dia = min(dia))

df_flowering_size_summary


# Flowering vs non-flowering size valleys --------------------------------------
df_flowering_valley_summary <- df_fl_valley %>%
  group_by(flower) %>%
  summarise(
    n = n(),
    median_dia_prev = median(dia_prev),
    median_dia = median(dia),
    median_dia_next = median(dia_next),
    mean_change_in = mean(log_change_in),
    median_change_in = median(log_change_in),
    mean_change_out = mean(log_change_out),
    median_change_out = median(log_change_out),
    mean_deviation_middle = mean(deviation_middle),
    median_deviation_middle = median(deviation_middle),
    prop_local_valley = mean(local_valley),
    .groups = "drop")

df_flowering_valley_summary


# Models confirming the general flowering pattern -----------------------------
mod_fl_dia <- lmer(
  log(dia) ~ log(dia_prev) + flower +
    (1 | year) +
    (1 | id),
  data = df_fl_valley)

summary(mod_fl_dia)

# Individual variance for this response is effectively zero, so year alone is
# sufficient for the middle-year deviation model.
mod_fl_valley <- lmer(
  deviation_middle ~ flower +
    (1 | year),
  data = df_fl_valley)

summary(mod_fl_valley)


# Small flowering observations -------------------------------------------------
df_small_flowering <- df_annual %>%
  filter(
    state_clean == "active",
    !is.na(dia),
    dia < 3,
    !is.na(scape),
    scape > 0) %>%
  arrange(dia, year) %>%
  select(
    site, pop, qu, plant, id, year,
    stage_clean,
    dia, scape,
    invol, nb, nr, herb,
    comment)

df_small_flowering %>%
  print(n = Inf, width = Inf)


# Strongest flowering-year size valleys ----------------------------------------
df_flowering_extreme <- df_fl_valley %>%
  filter(flower) %>%
  arrange(deviation_middle) %>%
  select(
    site, pop, qu, plant, id, year,
    stage_prev, stage_clean, stage_next,
    scape_prev, scape, scape_next,
    dia_prev, dia, dia_next,
    expected_dia,
    deviation_middle,
    herb, invol, nb, nr)

df_flowering_extreme %>%
  slice_head(n = 20) %>%
  print(n = Inf, width = Inf)


# Flags in annual data ----------------------------------------------------------
# These flags identify observations for sensitivity checks only.
# No diameter values are altered.

df_annual <- df_annual %>%
  mutate(
    flag_small_flowering =
      state_clean == "active" &
      !is.na(dia) &
      dia < 3 &
      !is.na(scape) &
      scape > 0,
    
    flag_flowering_below1 =
      state_clean == "active" &
      !is.na(dia) &
      dia < 1 &
      !is.na(scape) &
      scape > 0)


# Flag the four most extreme flowering-year valleys ----------------------------
# These are selected by their observed rank rather than by changing or assuming
# anything about the recorded measurements.

ids_extreme_flowering <- df_flowering_extreme %>%
  slice_head(n = 4) %>%
  select(id, year) %>%
  mutate(flag_extreme_flowering_valley = TRUE)

df_annual <- df_annual %>%
  left_join(
    ids_extreme_flowering,
    by = c("id", "year")) %>%
  mutate(
    flag_extreme_flowering_valley =
      replace_na(flag_extreme_flowering_valley, FALSE))


# Final flag counts -------------------------------------------------------------
df_annual %>%
  summarise(
    n_small_flowering = sum(flag_small_flowering),
    n_flowering_below1 = sum(flag_flowering_below1),
    n_extreme_flowering_valley =
      sum(flag_extreme_flowering_valley))


# Final flagged observations ---------------------------------------------------
df_annual %>%
  filter(
    flag_small_flowering |
      flag_extreme_flowering_valley) %>%
  arrange(
    desc(flag_extreme_flowering_valley),
    dia) %>%
  select(
    site, pop, qu, plant, id, year,
    stage_clean,
    dia, scape,
    invol, nb, nr, herb,
    flag_small_flowering,
    flag_flowering_below1,
    flag_extreme_flowering_valley,
    comment) %>%
  print(n = Inf, width = Inf)