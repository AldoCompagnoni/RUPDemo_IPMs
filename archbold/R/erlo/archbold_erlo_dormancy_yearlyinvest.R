warnings()

cat("\nSURVIVAL\n")
formula(mod_su_bestfit)
mods_su_dAICc

cat("\nGROWTH\n")
formula(mod_gr_bestfit)
mods_gr_dAICc

cat("\nDORMANCY\n")
formula(mod_do_bestfit)
mods_do_dAICc

cat("\nREACTIVATION\n")
formula(mod_ra_bestfit)
mods_ra_dAICc

cat("\nFLOWERING\n")
formula(mod_fl_bestfit)
mods_fl_dAICc

cat("\nSCAPES\n")
formula(mod_fl_n_bestfit)
mods_fl_n_dAICc

cat("\nRECRUITMENT\n")
formula(mod_re_bestfit)
mods_re_dAICc
re_year_diagnostic

c(
  recr_per_scape_0 = recr_per_scape_0,
  recr_per_scape_1 = recr_per_scape_1)



df_shift_check <- map_dfr(-2:2, function(shift_i) {
  
  df_i <- df_compare %>%
    select(model_year = year, asym_lambda, proj_lambda) %>%
    mutate(obs_year = model_year + shift_i) %>%
    left_join(
      df_compare %>% select(obs_year = year, obs_pgr),
      by = "obs_year") %>%
    drop_na()
  
  tibble(
    shift = shift_i,
    n = nrow(df_i),
    cor_asym = cor(df_i$asym_lambda, df_i$obs_pgr),
    cor_projected = cor(df_i$proj_lambda, df_i$obs_pgr),
    rmse_asym = sqrt(mean(
      (df_i$asym_lambda - df_i$obs_pgr)^2)),
    rmse_projected = sqrt(mean(
      (df_i$proj_lambda - df_i$obs_pgr)^2)))
})

df_shift_check


get_n_model_state <- function(year_i, site_i, pop_i, qu_i) {
  
  pars_i <- make_ipm_pars(year_i)
  
  sum(make_initial_n_quad(
    year0 = year_i,
    site_i = site_i,
    pop_i = pop_i,
    qu_i = qu_i,
    pars = pars_i))
}

df_compare_same_state <- df_compare_quad %>%
  mutate(
    n_t1_model_state = pmap_dbl(
      list(year + 1, site, pop, qu),
      get_n_model_state)) %>%
  group_by(year) %>%
  summarise(
    n_t0_model_state = sum(n_obs_model),
    n_t1_model_state = sum(n_t1_model_state),
    n_projected = sum(n_proj_model),
    obs_lambda_same_state =
      n_t1_model_state / n_t0_model_state,
    proj_lambda = n_projected / n_t0_model_state,
    .groups = "drop") %>%
  mutate(error = proj_lambda - obs_lambda_same_state)

df_compare_same_state %>%
  print(n = Inf)

df_compare_same_state %>%
  summarise(
    correlation = cor(obs_lambda_same_state, proj_lambda),
    mean_error = mean(error),
    rmse = sqrt(mean(error^2)))



#-------------------------------------------------------------------------------
df_compare_same_state %>%
  summarise(
    correlation = cor(
      obs_lambda_same_state, proj_lambda,
      use = "complete.obs"),
    mean_error = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)))

df_compare_quad %>%
  filter(!is.finite(n_proj_model)) %>%
  select(
    year, site, pop, qu, n_t0, n_t1, n_obs_model,
    disturbance_num, disturbance_prev_num) %>%
  print(n = Inf)

project_components <- function(
    yr, site_i, pop_i, qu_i, disturbance_y, disturbance_prev_y) {
  
  pars_y <- make_ipm_pars(yr)
  
  n0 <- make_initial_n_quad(
    yr, site_i, pop_i, qu_i, pars_y)
  
  if (sum(n0) <= 0) {
    return(tibble(
      year = yr, site = site_i, pop = pop_i, qu = qu_i,
      n_initial = sum(n0), proj_survival = NA_real_,
      proj_recruits = NA_real_))
  }
  
  k <- kernel(
    pars_y, disturbance = disturbance_y,
    disturbance_prev = disturbance_prev_y)
  
  n <- pars_y$mat_siz
  active0 <- n0[seq_len(n)]
  dormant0 <- n0[n + 1]
  
  proj_survival <-
    sum(k$Tmat %*% active0) +
    sum(k$Dorm_row %*% active0) +
    dormant0
  
  proj_recruits <- sum(k$Fmat %*% active0)
  
  tibble(
    year = yr, site = site_i, pop = pop_i, qu = qu_i,
    n_initial = sum(n0),
    proj_survival = proj_survival,
    proj_recruits = proj_recruits)
}

df_components <- pmap_dfr(
  df_compare_quad %>%
    transmute(
      yr = year,
      site_i = site,
      pop_i = pop,
      qu_i = qu,
      disturbance_y = disturbance_num,
      disturbance_prev_y = disturbance_prev_num),
  project_components)

df_components_year <- df_components %>%
  group_by(year) %>%
  summarise(
    n_initial = sum(n_initial, na.rm = TRUE),
    proj_survival = sum(proj_survival, na.rm = TRUE),
    proj_recruits = sum(proj_recruits, na.rm = TRUE),
    proj_total = proj_survival + proj_recruits,
    lambda_survival = proj_survival / n_initial,
    lambda_recruitment = proj_recruits / n_initial,
    lambda_total = proj_total / n_initial,
    .groups = "drop") %>%
  left_join(
    df_compare_same_state %>%
      select(year, n_t1_model_state, obs_lambda_same_state),
    by = "year")

df_components_year %>%
  filter(year %in% c(1991, 1992, 1995)) %>%
  print(n = Inf, width = Inf)