make_initial_n <- function(year0, pars){
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- (U - L) / n
  b <- L + c(0:n) * h
  
  df0 <- df %>%
    filter(year == year0,
           !is.na(logsize_t0),
           is.finite(logsize_t0))
  
  adult_counts <- hist(
    df0$logsize_t0,
    breaks = b,
    plot = FALSE
  )$counts
  
  adult_density <- adult_counts / h
  
  c(seedbank = pars$seedbank_size, adult_density)
}



# project through observed years------------------------------------------------
project_sequence <- function(years){
  
  n_vec <- make_initial_n(years[1], pars_mean)
  
  out <- data.frame(
    year = years[1],
    projected_adults = sum(n_vec[-1]) * ((pars_mean$U - pars_mean$L) / pars_mean$mat_siz),
    projected_total  = sum(n_vec)
  )
  
  for(i in 1:(length(years) - 1)){
    
    yr <- years[i]
    
    fire_y <- fire_lookup %>%
      filter(year == yr) %>%
      summarise(fire_num = max(fire_num, na.rm = TRUE)) %>%
      pull(fire_num)
    
    pars_y <- pars_yr[[match(yr, years_v)]]
    
    K <- kernel(pars_y, fire = fire_y)$k_yx
    
    n_vec <- K %*% n_vec
    
    out <- rbind(
      out,
      data.frame(
        year = years[i + 1],
        projected_adults = sum(n_vec[-1]) * ((pars_mean$U - pars_mean$L) / pars_mean$mat_siz),
        projected_total  = sum(n_vec)
      )
    )
  }
  
  out
}

proj <- project_sequence(1999:2024)

proj


# compare projected adults to observed adults ----------------------------------
obs_counts <- df %>%
  filter(!is.na(logsize_t0),
         is.finite(logsize_t0)) %>%
  group_by(year) %>%
  summarise(
    observed_adults = n(),
    .groups = "drop"
  )

proj_compare <- proj %>%
  left_join(obs_counts, by = "year") %>%
  mutate(
    projected_lambda = projected_adults / lag(projected_adults),
    observed_lambda  = observed_adults / lag(observed_adults)
  )

proj_compare


proj_compare %>%
  select(year, observed_adults, projected_adults, observed_lambda, projected_lambda)



# inspect annual projected recruits --------------------------------------------
project_sequence <- function(years){
  
  n_vec <- make_initial_n(years[1], pars_mean)
  h <- (pars_mean$U - pars_mean$L) / pars_mean$mat_siz
  
  out <- data.frame(
    year = years[1],
    projected_adults = sum(n_vec[-1]) * h,
    projected_total  = sum(n_vec),
    projected_recruits = NA_real_
  )
  
  for(i in 1:(length(years) - 1)){
    
    yr <- years[i]
    
    fire_y <- fire_lookup %>%
      filter(year == yr) %>%
      summarise(fire_num = max(fire_num, na.rm = TRUE)) %>%
      pull(fire_num)
    
    pars_y <- pars_yr[[match(yr, years_v)]]
    K <- kernel(pars_y, fire = fire_y)$k_yx
    
    projected_recruits <- sum(K[-1, 1] * n_vec[1]) * h
    
    n_vec <- K %*% n_vec
    
    out <- rbind(
      out,
      data.frame(
        year = years[i + 1],
        projected_adults = sum(n_vec[-1]) * h,
        projected_total  = sum(n_vec),
        projected_recruits = projected_recruits
      )
    )
  }
  
  out
}

proj <- project_sequence(1999:2024)

proj %>%
  select(year, projected_adults, projected_recruits)



# compare projection to observed counts ----------------------------------------
obs_counts <- df %>%
  filter(!is.na(logsize_t0),
         is.finite(logsize_t0)) %>%
  group_by(year) %>%
  summarise(
    observed_adults = n(),
    .groups = "drop"
  )

proj_compare <- proj %>%
  left_join(obs_counts, by = "year") %>%
  mutate(
    projected_lambda = projected_adults / lag(projected_adults),
    observed_lambda  = observed_adults / lag(observed_adults)
  )

proj_compare %>%
  select(year, observed_adults, projected_adults,
         observed_lambda, projected_lambda) %>% 
  ggplot() +
  geom_point(aes(x = projected_adults, y = observed_adults)) +
  theme_classic() +
  geom_abline(aes(intercept = 0, slope = 1))



# quantify whether the year effects --------------------------------------------
cor(
  proj_compare$observed_lambda,
  proj_compare$projected_lambda,
  use = "complete.obs")
lm(
  observed_lambda ~ projected_lambda,
  data = proj_compare) %>%
  summary()


# quantify whether recruitment is driving the observed dynamics ----------------
df_rec <- df %>%
  group_by(year) %>%
  summarise(
    recruits = sum(recruit, na.rm = TRUE),
    .groups = "drop")

proj_compare2 <- proj_compare %>%
  left_join(df_rec, by = "year")

cor(
  proj_compare2$observed_lambda,
  proj_compare2$recruits,
  use = "complete.obs")


lm(
  observed_lambda ~ recruits,
  data = proj_compare2) %>%
  summary()


# quantify how much recruitment varies among years -----------------------------
df_rec %>%
  summarise(
    mean_rec = mean(recruits),
    sd_rec   = sd(recruits),
    cv_rec   = sd_rec / mean_rec)
summary(df_rec$recruits)


# Does recruitment explain observed population growth after accounting for survival/growth quality? ----
df_lambda_T <- data.frame(
  year = years_v,
  lambda_T = sapply(
    seq_along(pars_yr),
    function(i) {
      
      yr <- years_v[i]
      
      fire_y <- fire_lookup %>%
        filter(year == yr) %>%
        summarise(fire_num = max(fire_num, na.rm = TRUE)) %>%
        pull(fire_num)
      
      Re(eigen(kernel(pars_yr[[i]], fire = fire_y)$Tmat)$values[1])
    }))

analysis_df <- proj_compare2 %>%
  left_join(df_lambda_T, by = "year") %>%
  select(year, observed_lambda, projected_lambda, lambda_T, recruits) %>%
  filter(!is.na(observed_lambda))

analysis_df

m1 <- lm(observed_lambda ~ lambda_T, data = analysis_df)
m2 <- lm(observed_lambda ~ recruits, data = analysis_df)
m3 <- lm(observed_lambda ~ lambda_T + recruits, data = analysis_df)

summary(m1)
summary(m2)
summary(m3)
summary(m1)$r.squared
summary(m2)$r.squared
summary(m3)$r.squared
summary(m3)$coefficients


plot(analysis_df$recruits, analysis_df$observed_lambda)

abline(m2)

plot(analysis_df$lambda_T, analysis_df$observed_lambda)

abline(m1)

pairs(
  analysis_df[, c("observed_lambda", "lambda_T", "recruits")])

