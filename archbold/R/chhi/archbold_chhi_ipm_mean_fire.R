# IPM mean - Archbold -  - Chrysopsis highlandsensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2026.04.20


# Setting the stage ------------------------------------------------------------
set.seed(100)
options(stringsAsFactors = FALSE)


# Packages ---------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(
  MASS,
  tidyverse,
  bbmle,
  patchwork,
  binom,
  skimr,
  lubridate,
  ipmr,
  janitor)


# Specification ----------------------------------------------------------------
v_head <- c('archbold')
v_species <- c('Chrysopsis highlandsensis')
v_custom_delimiter <- c()

v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

v_script_prefix <- str_c(v_head)
v_suffix <- ''

v_ggp_suffix <- paste(
  tools::toTitleCase(v_head), '-', v_species)

v_mod_set_gr <- c()
v_mod_set_su <- c()
v_mod_set_fl <- c()
v_mod_set_fr <- c()


# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_head))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

dir.create(dir_R,      showWarnings = FALSE, recursive = TRUE)
dir.create(dir_data,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_result, showWarnings = FALSE, recursive = TRUE)


# Functions --------------------------------------------------------------------
source('helper_functions/plot_binned_prop.R')
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')


# Data -------------------------------------------------------------------------
df_og <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_data.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    plant_id = identifier,
    year     = year0,
    survival = survival_1) %>%  
  mutate(
    plant_id = as.factor(plant_id)) %>%
  arrange(site, plant_id, year, survival)

df_og_extended <- df_og %>%
  mutate(recruit = ifelse(astg == 1, 1, 0))


# Fire data --------------------------------------------------------------------
df_fire <- read_csv(file.path(dir_data, 'chrysopsis_highlandsensis_fire.csv')) %>% 
  janitor::clean_names() %>% 
  rename(
    year = burn_yr,
    fire = treatment) %>% 
  select(!c(year0, notes))


# Mean data frame --------------------------------------------------------------
df <- df_og_extended %>%
  group_by(site , plant_id, year) %>%
  summarise(
    survives = if (all(is.na(survival))) NA_real_ else min(survival, na.rm = TRUE),
    size_t0  = if (all(is.na(dia)))      NA_real_ else max(dia,      na.rm = TRUE),
    size_t1  = if (all(is.na(dia_1)))    NA_real_ else max(dia_1,    na.rm = TRUE),
    flower   = if (all(is.na(hd)))       NA_real_ else max(hd,       na.rm = TRUE),
    recruit  = if (all(is.na(recruit)))  NA_real_ else min(recruit,  na.rm = TRUE),
    .groups = 'drop') %>% 
  mutate(
    logsize_t0   = log(size_t0),
    logsize_t1   = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3) %>%
  full_join(df_fire, by = c('site', 'year')) %>%
  mutate(
    fire = case_when(
      is.na(fire)    ~ 'No fire',
      fire == 'burn' ~ 'Fire',
      TRUE           ~ NA_character_),
    fire = factor(fire, levels = c('No fire', 'Fire')))


# Prepare data frames ----------------------------------------------------------
df_su <- df %>%
  filter(!is.na(survives), size_t0 != 0, !is.na(fire)) %>%
  select(site, plant_id, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire)


# Survival plotting ------------------------------------------------------------
df_su_binned <- df_su %>%
  group_split(fire) %>%
  purrr::map_df(~ plot_binned_prop(.x, 10, logsize_t0, survives) %>%
                  mutate(fire = unique(.x$fire)))

fig_su_overall <- ggplot(df_su_binned, aes(x = logsize_t0, y = survives, color = fire)) +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives, color = fire),
              position = position_jitter(width = 0.1, height = 0.3), alpha = 0.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c('No fire' = 'black', 'Fire' = 'red')) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = 'top') +
  labs(
    title = 'Survival',
    subtitle = v_ggp_suffix,
    x = expression('log(size)'[t0]),
    y = 'Survival Probability')

fig_su_overall


# Survival models --------------------------------------------------------------
mod_su_0 <- glm(survives ~ fire, data = df_su, family = 'binomial') 
mod_su_1 <- glm(survives ~ logsize_t0 + fire, data = df_su, family = 'binomial') 
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2 + fire, data = df_su, family = 'binomial')  
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire, data = df_su, family = 'binomial')  

mods_su <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- AICtab(mods_su, weights = TRUE, sort = FALSE)$dAIC
mods_su_sorted <- order(mods_su_dAIC)

mod_su_bestfit <- mods_su[[mods_su_sorted[1]]]


# Prediction -------------------------------------------------------------------
df_su_pred <- data.frame(
  logsize_t0   = rep(seq(min(df_su$logsize_t0), max(df_su$logsize_t0), length.out = 100), 2),
  logsize_t0_2 = rep(seq(min(df_su$logsize_t0_2), max(df_su$logsize_t0_2), length.out = 100), 2),
  fire = rep(c('No fire', 'Fire'), each = 100)) %>%
  mutate(
    fire = factor(fire, levels = c('No fire', 'Fire')),
    survives = predict(mod_su_bestfit, newdata = ., type = 'response'))


# Binned data ------------------------------------------------------------------
df_su_binned <- bind_rows(
  plot_binned_prop(filter(df_su, fire == 'No fire'), 10, logsize_t0, survives) %>%
    mutate(fire = 'No fire'),
  plot_binned_prop(filter(df_su, fire == 'Fire'), 10, logsize_t0, survives) %>%
    mutate(fire = 'Fire'))


# Survival plot ----------------------------------------------------------------
# Plot 1
fig_su_line_combined <- ggplot() +
  geom_jitter(data = df_su, aes(x = logsize_t0, y = survives, color = fire),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives, color = fire),
            linewidth = 0.9) +  
  scale_color_manual(values = c('No fire' = 'black', 'Fire' = 'red')) +
  theme_bw() +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Survival to time t1') +
  theme(legend.position = 'none')

# Plot 2 
fig_su_bin_combined <- ggplot() +
  geom_point(data = df_su_binned, aes(x = logsize_t0, y = survives, color = fire)) +
  geom_errorbar(data = df_su_binned, aes(x = logsize_t0, ymin = lwr, ymax = upr, color = fire),
                width = 0.2) +
  geom_line(data = df_su_pred, aes(x = logsize_t0, y = survives, color = fire),
            linewidth = 0.9) +
  scale_color_manual(values = c('No fire' = 'black', 'Fire' = 'red')) +
  theme_bw() +
  ylim(0, 1) +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Survival to time t1') +
  theme(legend.title = element_blank(), legend.position = 'top')

# Combine
fig_su_all <- fig_su_line_combined + fig_su_bin_combined +
  plot_annotation(
    title = 'Survival',
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = 'bold'),
      plot.subtitle = element_text(size = 10, face = 'italic'))  )

fig_su_all


# Growth data ------------------------------------------------------------------
df_gr <- df %>%
  filter(size_t0 != 0,
         !is.na(size_t1),
         size_t1 > 0,
         !is.na(fire)) %>%
  mutate(
    logsize_t0   = log(size_t0),
    logsize_t1   = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3
  ) %>%
  select(plant_id, year, size_t0, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         fire)

fig_gr_overall <- ggplot(df_gr, aes(x = logsize_t0, y = logsize_t1, color = fire)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 8),
    title           = element_text(size = 10),
    plot.subtitle   = element_text(size = 8),
    legend.title    = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title    = "Growth",
    subtitle = v_ggp_suffix,
    x        = expression('log(size)'[t0]),
    y        = expression('log(size)'[t1])
  )

fig_gr_overall



# Growth model -----------------------------------------------------------------
mod_gr_0 <- lm(logsize_t1 ~ fire, data = df_gr)

mod_gr_1 <- lm(logsize_t1 ~ logsize_t0 + fire, data = df_gr)

mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + fire, data = df_gr)

mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire, 
               data = df_gr)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- AICtab(mods_gr, weights = TRUE, sort = FALSE)$dAIC
mods_gr_sorted <- order(mods_gr_dAIC)

if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr + 1
  v_mod_gr_index       <- v_mod_set_gr
}

mod_gr_bestfit <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef   <- coef(mod_gr_bestfit)


# Prediction -------------------------------------------------------------------
df_gr_pred <- data.frame(
  logsize_t0 = rep(
    seq(min(df_gr$logsize_t0, na.rm = TRUE),
        max(df_gr$logsize_t0, na.rm = TRUE), length.out = 100), 2),
  fire = rep(c("No fire", "Fire"), each = 100)
) %>%
  mutate(
    fire = factor(fire, levels = c("No fire", "Fire")),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3   # ✅ ADD THIS
  )

df_gr_pred$logsize_t1 <- predict(mod_gr_bestfit, newdata = df_gr_pred)


# Growth plots ----------------------------------------------------------------- 
# Plot 1: Observed points + prediction lines
fig_gr_line_combined <- ggplot(df_gr, aes(x = logsize_t0, y = logsize_t1, color = fire)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_line(data = df_gr_pred, aes(x = logsize_t0, y = logsize_t1, color = fire),
            linewidth = 1) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  theme_bw() +
  labs(x = expression('log(size)'[t0]),
       y = expression('log(size)'[t1])) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()  )

# Plot 2: Predicted vs observed
fig_gr_pred_combined <- ggplot(df_gr, aes(x = predict(mod_gr_bestfit, newdata = df_gr),
                                          y = logsize_t1, color = fire)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = 'black',
              linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  theme_bw() +
  labs(x = "Predicted", y = "Observed") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()  )

# Combine plots
fig_gr_all_combined <- fig_gr_line_combined + fig_gr_pred_combined +
  plot_annotation(
    title = "Growth Prediction",
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic")))

fig_gr_all_combined



# Growth variance --------------------------------------------------------------
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)  
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2  
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE)) 




# Flowering data ----------------------------------------------------------------
df_fl <- df %>%
  filter(!is.na(flower),
         !is.na(fire),
         !is.na(size_t0),
         !is.na(survives)) %>%
  select(plant_id, year, size_t0, flower, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire) %>%
  mutate(
    flower = if_else(flower > 0, 1, flower)
  )

# Create binned flowering probability by fire group
df_fl_binned <- df_fl %>%
  group_split(fire) %>%
  map_df(~ plot_binned_prop(.x, 10, logsize_t0, flower) %>%
           mutate(fire = unique(.x$fire)))

# Plot overlapped flowering probability
fig_fl_overall <- ggplot(df_fl_binned, aes(x = logsize_t0, y = flower, color = fire)) +
  geom_jitter(
    data = df_fl,
    aes(x = logsize_t0, y = flower, color = fire),
    position = position_jitter(width = 0.1, height = 0.3),
    alpha = 0.1
  ) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9), limits = c(0, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "Flowering probability by fire status",
    subtitle = v_ggp_suffix,
    x = expression('log(size)'[t0]),
    y = "Flowering Probability"
  )

fig_fl_overall


# Flower model -----------------------------------------------------------------
mod_fl_0 <- glm(flower ~ fire, data = df_fl, family = 'binomial') 

mod_fl_1 <- glm(flower ~ logsize_t0 + fire, data = df_fl, family = 'binomial') 

mod_fl_2 <- glm(flower ~ logsize_t0 + logsize_t0_2 + fire,
                data = df_fl, family = 'binomial')  

mod_fl_3 <- glm(flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire,
                data = df_fl, family = 'binomial')  

mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- AICtab(mods_fl, weights = TRUE, sort = FALSE)$dAIC
mods_fl_sorted <- order(mods_fl_dAIC)

if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index       <- mod_fl_index_bestfit - 1 
} else {
  mod_fl_index_bestfit <- v_mod_set_fl + 1
  v_mod_fl_index       <- v_mod_set_fl
}

mod_fl_bestfit <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef   <- coef(mod_fl_bestfit)


# Prediction -------------------------------------------------------------------
df_fl_pred <- data.frame(
  logsize_t0 = rep(seq(min(df_fl$logsize_t0, na.rm = TRUE),
                       max(df_fl$logsize_t0, na.rm = TRUE), length.out = 100), 2),
  fire = rep(c("No fire", "Fire"), each = 100)
) %>%
  mutate(
    fire = factor(fire, levels = c("No fire", "Fire")),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3
  )

df_fl_pred$flower <- predict(mod_fl_bestfit, newdata = df_fl_pred, type = "response")



# Binned observed data for both fire levels
df_fl_binned <- bind_rows(
  plot_binned_prop(filter(df_fl, fire == "No fire"), 10, logsize_t0, flower) %>%
    mutate(fire = "No fire"),
  plot_binned_prop(filter(df_fl, fire == "Fire"), 10, logsize_t0, flower) %>%
    mutate(fire = "Fire")
)



# Plot 1: Raw jitter + prediction lines
fig_fl_line_combined <- ggplot() +
  geom_jitter(data = df_fl, aes(x = logsize_t0, y = flower, color = fire),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower, color = fire),
            linewidth = 0.9) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  theme_bw() +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Flowering Probability') +
  theme(legend.position = "none")



# Plot 2: Binned + prediction
fig_fl_bin_combined <- ggplot() +
  geom_point(data = df_fl_binned, aes(x = logsize_t0, y = flower, color = fire)) +
  geom_errorbar(data = df_fl_binned, aes(x = logsize_t0, ymin = lwr, ymax = upr, color = fire),
                width = 0.2) +
  geom_line(data = df_fl_pred, aes(x = logsize_t0, y = flower, color = fire),
            linewidth = 0.9) +
  scale_color_manual(values = c("No fire" = "black", "Fire" = "red")) +
  theme_bw() +
  ylim(0, 1) +
  labs(title = NULL, x = 'Size at time t0 (log())', y = 'Flowering Probability') +
  theme(legend.title = element_blank(), legend.position = "top")


# Combine ----------------------------------------------------------------------
fig_fl_all <- fig_fl_line_combined + fig_fl_bin_combined +
  plot_annotation(
    title = "Flowering",
    subtitle = v_ggp_suffix,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic")
    )
  )

fig_fl_all



# Number of flowers (conditional on flowering) 
mod_fl_count <- MASS::glm.nb(flower ~ logsize_t0 + logsize_t0_2 + fire, 
                             data = df_fl %>% filter(flower > 0)) 
flc_coefs <- coef(mod_fl_count)


# Recruitment ------------------------------------------------------------------
df_re <- df %>%
  group_by(year, site) %>%
  summarise(
    fire = if_else(any(fire == "Fire", na.rm = TRUE), "Fire", "No fire"),
    tot_p_area = sum(size_t0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  {
    df_quad <- .
    
    df_group <- df_quad %>%
      group_by(year) %>%
      summarise(g_cov = mean(tot_p_area), .groups = "drop")
    
    df_cover <- left_join(df_quad, df_group, by = "year") %>%
      mutate(year = as.integer(year + 1)) %>%
      drop_na()
    
    df_recruit <- df %>%
      group_by(year, site) %>%
      summarise(nr_quad = sum(recruit, na.rm = TRUE), .groups = "drop")
    
    left_join(df_cover, df_recruit, by = c("year", "site"))
  } %>%
  mutate(fire = factor(fire, levels = c("No fire", "Fire")))

# Plot: recruitment vs parent area
ggplot(df_re, aes(x = tot_p_area, y = nr_quad, colour = fire)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1) +  
  theme_bw() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))


# Density dependence -----------------------------------------------------------
df_re_qd <- df %>% 
  group_by(site, year) %>%
  summarise(rec_qd_t1 = sum(recruit, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    df %>% 
      group_by(site, year) %>% 
      summarise(
        nr_ind = sum(!is.na(size_t0)),
        fire = if_else(any(fire == "Fire", na.rm = TRUE), "Fire", "No fire"),
        .groups = "drop") %>% 
      mutate(year = year - 1),
    by = c('site', 'year')) %>%
  mutate(fire = factor(fire, levels = c("No fire", "Fire")))

fig_re_dens <- ggplot(data = df_re_qd) + 
  geom_jitter(aes(y = rec_qd_t1, x = nr_ind, color = fire)) + 
  geom_smooth(aes(y = rec_qd_t1, x = nr_ind, color = fire), method = 'lm') + 
  theme_bw() + 
  labs(title    = 'Recruitment - density dependence: Quad level',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))

fig_re_dens



# Recruitment model ------------------------------------------------------------
df_re_mod <- df_re %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
mod_rec <- MASS::glm.nb(nr_quad ~ fire, data = df_re_mod)

# Generate predictions for recruitment
df_re_mod <- df_re_mod %>% 
  mutate(mod_pred = predict(mod_rec, type = 'response')) 


# Per-capita reproduction ------------------------------------------------------
df_repr_pc <- df %>%
  filter(!is.na(size_t0)) %>% 
  summarize(n_adults = n()) %>%
  bind_cols(
    df_re_mod %>%
      summarize(nr_quad = sum(nr_quad, na.rm = TRUE),
                mod_pred = sum(mod_pred, na.rm = TRUE))) %>%
  mutate(
    repr_pc_mean = mod_pred / n_adults,
    repr_pc_obs = nr_quad / n_adults) %>%
  drop_na()

# overall level
df %>% 
  filter(!is.na(size_t0)) %>% 
  summarize(n_adults = n()) %>% 
  bind_cols(df %>% 
              filter(!is.na(recruit)) %>%
              summarize(n_rec = n())) %>% 
  mutate(rp_pc_m = n_rec / n_adults)


# site level
df %>% 
  filter(!is.na(size_t0)) %>%
  group_by(site) %>% 
  summarize(n_adults = n()) %>% 
  left_join(df %>% 
              filter(!is.na(recruit)) %>%
              group_by(site) %>%
              summarize(n_rec = n()),
            by = 'site') %>% 
  mutate(n_rec = ifelse(is.na(n_rec), 0, n_rec)) %>%
  mutate(repr_pc_mean = n_rec / n_adults) %>% 
  summarise(n_adults = sum(n_adults),
            n_rec    = sum(n_rec),
            rp_pc_m  = mean(repr_pc_mean))


# Extracting parameter estimates -----------------------------------------------
# Survival
coef_su_fe  <- data.frame(coefficient = names(coef(mod_su_bestfit)),
                          value       =       coef(mod_su_bestfit))

coef_su <- Reduce(function(...) rbind(...), list(coef_su_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Growth
coef_gr_fe  <- data.frame(coefficient = names(coef(mod_gr_bestfit)),
                          value       =       coef(mod_gr_bestfit))
coef_gr_var <- data.frame(coefficient = names(coef(mod_gr_var)),
                          value       =       coef(mod_gr_var))

coef_gr <- Reduce(function(...) rbind(...), list(coef_gr_fe, coef_gr_var)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Flower
coef_fl_fe  <- data.frame(coefficient = names(coef(mod_fl_bestfit)),
                          value       =       coef(mod_fl_bestfit))

coef_fl <- Reduce(function(...) rbind(...), list(coef_fl_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# Recruitment 
df_re_size <- df %>% subset(recruit == 1)

# Miscellany
coef_misc   <- data.frame(coefficient = c('rec_siz', 'rec_sd',
                                          'max_siz', 'min_siz'),
                          value       = c(mean(log(df_re_size$size_t0), na.rm = T), 
                                          sd(  log(df_re_size$size_t0), na.rm = T),
                                          df_gr$logsize_t0 %>% max, 
                                          df_gr$logsize_t0 %>% min))

extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = v_script_prefix,
  species = v_species,
  surv_b0 = extr_value(coef_su, 'b0'),
  surv_b1 = extr_value(coef_su, 'logsize_t0'),
  surv_b2 = extr_value(coef_su, 'logsize_t0_2'),
  surv_b3 = extr_value(coef_su, 'logsize_t0_3'),
  surv_bf = extr_value(coef_su, 'fire'),
  grow_b0 = extr_value(coef_gr, 'b0'),
  grow_b1 = extr_value(coef_gr, 'logsize_t0'),
  grow_b2 = extr_value(coef_gr, 'logsize_t0_2'),
  grow_b3 = extr_value(coef_gr, 'logsize_t0_3'),
  grow_bf = extr_value(coef_gr, 'fire'),
  a       = extr_value(coef_gr, 'a'),
  b       = extr_value(coef_gr, 'b'),
  fl_b0   = extr_value(coef_fl, 'b0'),
  fl_b1   = extr_value(coef_fl, 'logsize_t0'),
  fl_b2   = extr_value(coef_fl, 'logsize_t0_2'),
  fl_b3   = extr_value(coef_fl, 'logsize_t0_3'),
  fl_bf   = extr_value(coef_fl, 'fire'),
  recr_sz = extr_value(coef_misc, 'rec_siz'),
  recr_sd = extr_value(coef_misc, 'rec_sd'),
  L       = extr_value(coef_misc, 'min_siz'),
  U       = extr_value(coef_misc, 'max_siz'),
  mat_siz = 200,
  mod_su_index = v_mod_su_index,
  mod_gr_index = v_mod_gr_index,
  mod_gr_index = v_mod_fl_index))


write.csv(pars, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_', v_sp_abb, '_pars.csv'))

# Function describing standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- 0
  for (i in 0:num_pars) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}

# Function describing the invert logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}


# Survival of x-sized individual to time t1
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  for (i in 1:num_pars) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^(i)
    }
  }
  return(inv_logit(survival_value))
}

# Function describing the transition kernel
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Function describing the recruitment 
fy <- function(y, pars, h){
  n_recr  <- pars$fecu_b0
  recr_y  <- dnorm(y, pars$recr_sz, max(h/10, pars$recr_sd)) * h
  recr_y  <- recr_y / sum(recr_y)
  f       <- n_recr * recr_y
  return(f)
}

# Kernel
kernel <- function(pars) {
  
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
  Fmat        <- matrix(0, n, n)
  Fmat[]      <- matrix(fy(y, pars, h), n, n)
  
  # Survival vector
  Smat   <- c()
  Smat   <- sx(y, pars)
  
  # Growth matrix
  Gmat   <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Growth/survival transition matrix
  Tmat   <- matrix(0, n, n)
  
  # Correct for eviction of offspring
  for(i in 1:(n / 2)) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for(i in (n / 2 + 1):n) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrices
  k_yx <- Fmat + Tmat
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y))
}

lambda_ipm <- function(i) {
  return(Re(eigen(kernel(i)$k_yx)$value[1]))
}

# mean population growth rate
lam_mean <- lambda_ipm(pars)
lam_mean

# observed population growth rate
# Population counts at time t0
pop_counts_t0 <- df %>%
  group_by(year, site) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  group_by(year, site) %>%
  summarize(n_t1 = n()) %>% 
  ungroup 

# Calculate observed population growth rates, 
# accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
  # by dropping NAs, we remove gaps in sampling!
  drop_na %>% 
  group_by(year) %>% 
  summarise(n_t0 = sum(n_t0),
            n_t1 = sum(n_t1)) %>% 
  ungroup %>% 
  mutate(obs_pgr = n_t1 / n_t0)

# Geometric mean of yearly population growth rates
lam_mean_count <- exp(mean(log(pop_counts$obs_pgr), na.rm = T))

# Overall (aggregated) population growth rate
lam_mean_overall <- sum(pop_counts$n_t1) / sum(pop_counts$n_t0)


# Building the IPM with ipmr ---------------------------------------------------
proto_ipm_p <- init_ipm(sim_gen   = 'simple',
                        di_dd     = 'di',
                        det_stoch = 'det') %>% 
  define_kernel(
    name      = 'P',
    family    = 'CC',
    formula   = s * g,
    s         = plogis(
      surv_b0 + 
        (if (mod_su_index >= 1) surv_b1 * size_1   else 0) +
        (if (mod_su_index >= 2) surv_b2 * size_1^2 else 0) +
        (if (mod_su_index >= 3) surv_b3 * size_1^3 else 0)),
    
    mu_g      = grow_b0 + 
      (if (mod_gr_index >= 1) grow_b1 * size_1   else 0) +
      (if (mod_gr_index >= 2) grow_b2 * size_1^2 else 0) +
      (if (mod_gr_index >= 3) grow_b3 * size_1^3 else 0),
    
    g         = dnorm(size_2, mu_g, grow_sig),
    grow_sig  = sqrt(a * exp(b * size_1)),
    data_list = pars,
    states    = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions(fun = 'norm', target = 'g')
  ) %>% 
  
  define_kernel(
    name      = 'F',
    family    = 'CC',
    formula   = fecu_b0 * r_d,
    r_d       = dnorm(size_2, recr_sz, recr_sd),
    data_list = pars,
    states    = list(c('size')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'r_d')
  ) %>% 
  
  define_impl(
    make_impl_args_list(
      kernel_names = c(  'P', 'F'),
      int_rule     = rep('midpoint', 2),
      state_start  = rep('size', 2),
      state_end    = rep('size', 2))
  ) %>%
  
  define_domains(
    size = c(pars$L,
             pars$U,
             pars$mat_siz)
  ) %>%
  
  define_pop_state(
    n_size = rep(1 / 200, 200)
  )

ipmr_p <- make_ipm(proto_ipm  = proto_ipm_p, iterations = 200)

lam_mean_ipmr <- lambda(ipmr_p)

lam_out       <- data.frame(coefficient = names(lam_mean_ipmr), 
                            value       = lam_mean_ipmr)
write.csv(lam_out, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_',v_sp_abb, '_lambda_vec.csv'))

lam_out_wide  <- as.list(pivot_wider(lam_out, 
                                     names_from  = 'coefficient', 
                                     values_from = 'value'))

write.csv(lam_out_wide, row.names = F, paste0(
  dir_data, '/', v_script_prefix, '_',v_sp_abb, '_lambda.csv'))
