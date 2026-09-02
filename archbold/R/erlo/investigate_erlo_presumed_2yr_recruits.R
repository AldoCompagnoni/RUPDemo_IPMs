# Investigate presumed 2-year recruits ------------------------------------------
# Eriogonum longifolium
#
# Compares:
# 1) presumed 2-year recruits = newly observed adults with stage_entry == 2
# 2) known recruits          = observed seedlings
# 3) adults                  = established active individuals that are not
#                              entering as recruits in that census
#
# The original workdata script creates df_recruit_entry and df_transition.

library(tidyverse)
library(patchwork)

# Run the workdata preparation first -------------------------------------------
source("C:/code/RUPDemo_IPMs/archbold/R/erlo/archbold_erlo_workdata_260819.R")


# Recruit identities -----------------------------------------------------------
df_recruit_groups <- df_recruit_entry %>%
  filter(recruit_class %in% c("observed_seedling", "presumed_2yr")) %>%
  transmute(
    id,
    year,
    group = recode(
      recruit_class,
      observed_seedling = "Known recruit",
      presumed_2yr = "Presumed 2yr new adult"))


# Established adults -----------------------------------------------------------
# Active plants that are not entering as either recruit group in that census.
df_adults <- df_transition %>%
  filter(
    state_clean == "active",
    size_t0 > 0,
    is.finite(logsize_t0)) %>%
  anti_join(
    df_recruit_groups %>%
      select(id, year),
    by = c("id", "year")) %>%
  transmute(
    id, year,
    group = "Adult")


# Combine groups with demographic transitions ---------------------------------
df_investigate <- bind_rows(
  df_recruit_groups,
  df_adults) %>%
  left_join(
    df_transition %>%
      select(
        id, year, site, pop, qu, plant,
        state_clean, state_t1,
        stage_clean, stage_t1,
        size_t0, size_t1,
        logsize_t0, logsize_t1,
        annual_transition),
    by = c("id", "year")) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Presumed 2yr new adult",
        "Known recruit",
        "Adult")))


# Basic counts -----------------------------------------------------------------
df_investigate %>%
  count(group)

df_investigate %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_size_t0 = sum(!is.na(size_t0)),
    n_growth = sum(!is.na(size_t0) & !is.na(size_t1)),
    mean_size_t0 = mean(size_t0, na.rm = TRUE),
    median_size_t0 = median(size_t0, na.rm = TRUE),
    min_size_t0 = min(size_t0, na.rm = TRUE),
    max_size_t0 = max(size_t0, na.rm = TRUE),
    .groups = "drop")


# Size distribution ------------------------------------------------------------
df_size <- df_investigate %>%
  filter(
    size_t0 > 0,
    is.finite(logsize_t0))

fig_size_density <- ggplot(
  df_size,
  aes(x = logsize_t0, color = group, fill = group)) +
  geom_density(alpha = 0.18, linewidth = 1) +
  labs(
    title = "Entry-size distribution",
    subtitle = "Presumed 2yr recruits vs known recruits vs established adults",
    x = expression("log(diameter)"[t0]),
    y = "Density",
    color = NULL,
    fill = NULL) +
  theme_bw()

fig_size_density


fig_size_raw <- ggplot(
  df_size,
  aes(x = group, y = size_t0, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.4) +
  scale_y_log10() +
  labs(
    title = "Observed size at entry",
    x = NULL,
    y = expression("Diameter"[t0]~"(log scale)"),
    color = NULL) +
  theme_bw() +
  theme(legend.position = "none")

fig_size_raw


# Growth to next year -----------------------------------------------------------
df_growth <- df_investigate %>%
  filter(
    annual_transition,
    state_t1 == "active",
    size_t0 > 0,
    size_t1 > 0,
    is.finite(logsize_t0),
    is.finite(logsize_t1))

df_growth %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_log_growth = mean(logsize_t1 - logsize_t0),
    median_log_growth = median(logsize_t1 - logsize_t0),
    mean_ratio = mean(size_t1 / size_t0),
    median_ratio = median(size_t1 / size_t0),
    .groups = "drop")


# Same basic form as the IPM growth plot, now colored by demographic origin.
fig_growth <- ggplot(
  df_growth,
  aes(
    x = logsize_t0,
    y = logsize_t1,
    color = group)) +
  geom_point(alpha = 0.45) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2) +
  labs(
    title = "Growth to the next annual census",
    subtitle = "Separate observed growth relationships for the three groups",
    x = expression("log(diameter)"[t0]),
    y = expression("log(diameter)"[t1]),
    color = NULL) +
  theme_bw()

fig_growth


# Direct growth change ----------------------------------------------------------
df_growth <- df_growth %>%
  mutate(
    log_growth = logsize_t1 - logsize_t0,
    size_ratio = size_t1 / size_t0)

fig_growth_change <- ggplot(
  df_growth,
  aes(x = group, y = log_growth, color = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.4) +
  labs(
    title = "Annual change in log size",
    x = NULL,
    y = expression("log(diameter)"[t1] - "log(diameter)"[t0]),
    color = NULL) +
  theme_bw() +
  theme(legend.position = "none")

fig_growth_change


# Combined diagnostic -----------------------------------------------------------
fig_investigate <- (
  fig_size_density + fig_size_raw
) / (
  fig_growth + fig_growth_change
) +
  plot_annotation(
    title = "Investigation of presumed 2-year recruits")

fig_investigate


# Inspect the presumed 2-year individuals directly -----------------------------
df_presumed_2yr <- df_investigate %>%
  filter(group == "Presumed 2yr new adult") %>%
  arrange(year, site, pop, qu, plant)

df_presumed_2yr %>%
  select(
    site, pop, qu, plant, id, year,
    stage_clean, stage_t1,
    size_t0, size_t1,
    state_t1)

# Individuals with usable next-year growth
df_presumed_2yr %>%
  filter(!is.na(size_t1)) %>%
  select(
    site, pop, qu, plant, id, year,
    size_t0, size_t1,
    stage_clean, stage_t1)
