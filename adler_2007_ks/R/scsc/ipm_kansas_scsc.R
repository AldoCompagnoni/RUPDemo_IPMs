# IPM for alder 2007 & chu 2016; kansas; Schizachyrium scoparium

# Niklas Neisse
# 2024.10.07

# reading in, and cleaning the data
# exploring the overall-years rates
# setting up the vital rate data-frames for the year specific 

# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source( 'helper_functions/load_packages.R' )
load_packages( tidyverse, patchwork, skimr )


# Data -------------------------------------------------------------------------
# Define the species variable
species <- "Schizachyrium scoparium"
sp_abb  <- tolower(gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                                       collapse = "")))
## Species data frame
df <- read.csv('https://www.dropbox.com/scl/fi/dkrh3vd2d2iku6837a7kh/KS_SCSC_ANGE.csv?rlkey=52e6iouahfhzlpx4mhkm7lpjf&dl=1') %>% 
  filter(Species == species) %>%
  select(-c(Suspect, nearEdge, Site)) %>% # geometry, 
  mutate(across(c(Quad), as.factor)) %>% 
  rename(species  = Species, 
         size_t0  = basalArea_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1,
         year     = Year,
         quad     = Quad,
         track_id = trackID) %>% 
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),
         logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3)
skim(df)


# Data exploration -------------------------------------------------------------

## Quadrat inventory
# grouping by quard, year and counting the individuals
inv_plot_per_year <- df %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>% 
  # ploting year against quadrat
  ggplot() + 
  geom_point(aes(x = year, y = quad)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 5))
inv_plot_per_year

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_inventory_quadrat_per_year.png"),  
       plot = inv_plot_per_year,
       width = 6, height = 4, dpi = 150)


# Getting the dfs ready --------------------------------------------------------

## Survival
surv_df <- 
  subset(df, !is.na(survives)) %>%
  subset(size_t0 != 0) %>%
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

## Growth
grow_df <- 
  df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

## Total area 
quad_df <- 
  df %>%  
  group_by (quad, year) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

group_df <- 
  quad_df %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

cover_df <- 
  left_join(quad_df, group_df) %>% 
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

## Number of recruits
recr_df <- 
  df %>%   
  group_by (year, quad) %>% 
  summarise(nr_quad = sum(recruit, na.rm = T)) %>% 
  ungroup

recr_df <- left_join(cover_df, recr_df)

# Create the data folder for the species if it does not exits already
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))}
# Save all data
write.csv(df,      paste0("adler_2007_ks/data/", sp_abb, "/data_df.csv"))
write.csv(surv_df, paste0("adler_2007_ks/data/", sp_abb, "/survival_df.csv"))
write.csv(grow_df, paste0("adler_2007_ks/data/", sp_abb, "/growth_df.csv"))
write.csv(recr_df, paste0("adler_2007_ks/data/", sp_abb, "/recruitment_df.csv"))

# Plotting the data --------------------------------------------------------

## sizes at log scale 
hist_t0 <- 
  ggplot(df, aes(x = logsize_t0)) +
  geom_histogram(binwidth = 0.2, fill = "grey", color = "black") +
  labs(title = "Histogram of size at time t0", x = "Size at time t0") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

hist_t1 <- 
  ggplot(df, aes(x = logsize_t1)) +
  geom_histogram(binwidth = 0.2, fill = "white", color = "black") +
  labs(title = "Histogram of size at time t1", x = "Size at time t1") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

hist_sizes_log <- hist_t0 + hist_t1

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_hist_sizes_log.png"), 
       plot = hist_sizes_log, 
       width = 8, height = 3, units = "in", dpi = 150)

## Survival
# function removes NAs, and provides standard errors, 
source('helper_functions/plot_binned_prop.R')

surv_overall <- 
  ggplot(data = plot_binned_prop(df, 10, logsize_t0, survives)) +
  geom_point(aes(x = logsize_t0, 
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(x = expression('log(size)'[t0]),
       y = expression('Survival to time t1'))


ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_surv.png"), 
       plot = surv_overall, 
       width = 4, height = 3, units = "in", dpi = 150)

## Growth
gr_overall <-
  ggplot(data  = grow_df, aes( x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(x = expression('log(size) ' [t0]),
       y = expression('log(size)  '[t1]))

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_gr.png"), 
       plot = gr_overall, 
       width = 4, height = 3, units = "in", dpi = 150)

## Recruitment
rec_overall <- 
  ggplot(recr_df, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +
  theme_bw() +
  labs(x = expression('Total parent plant area '[t0]),
       y = expression('Number of recruits '     [t1]))

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_rec.png"), 
       plot = rec_overall, 
       width = 4, height = 3, units = "in", dpi = 150)



# Fitting vital rate models for the mean IPM -----------------------------------

# Growth
gr_mod_mean  <- lm(logsize_t1 ~ logsize_t0, data = grow_df)
grow_df$pred <- predict( gr_mod_mean, type = "response" )

grow_line <- 
  ggplot(grow_df, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point() +
  geom_abline(aes(intercept = coef(gr_mod_mean)[1],
                  slope     = coef(gr_mod_mean)[2]),
              color= 'red', lwd = 2) +
  theme_bw()

grow_pred <- 
  ggplot(grow_df, aes(x = pred, y = logsize_t1)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), 
              color = "red", lwd = 2) + 
  theme_bw()

grow_overall_pred <- grow_line + grow_pred + plot_layout() 

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_grow_pred.png"), 
       plot = grow_overall_pred, 
       width = 8, height = 4, units = "in", dpi = 150)

x         <- fitted(gr_mod_mean)
y         <- resid(gr_mod_mean)^2
gr_var_m  <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0))


# Survival
su_mod_mean   <- glm(survives ~ logsize_t0, 
                     data = surv_df, family = "binomial" )
su_mod_mean_2 <- glm(survives ~ logsize_t0 + logsize_t0_2, 
                     data = surv_df, family = "binomial" )
su_mod_mean_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
                     data = surv_df, family = "binomial" )

surv_x <- seq(min(surv_df$logsize_t0, na.rm = T), 
              max(surv_df$logsize_t0, na.rm = T), length.out = 100)
surv_pred <- boot::inv.logit(coef(su_mod_mean)[1] + coef(su_mod_mean)[2] * surv_x)

surv_pred_df <- data.frame(logsize_t0 = surv_x, survives = surv_pred)

surv_line <- 
  ggplot() +
  geom_jitter(data = surv_df, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0, height = 0.25) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2 ) +
  theme_bw()

surv_bin <- 
  ggplot() +
  geom_point(data =  plot_binned_prop(df, 10, 
                                      logsize_t0, survives), 
             aes(x = logsize_t0, 
                 y = survives) ) +
  geom_errorbar(data =  plot_binned_prop(df, 10, 
                                         logsize_t0, survives), 
                aes(x = logsize_t0, 
                    ymin = lwr,
                    ymax = upr) ) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2) +
  theme_bw()

surv_overall_pred <- surv_line + surv_bin + plot_layout()

ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_surv_pred.png"), 
       plot = surv_overall_pred, 
       width = 8, height = 3, units = "in", dpi = 150) 

## Recruitment
recr_nona_nr_quad <- recr_df %>% filter(!is.na(nr_quad))
rec_mod_mean <- MASS::glm.nb(nr_quad ~ 1, data = recr_nona_nr_quad)

recr_nona_nr_quad <- 
  recr_nona_nr_quad %>% 
  mutate(pred_mod_mean = predict(rec_mod_mean, type = "response")) 

rec_sums_df_m <- 
  recr_nona_nr_quad %>%
  summarize(nr_quad = sum(nr_quad),
            pred_mod_mean = sum(pred_mod_mean))

indiv_m <- surv_df %>%
  summarize(n_adults = n())

repr_pc_m <- indiv_m %>%
  bind_cols(rec_sums_df_m) %>%
  mutate(repr_pc_mean = pred_mod_mean / n_adults) %>%
  mutate(repr_pc_obs = nr_quad / n_adults) %>%
  drop_na

repr_pc_m
