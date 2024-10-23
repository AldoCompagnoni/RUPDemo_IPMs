# IPM for alder 2007 & chu 2016; kansas; Bouteloua curtipendula

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com
# Date: 2024.10.23

# Process plant growth data to explore vital rates for 
# the species Sporobolus asper and 
# generates visualizations and model fits for IPM analysis.

# Setting the stage ------------------------------------------------------------
# Clear the workspace by removing all objects in the global environment
rm(list = ls()) 
# Set a random seed for reproducibility of results
set.seed(100)
# Ensure strings are treated as characters rather than factors
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork, skimr, ipmr)


# Data -------------------------------------------------------------------------
# Specify the species
species <- 'Bouteloua curtipendula'

# Create a unique species abbreviation for file naming
sp_abb  <- tolower(
  gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                      collapse = "")))

# Load the plant tracking data (this line is commented out for now)
# source(paste0('adler_2007_ks/R/', sp_abb, '/kansas_tracker_', sp_abb, '.R'))

## Read and clean the species data
df <- read.csv(paste0('adler_2007_ks/data/', 
                      sp_abb, '/ks_', sp_abb, '.csv')) %>% 
  filter(Species == species) %>%
  select(-c(Suspect, nearEdge, Site)) %>%
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
# Provide a summary of the cleaned data
skim(df)

# Data exploration -------------------------------------------------------------
# Analyze the quadrat inventory by year
# Group data by quadrant and year, counting the number of individuals
inv_plot_per_year <- df %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # Create a scatter plot of quadrat counts over the years
  ggplot() + 
  geom_point(aes(x = year, y = quad)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 5))
inv_plot_per_year

# Create results directory if it doesn't exist
if (!dir.exists(paste0("adler_2007_ks/results/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/results/", sp_abb))
}
# Save the inventory plot to a file
ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, "/overall_inventory_quadrat_per_year.png"),  
       plot = inv_plot_per_year,
       width = 6, height = 4, dpi = 150)

# Prepare data frames for analysis ---------------------------------------------
# Survival data frame
surv_df <- 
  subset(df, !is.na(survives)) %>%
  subset(size_t0 != 0) %>%
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Growth data frame
grow_df <- 
  df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Total area data frame
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

# Recruitment data frame
recr_df <- 
  df %>%   
  group_by (year, quad) %>% 
  summarise(nr_quad = sum(recruit, na.rm = T)) %>% 
  ungroup

recr_df <- left_join(cover_df, recr_df)

# Create the data folder for the species if it does not exist already
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))}

write.csv(df,      paste0("adler_2007_ks/data/", sp_abb, "/data_df.csv"))
write.csv(surv_df, paste0("adler_2007_ks/data/", sp_abb, "/survival_df.csv"))
write.csv(grow_df, paste0("adler_2007_ks/data/", sp_abb, "/growth_df.csv"))
write.csv(recr_df, paste0("adler_2007_ks/data/", sp_abb, "/recruitment_df.csv"))

# Plotting the data --------------------------------------------------------
# Create histograms for log-transformed sizes at t0 and t1
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

# Survival analysis
# Load custom function for plotting binned proportions
source('helper_functions/plot_binned_prop.R')

# Generate a plot for overall survival based on size at t0
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

# Save the survival plot to a file
ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_surv.png"), 
       plot = surv_overall, 
       width = 4, height = 3, units = "in", dpi = 150)

# Growth analysis
# Create a scatter plot of size at t0 versus size at t1
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

# Recruitment analysis
# Create a scatter plot showing the relationship between 
# total parent plant area and number of recruits
rec_overall <- 
  ggplot(recr_df, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(x = expression('Total parent plant area '[t0]),   
       y = expression('Number of recruits '     [t1]))  

# Save the recruitment plot as a PNG file
ggsave(paste0("adler_2007_ks/results/", sp_abb, "/overall_rec.png"), 
       plot = rec_overall, 
       width = 4, height = 3, units = "in", dpi = 150)


# Fit vital rate models for the mean IPM -----------------------------------
# Fit growth models to predict size at time t1 based on size at time t0
# Linear model
gr_mod_mean   <- 
  lm(logsize_t1 ~ logsize_t0, data = grow_df)
# Quadratic model
gr_mod_mean_2 <- 
  lm(logsize_t1 ~ logsize_t0 + logsize_t0_2, data = grow_df)  
# Cubic model
gr_mod_mean_3 <- 
  lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, data = grow_df)  
# Compare models using AIC
gr_mods       <- list(gr_mod_mean, gr_mod_mean_2, gr_mod_mean_3)
gr_aic        <- sapply(gr_mods, AIC)
# Assign the best model
gr_mod_bestfit_index <- which.min(gr_aic)
gr_mod_bestfit       <- gr_mods[[gr_mod_bestfit_index]]
gr_ranef             <- coef(gr_mod_bestfit)

# Predict size at time t1 using the mean growth model
grow_df$pred <- predict(gr_mod_bestfit, type = "response")

# Plot observed size at time t1 against size at time t0 with the fitted line
source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

grow_line <- 
  ggplot(grow_df, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, gr_ranef), 
                color = line_color_pred_fun(gr_ranef), 
                lwd = 2) +
  theme_bw()

# Plot predicted versus observed size at time t1
grow_pred <- 
  ggplot(grow_df, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = "red", lwd = 2) + 
  theme_bw()

# Combine growth line and prediction plots
grow_overall_pred <- grow_line + grow_pred + plot_layout() 

# Save the growth prediction plot
ggsave(paste0("adler_2007_ks/results/", sp_abb, 
              "/overall_grow_pred_logs", gr_mod_bestfit_index, ".png"), 
       plot = grow_overall_pred, 
       width = 8, height = 4, units = "in", dpi = 150)

# Fit a model to assess variance in growth
# Fitted values from growth model
x         <- fitted(gr_mod_mean)  
# Squared residuals
y         <- resid(gr_mod_mean)^2  
# Non-linear model for variance
gr_var_m  <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0))  

# Survival
# Fit models to predict survival based on size at time t0
# Logistic regression
su_mod_mean   <- glm(survives ~ logsize_t0, 
                     data = surv_df, family = "binomial") 
# Quadratic logistic model
su_mod_mean_2 <- glm(survives ~ logsize_t0 + logsize_t0_2, 
                     data = surv_df, family = "binomial")  
# Cubic logistic model
su_mod_mean_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
                     data = surv_df, family = "binomial")  
# Compare models using AIC
su_mods       <- list(su_mod_mean, su_mod_mean_2, su_mod_mean_3)
su_aic        <- sapply(su_mods, AIC)
# Assign the best model
su_mod_bestfit_index <- which.min(su_aic)
su_mod_bestfit       <- su_mods[[su_mod_bestfit_index]]
su_ranef             <- coef(su_mod_bestfit)

# Generate predictions for survival across a range of sizes
surv_x <- seq(min(surv_df$logsize_t0, na.rm = T), 
              max(surv_df$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
surv_pred_df <- predictor_fun(surv_x, su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = surv_x, survives = .)

# Plot observed survival with fitted line
surv_line <- 
  ggplot() +
  geom_jitter(data = surv_df, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0, height = 0.25) + 
  geom_line(data = surv_pred_df, aes(x = logsize_t0, 
                                     y = survives),
            color = line_color_pred_fun(su_ranef), 
            lwd   = 2) +  
  theme_bw()

# Plot binned survival proportions with error bars
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

# Combine survival plots
surv_overall_pred <- surv_line + surv_bin + plot_layout()

# Save the survival prediction plot
ggsave(paste0("adler_2007_ks/results/", sp_abb, 
              "/overall_surv_pred_logs", su_mod_bestfit_index, ".png"), 
       plot = surv_overall_pred, 
       width = 8, height = 3, units = "in", dpi = 150) 

## Recruitment
# Filter recruitment data to exclude NAs
recr_nona_nr_quad <- recr_df %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
rec_mod_mean <- MASS::glm.nb(nr_quad ~ 1, data = recr_nona_nr_quad)

# Generate predictions for recruitment
recr_nona_nr_quad <- 
  recr_nona_nr_quad %>% 
  mutate(pred_mod_mean = predict(rec_mod_mean, type = "response")) 

# Summarize total number of recruits and predictions
rec_sums_df_m <- 
  recr_nona_nr_quad %>%
  summarize(nr_quad = sum(nr_quad),
            pred_mod_mean = sum(pred_mod_mean))

# Count number of adult individuals
indiv_m <- surv_df %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
repr_pc_m <- indiv_m %>%
  bind_cols(rec_sums_df_m) %>%
  mutate(repr_pc_mean = pred_mod_mean / n_adults) %>%
  mutate(repr_pc_obs = nr_quad / n_adults) %>%
  drop_na 

# Output the reproduction per capita summary
repr_pc_m