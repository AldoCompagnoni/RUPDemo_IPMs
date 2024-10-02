# IPM for alder 2007 & chu 2016; kansas; bocu

# Niklas Neisse
# 2024.09.30

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

# Define CRAN packages
.cran_packages <- c('tidyverse','patchwork','skimr') 
# Check if CRAN packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}
# Load required packages
sapply(.cran_packages, require, character.only = TRUE) 


# Data -------------------------------------------------------------------------
## Bouteloua  curtipendula (bocu) data frame
df <- read.csv("adler_2007_ks/data/ks_grasses.csv") %>%  
  filter(Species == "Bouteloua curtipendula") %>%
  select(-c(geometry, Suspect, nearEdge, Species, Site)) %>% 
  mutate(across(c(Quad), as.factor)) %>% 
  rename(size_t0 = basalArea_genet,
         size_t1 = size_tplus1,
         survives = survives_tplus1,
         year = Year,
         quad = Quad,
         track_id = trackID) %>% 
  mutate(logsize_t0 = log(size_t0),
         logsize_t1 = log(size_t1),
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
  theme_bw()
inv_plot_per_year

ggsave("adler_2007_ks/results/bocu/overall_inventory_quadrat_per_year.png",  
       plot = inv_plot_per_year,
       width = 6, height = 4, dpi = 150)


# Getting the dfs ready --------------------------------------------------------

## Survival
surv_df <- 
  subset(df, !is.na(survives)) %>%
  subset(size_t0 != 0) %>%
  select(quad, track_id, year, 
         size_t0, survives, size_t1, 
         logsize_t0, logsize_t1,
         logsize_t0_2, logsize_t0_3)

## Growth
grow_df <- 
  df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(quad, track_id, year, 
         size_t0, survives, size_t1, 
         logsize_t0, logsize_t1,
         logsize_t0_2, logsize_t0_3)

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

recr_df <- 
  left_join(cover_df, recr_df)

write.csv(df,      "adler_2007_ks/data/bocu/data_df.csv")
write.csv(surv_df, "adler_2007_ks/data/bocu/survival_df.csv")
write.csv(grow_df, "adler_2007_ks/data/bocu/growth_df.csv")
write.csv(recr_df, "adler_2007_ks/data/bocu/recruitment_df.csv")

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

ggsave('adler_2007_ks/results/bocu/overall_hist_sizes_log.png', 
       plot = hist_sizes_log, 
       width = 8, height = 3, units = "in", dpi = 150)

## Survival
n_bins <- 10
h    <- (max(df$logsize_t0, na.rm = T) - min(df$logsize_t0, na.rm = T)) / n_bins
lwr  <- min(df$logsize_t0, na.rm = T) + (h * c(0:(n_bins - 1)))
upr  <- lwr + h
mid  <- lwr + (1/2 * h)

binned_prop <- function(lwr_x, upr_x, response){
  
  id  <- which(df$logsize_t0 > lwr_x & df$logsize_t0 < upr_x) 
  tmp <- df[id,]
  
  if(response == 'prob'){  return(sum (tmp$survives, na.rm = T) / nrow(tmp))}
  if(response == 'n_size'){return(nrow(tmp))}
  
}

# source( 'C:/Users/tn75utid/Desktop/R_templates/plot_binned_prop.R' )
# plot_binned_prop(df, 10, logsize, survives)

# function removes NAs, and provides standard errors, 
source('plot_binned_prop.R')

# y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
# x_binned <- mid
# y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
# 
# surv_binned <- data.frame(xx  = x_binned, 
#                           yy  = y_binned,
#                           nn  = y_n_size) %>% 
#   setNames(c('logsize_t0', 'survives', 'n_size'))
# 
# 
# surv_overall <- 
#   ggplot(data  = surv_binned, aes(x = logsize_t0, y = survives)) +
#   geom_point(alpha = 1,
#              pch   = 16,
#              size  = 1,
#              color = 'red' ) +
#   scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
#   ylim( 0, 1 ) +
#   theme_bw( ) +
#   theme( axis.text = element_text( size = 8 ),
#          title     = element_text( size = 10 ) ) +
#   labs( x = expression( 'log(size)'[t0] ),
#         y = expression( 'Survival to time t1' ) )

surv_overall <- ggplot(data  = plot_binned_prop(df, 10, 
                                                logsize_t0, survives) 
                       ) +
  geom_point(aes(x = logsize_t0, 
                 y = survives ),
             alpha = 1,
             pch   = 16,
             color = 'red' ) +
  geom_errorbar( aes(x = logsize_t0, 
                     ymin = lwr,
                     ymax = upr),
                 size = 0.5,
                 width = 0.5) +
  scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
  ylim( 0, 1 ) +
  theme_bw( ) +
  theme( axis.text = element_text( size = 8 ),
         title     = element_text( size = 10 ) ) +
  labs( x = expression( 'log(size)'[t0] ),
        y = expression( 'Survival to time t1' ) )


ggsave('adler_2007_ks/results/bocu/overall_surv.png', 
       plot = surv_overall, 
       width = 4, height = 3, units = "in", dpi = 150)

## Growth
gr_overall <-
  ggplot(data  = grow_df, aes( x = logsize_t0, y = logsize_t1)) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 0.7,
              color = 'red' ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ) ) +
  labs( x = expression( 'log( size )'[t0] ),
        y = expression( 'log( size )'[t1] ) )

ggsave('adler_2007_ks/results/bocu/overall_gr.png', 
       plot = gr_overall, 
       width = 4, height = 3, units = "in", dpi = 150)

## Recruitment
rec_overall <- 
  ggplot(recr_df, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +
  theme_bw() +
  labs(x = expression('Total parent plant area'[t0]),
       y = expression('Number of recruits'[t1]))

ggsave('adler_2007_ks/results/bocu/overall_rec.png', 
       plot = rec_overall, 
       width = 4, height = 3, units = "in", dpi = 150)



# Fitting vital rate models for the mean IPM -----------------------------------

# Growth
gr_mod_mean <- lm(logsize_t1 ~ logsize_t0, data = grow_df)

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
ggsave('adler_2007_ks/results/bocu/overall_grow_pred.png', 
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
  geom_jitter(data = surv_df, aes(x = logsize_t0, y = survives),
              alpha = 0.25, width = 0, height = 0.25) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, y = survives),
            color = 'red', lwd   = 2 ) +
  theme_bw()

surv_bin <- 
  ggplot() +
  geom_point(data =  plot_binned_prop(df, 10, logsize_t0, survives), 
             aes(x = logsize_t0, 
                 y = survives) ) +
  geom_errorbar(data =  plot_binned_prop(df, 10, logsize_t0, survives), 
             aes(x = logsize_t0, 
                 ymin = lwr,
                 ymax = upr) ) +
  geom_line(data = surv_pred_df, aes(x = logsize_t0, y = survives),
            color = 'red', lwd   = 2) +
  theme_bw()

surv_overall_pred <- surv_line + surv_bin + plot_layout()

ggsave('adler_2007_ks/results/bocu/overall_surv_pred.png', 
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
