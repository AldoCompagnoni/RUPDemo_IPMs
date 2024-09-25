# IPM for kansas bocu

# Niklas Neisse
# 2024.09.23

# 

rm( list = ls() )
options( stringsAsFactors = F )
library( tidyverse )
library( ggplot2 )
library( patchwork )
library(skimr)


# Data -------------------------------------------------------------------------

## 
# x and y are locations within the quadrat
# area refers to the individual in cm2
bocu_buf5_dorm1 <- read_csv("data/BOCU_buf5_dorm1.csv")
skim(bocu_buf5_dorm1)
# Dormancy: mean_sd =  0.11_0.31 -> few individuals go into dormancy
  # 1 = dormancy; for unknown time (1 census?)


## growDnoNA.csv
growDnoNA <- read_csv("data/growDnoNA.csv") %>% 
  select(-c(1))
skim(growDnoNA)

## Recruits
recArea <- read_csv("data/recArea.csv") %>% 
  select(-c(1))
skim(recArea)
recSize <- read_csv("data/recSize.csv") %>% 
  select(-c(1))
skim(recSize)

# Data exploration -------------------------------------------------------------
## Quadrat inventory
# grouping by quard, year and counting the individuals
bocu_buf5_dorm1 %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(SP_ID)) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>% 
# ploting year against quadrat
  ggplot() + 
  geom_point(aes(x = year, y = quad))

growDnoNA %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(trackID)) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>% 
  # ploting year against quadrat
  ggplot() + 
  geom_point(aes(x = year, y = quad))


## fixing dormancy
bocu_buf5_dorm1 %>% count(.$year == max(.$year) & is.na(.$survives))
# bocu_buf5_dorm1 %>% 
#   mutate(survives = replace(
#     survives, is.na(survives) & year == max(year), 0))


# Getting the dfs ready --------------------------------------------------------

## Full data set
df <- bocu_buf5_dorm1 %>% full_join(growDnoNA, 
                                    by = c('year', 'quad', 'trackID', 'distEdgeMin',
                                           'age')) %>% 
  mutate(logsize = log(area.t0))


## Survival
surv <- subset(df, !is.na(survives)) %>%
  subset(area != 0) %>%
  select(quad, year, trackID, Group, Grazing,
         area.t0, logsize,
         survives, area.t1)

## Growth
grow <- df %>% 
  subset(area.t0 != 0) %>%
  subset(area.t1 != 0) %>% 
  select(quad, year, trackID, Group, Grazing,
         area.t0, logsize,
         survives, area.t1)

## Total area 
quad_df <- df %>%  
  group_by(quad, year) %>% 
  summarise(tot_p_area = sum(area.t0, na.rm = T)) %>% 
  ungroup

group_df <- quad_df %>% 
  group_by(year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

cover_df <- left_join(quad_df, group_df) %>% 
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

## Number of recruits
recr_df <- recArea

recr <- recr_df %>% full_join(cover_df)

write.csv( df, "data/data_df_df.csv" )
write.csv( surv, "data/survival_df.csv" )
write.csv( grow, "data/growth_df.csv" )
write.csv( recr, "data/recruitment_df.csv" )

# Plotting the data --------------------------------------------------------
# Raw sizes 
png( 'results/histograms_log.png', width = 8, height = 3, units = "in", res = 150 )
par( mfrow = c( 1, 2 ), mar = c( 3.5, 3.5, 1, 0.2 ), mgp = c( 2, 0.7, 0 ), cex = 0.8 )
hist(log(df$area.t0), main = "Histogram of size at time t0", xlab = "Size at time t0")
hist(log(df$area.t1), main = "Histogram of size at time t1", xlab = "Size at time t1" )
dev.off()

## Survival
n_bins <- 10
h    <- (max(df$logsize, na.rm = T) - min(df$logsize, na.rm = T)) / n_bins
lwr  <- min(df$logsize, na.rm = T) + (h * c(0:(n_bins - 1)))
upr  <- lwr + h
mid  <- lwr + (1/2 * h)

binned_prop <- function(lwr_x, upr_x, response){
  
  id  <- which(df$logsize > lwr_x & df$logsize < upr_x) 
  tmp <- df[id,]
  
  if(response == 'prob'){  return(sum (tmp$survives, na.rm = T) / nrow(tmp))}
  if(response == 'n_size'){return(nrow(tmp))}
  
}

# source( 'C:/Users/tn75utid/Desktop/R_templates/plot_binned_prop.R' )
# plot_binned_prop(df, 10, logsize, survives)

df$logsize

y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
x_binned <- mid
y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist

surv_binned <- data.frame(xx  = x_binned, 
                          yy  = y_binned,
                          nn  = y_n_size) %>% 
  setNames(c('logsize', 'survives_tplus1', 'n_size'))


png( 'results/survival_binned.png', width = 6, height = 4, units = "in", res = 150 )
ggplot(data  = surv_binned, aes(x = logsize, y = survives_tplus1)) +
  geom_point(alpha = 1,
             pch   = 16,
             size  = 1,
             color = 'red' ) +
  scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
  
  theme_bw( ) +
  theme( axis.text = element_text( size = 8 ),
         title     = element_text( size = 10 ) ) +
  labs( x = expression( 'log(size)'[t0] ),
        y = expression( 'Survival to time t1' ) )
dev.off()


## Growth
png( 'results/growth.png', width = 6, height = 4, units = "in", res = 150 )

ggplot(data  = grow, aes( x = logsize, y = log(area.t1))) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 0.7,
              color = 'red' ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ) ) +
  labs( x = expression( 'log( size )'[t0] ),
        y = expression( 'log( size )'[t1] ) )
dev.off()

## Recruitment
png('results/recruit.png', width = 6, height = 4, units = "in", res = 150 )

ggplot( recr, aes( x = totParea, y = NRquad ) ) + 
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  theme_bw( ) +
  labs( x = expression( 'Total parent plant area'[t0] ),
        y = expression( 'Number of recruits'[t1] ) )
dev.off()


# Fitting vital rate models for the mean IPM -----------------------------------
## Growth
grow_df      <- grow %>% 
  mutate( logarea.t0 = log( area.t0 ),
          logarea.t1 = log( area.t1 ) )

gr_mod_mean <- lm( logarea.t1 ~ logarea.t0, data = grow_df)

grow_df$pred <- predict( gr_mod_mean, type = "response" )

grow_line <- ggplot( grow_df, aes( x = logarea.t0, y = logarea.t1 ) ) +
  geom_point( ) +
  geom_abline( aes( intercept = coef( gr_mod_mean )[1],
                    slope     = coef( gr_mod_mean )[2] ),
               color     = 'red',
               lwd       = 2 )

grow_pred <- ggplot( grow_df, aes( x = logarea.t1, y = pred ) ) +
  geom_point( ) +
  geom_abline( aes( intercept = 0,
                    slope = 1 ),
               color = "red",
               lwd = 2 )

grow_line + grow_pred + plot_layout( )

x         <- fitted( gr_mod_mean )
y         <- resid( gr_mod_mean )^2
gr_var_m  <- nls( y ~ a * exp( b * x ), start = list( a = 1, b = 0 ) )


## Survival
surv_df      <- surv %>% 
  mutate( logarea = log( area.t0 ),
          logarea_t0_2 = logarea^2,
          logarea_t0_3 = logarea^3,
          )

su_mod_mean <- glm( survives ~ logarea, data = surv_df, family = "binomial" )
su_mod_mean_2 <- glm( survives ~ logarea + logarea^2, data = surv_df, family = "binomial" )
su_mod_mean_3 <- glm( survives ~ logarea + logarea^2 + logarea^3, data = surv_df, family = "binomial" )


surv_x <- seq( min( surv_df$logarea, na.rm = T), max( surv_df$logarea, na.rm = T), length.out = 100)
surv_pred <- boot::inv.logit( coef( su_mod_mean )[1] + coef( su_mod_mean )[2] * surv_x )

surv_pred_df <- data.frame( logarea = surv_x, survives_tplus1 = surv_pred )

surv_line <- ggplot( ) +
  geom_jitter( data = surv_df, aes( x        = logarea, 
                                    y        = survives ), 
               alpha    = 0.25, 
               width    = 0, 
               height   = 0.25 ) +
  geom_line( data = surv_pred_df, aes( x     = logarea,
                                       y     = survives_tplus1 ),
             color = 'red',
             lwd   = 2 )

surv_bin <- ggplot( ) +
  geom_point( data = surv_binned, aes( x     = logsize, 
                                       y     = survives_tplus1 ) ) +
  geom_line( data = surv_pred_df, aes( x     = logarea,
                                       y     = survives_tplus1 ),
             color = 'red',
             lwd   = 2 )

# png( 'results/Bou_gra/survival_pred.png', width = 10, height = 8, units = "in", res = 150 )

surv_line + surv_bin + plot_layout( )


## Recruitment
rec_mod_mean <- MASS::glm.nb( NRquad ~ 1, data = recr_df )

recr_df        <- recr_df %>% 
  mutate( pred_mod_mean = predict( rec_mod_mean, type = "response" ) ) 

rec_sums_df_m <- recr_df %>%
  summarize( NRquad    = sum( NRquad ),
             pred_mod_mean = sum( pred_mod_mean ) )

indiv_m <- surv_df %>%
  summarize( n_adults = n( ) )


repr_pc_m <- indiv_m %>%
  bind_cols( rec_sums_df_m ) %>%
  mutate( repr_pc_mean   = pred_mod_mean / n_adults ) %>%
  mutate( repr_pc_obs    = NRquad / n_adults ) %>%
  drop_na

repr_pc_m
