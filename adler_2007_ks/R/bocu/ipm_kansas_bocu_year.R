# IPM for kansas bocu year specific

# Niklas Neisse
# 2024.09.23

# 

rm( list = ls() )
options( stringsAsFactors = F )
library( tidyverse )
library( ggplot2 )
library( patchwork )
library(skimr)
library(lme4)
library(bbmle)


grow <- read.csv( "data/growth_df.csv" )
surv <- read.csv( "data/survival_df.csv" )
recr <- read.csv( "data/recruitment_df.csv" )
df <- read.csv( "data/data_df_df.csv" )


df$logsize_t1 <- log( df$area.t1)

df_long <- pivot_longer( df, cols = c( logsize, logsize_t1 ), 
                         names_to = "size", values_to = "size_value" ) %>% 
  select(c(area.t0, year, size,  size_value))

df_long$size     <- as.factor( df_long$size)
df_long$Year_fac <- as.factor( df_long$year )



size_labs <- c( "at time t0", "at time t1" )
names(size_labs) <- c( "logsize", "logsize_t1" )

png( 'results/yr/histograms.png', width = 6, height = 13, units = "in", res = 150 )
df_long %>% ggplot( aes( x = size_value ) ) +
  geom_histogram( binwidth = 1 ) +
  facet_grid( Year_fac ~ size, 
              scales = "free_y",
              labeller = labeller( size = size_labs ) ) +
  labs( x = "log( size )",
        y = "Frequency" )
dev.off()


## survival 

source( 'C:/Users/tn75utid/Desktop/R_templates/df_binned_prop_year.R' )

surv_yrs       <- data.frame( year = surv$year %>% unique %>% sort )
surv_bin_yrs   <- lapply( 1:nrow(surv_yrs), df_binned_prop_year, df, 15, 
                          logsize, survives, surv_yrs )

surv_bin_yrs <- Filter(function(df) nrow(df) > 0, surv_bin_yrs)

surv_yr_pan_df <- bind_rows( surv_bin_yrs ) %>% 
  mutate( transition = paste( paste0( year ),
                              substr( paste0( year + 1 ), 3, 4 ),
                              sep = '-' ) ) %>% 
  mutate( year       = as.integer( year - surv_yrs[1,] ) )


png( 'results/yr/survival_binned_yr.png', width = 10, height = 6, units = "in", res = 150 )
ggplot( data   = surv_yr_pan_df, 
        aes( x = logsize, 
             y = survives ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
  # split in panels
  facet_wrap( .~ transition, nrow = 4 ) +
  theme_bw( ) +
  theme( axis.text = element_text( size = 8 ),
         title     = element_text( size = 10 ),
         strip.text.y  = element_text( size   = 5,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.text.x  = element_text( size   = 5,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing         = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'log(size)'[t0] ),
        y = expression( 'Survival to time t1' ) )
dev.off()

## Growth

grow_yr_pan_df <- grow %>%
  mutate( transition = paste( paste0( year ),
                              substr( paste0( year + 1 ), 3, 4 ),
                              sep = '-' ) ) %>% 
  mutate( year       = as.integer( year - 1996 ) )

png( 'results/yr/growth_yr.png', width = 10, height = 6, units = "in", res = 150 )
ggplot(data  = grow_yr_pan_df, aes( x = logsize, y = log( area.t1 ) ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 0.7,
              color = 'red' ) +
  # split in panels
  facet_wrap( .~ transition, nrow = 4 ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ),
         strip.text.y  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.text.x  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing         = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'log( size )'[t0] ),
        y = expression( 'log( size )'[t1] ) )
dev.off()


## Recruits
indiv_qd <- surv %>%
  group_by( quad ) %>%
  count( year ) %>% 
  rename( n_adults = n ) %>% 
  mutate( year = year + 1 )

repr_yr <- indiv_qd %>% 
  left_join( recr ) %>%
  mutate( repr_pc    = NRquad / n_adults ) %>% 
  mutate( year = year - 1 ) %>% 
  drop_na

png( 'results/yr/recruit_yr.png', width = 10, height = 6, units = "in", res = 150 )
repr_yr %>% 
  filter( NRquad != max( repr_yr$NRquad ) ) %>% 
  filter( n_adults != max( repr_yr$n_adults ) ) %>% 
  ggplot( aes( x = n_adults, y = NRquad ) ) +
  geom_point( alpha = 1,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  facet_wrap( .~ year, nrow = 4 ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ),
         strip.text.y  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.text.x  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing         = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'Number of adults '[ t0] ),
        y = expression( 'Number of recruits '[ t1] ) )
dev.off()


## Recruitment size

# !!!
  # This is impossible, we dont have the data for the recruit
# !!!

# recSize <- df %>% subset( recruit == 1)
# 
# recSize$Year_fac <- as.factor( recSize$Year )
# 
# # png( 'results/Bou_gra_yr/recr_histograms.png', width = 10, height = 6, units = "in", res = 150 )
# 
# recSize %>% ggplot( aes( x = logsize ) ) +
#   geom_histogram( ) +
#   facet_wrap( Year_fac ~ ., 
#               scales = "free_y",
#               nrow = 4 ) +
#   labs( x = expression('log( size )'[t0]),
#         y = "Frequency" )


# Fitting vital rate models with the random effect of year ---------------------
## Survival
surv_df      <- surv %>% 
  mutate( logarea = log( area.t0 ),
          logarea_2 = logarea^2,
          logarea_3 = logarea^3) %>% 
  drop_na()

su_mod_yr <- glmer( survives ~ logarea + ( logarea | year ), data = surv_df, family = binomial )
ranef_su <- data.frame( coef( su_mod_yr )[1] )


su_mod_yr_2 <- glmer( survives ~ logarea + logarea_2 + ( logarea | year ), data = surv_df, family = binomial )
  # Model failed to converge with max|grad| = 0.00247821
ranef_su_2 <- data.frame( coef( su_mod_yr_2 )[1] )

su_mod_yr_3 <- glmer( survives ~ logarea + logarea_2 + logarea_3 + ( logarea | year ), data = surv_df, family = binomial )
  # Model is nearly unidentifiable: very large eigenvalue
ranef_su_3 <- data.frame( coef( su_mod_yr_3 )[1] )

s_mods <- c( su_mod_yr, su_mod_yr_2, su_mod_yr_3 )
AICtab( s_mods, weights = T )
  # model2: incl. quadratic term has the lowest AIC


years_v  <- c( surv_bin_yrs[[1]][1,'year']:surv_bin_yrs[[length(surv_bin_yrs)]][1,'year'] )

v <- vector(rep(NA,length(surv_bin_yrs)))
for (ii in 1:length(surv_bin_yrs)) {
  v[ii] <- surv_bin_yrs[[ii]][1,'year']
}


surv_yr_plots_2 <- function( i ){
  surv_temp   <- as.data.frame( surv_bin_yrs[[i]] )
  x_temp      <- seq( min( surv_temp$logsize, na.rm = T ), 
                      max( surv_temp$logsize, na.rm = T ), 
                      length.out = 100)
  pred_temp   <- boot::inv.logit( ranef_su_2[i,1] + ranef_su_2[i,2] * x_temp ) 
  pred_temp_df <- data.frame( logarea = x_temp, survives = pred_temp )
  temp_plot <- surv_temp %>% ggplot( ) +
    geom_point( aes( x = logsize, y = survives ) ) +
    geom_line( data = pred_temp_df, aes( x     = logarea,
                                         y     = survives ),
               color = 'green',
               lwd   = 1  ) +
    labs( title = paste0( years_v[i] ),
          x = expression( 'log( size )'[t0] ),
          y = expression( 'Survival probability  '[ t1] ) )
  if( i %in% c(setdiff(1:length(years_v), seq(1,length(years_v), by = 4))) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
surv_yrs_2 <- lapply( 1:length(surv_bin_yrs), surv_yr_plots_2)
surv_years_2 <- wrap_plots( surv_yrs_2 ) + plot_layout( nrow = 4 )

png( 'results/yr/survival_pred_2.png', width = 10, height = 8, units = "in", res = 150 )
surv_years_2
dev.off()


## Growth
grow_df <- grow %>% 
  mutate( logarea_t0 = log( area.t0 ),
          logarea_t1 = log( area.t1 ),
          logarea_t0_2 = logarea_t0^2,
          logarea_t0_3 = logarea_t0^3) %>% 
  select(- c(X,logsize)) %>% 
  # lmer does not converge, therefor we reduce the dataset to years with loads of data
  filter(year <= 55)

gr_mod_yr <- lmer( logarea_t1 ~ logarea_t0 + ( logarea_t0 | year ), data = grow_df )
gr_mod_yr_2 <- lmer( logarea_t1 ~ logarea_t0 + logarea_t0_2 + ( logarea_t0 | year ), data = grow_df )
gr_mod_yr_3 <- lmer( logarea_t1 ~ logarea_t0 + logarea_t0_2 + logarea_t0_3 + ( logarea_t0 | year ), data = grow_df )

g_mods <- c( gr_mod_yr, gr_mod_yr_2, gr_mod_yr_3 )
AICtab( g_mods, weights = T )
# model2: incl. quadratic term has the lowest AIC, with weight = 0.99 

ranef_gr <- data.frame( coef( gr_mod_yr )[1] )
ranef_gr_2 <- data.frame( coef( gr_mod_yr_2 )[1] )
ranef_gr_3 <- data.frame( coef( gr_mod_yr_3 )[1] )

grow_yr_plots <- function( i ){
  temp_plot <- grow_df %>% filter( year == i ) %>% ggplot( ) +
    geom_point( aes( x = logarea_t0, 
                     y = logarea_t1 ) ) +
    geom_abline( aes( intercept = ranef_gr[which(rownames( ranef_gr_2 ) == i ),1],
                      slope     = ranef_gr[which(rownames( ranef_gr_2 ) == i ),2] ),
                 color = "red",
                 lwd   = 1 ) +
    labs( title = paste0( i ),
          x = expression( 'log( size ) '[ t0] ),
          y = expression( 'log( size ) '[ t1] ) )
  if( i %in% c(c(setdiff(1:length(unique(grow_df$year)), seq(1,length(unique(grow_df$year)), by = 4)))) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
grow_yrs <- lapply(unique(grow_df$year), grow_yr_plots )
grow_years <- wrap_plots( grow_yrs ) + plot_layout( nrow = 4 )

png( 'results/yr/grow_pred_2.png', width = 10, height = 6, units = "in", res = 150 )
grow_years
dev.off()




