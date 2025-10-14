# Visualization of the "misfits" from various buffers - Bouteloua eriopoda

library( tidyverse )
library( sf )
library( binom )
library( patchwork )

# Data -------------------------------------------------------------------------

mf01 <- read.csv("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers_misfits/anderson_2016_az_boer_mistfits_b_01.csv")
mf02 <- read.csv("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers_misfits/anderson_2016_az_boer_mistfits_b_02.csv")
mf03 <- read.csv("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers_misfits/anderson_2016_az_boer_mistfits_b_03.csv")
mf04 <- read.csv("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers_misfits/anderson_2016_az_boer_mistfits_b_04.csv")
mf05 <- read.csv("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers_misfits/anderson_2016_az_boer_mistfits_b_05.csv")

buff01 <- read_rds("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers/anderson16az_boer_raw_01.rds")
buff02 <- read_rds("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers/anderson16az_boer_raw_02.rds")
buff03 <- read_rds("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers/anderson16az_boer_raw_03.rds")
buff04 <- read_rds("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers/anderson16az_boer_raw_04.rds")
buff05 <- read_rds("C:/Users/Aspen Workman/Dropbox/PHD/RUPDemo_IPMs/anderson_2016_az/data/boer/buffers/anderson16az_boer_raw_05.rds")


# Formatting misfits -----------------------------------------------------------

mf01_sum <- mf01 %>% group_by( group ) %>% summarise( group_sums = length( group ) )
mf01_sum[7,1] <- "total"
mf01_sum[7,2] <- sum( mf01_sum[1:6,2] )
mf01_sum$buffer <- 1

mf02_sum <- mf02 %>% group_by( group ) %>% summarise( group_sums = length( group ) )
mf02_sum[6,1] <- "total"
mf02_sum[6,2] <- sum( mf02_sum[1:5,2] )
mf02_sum$buffer <- 2

mf03_sum <- mf03 %>% group_by( group ) %>% summarise( group_sums = length( group ) )
mf03_sum[5,1] <- "total"
mf03_sum[5,2] <- sum( mf03_sum[1:4,2] )
mf03_sum$buffer <- 3

mf04_sum <- mf04 %>% group_by( group ) %>% summarise( group_sums = length( group ) )
mf04_sum[5,1] <- "total"
mf04_sum[5,2] <- sum( mf04_sum[1:4,2] )
mf04_sum$buffer <- 4

mf05_sum <- mf05 %>% group_by( group ) %>% summarise( group_sums = length( group ) )
mf05_sum[7,1] <- "total"
mf05_sum[7,2] <- sum( mf05_sum[1:6,2] )
mf05_sum$buffer <- 5

buff <- mf01_sum %>% bind_rows( mf02_sum ) %>% bind_rows( mf03_sum ) %>% bind_rows( mf04_sum ) %>% bind_rows( mf05_sum )

buff[which(buff$group == "Other..."),'group'] <- "Other"


# Plotting misfits -------------------------------------------------------------

buff %>% filter( group != "total" ) %>% ggplot( aes( x = buffer,
                                                     y = group_sums,
                                                     fill = group ) ) +
  geom_bar( position = "stack", stat = "identity", color = "white" ) +
  geom_bar( data = buff[which( buff$group == "total" ),], position = "dodge",
            stat = "identity", color = "black", fill = NA ) +
  scale_fill_manual( values = c( "lightgreen", "darkgreen", "lightblue", "darkblue", "plum", "orange", "firebrick" ) ) +
  labs( x = "Buffer width",
        y = "Misfits",
        fill = "Category" ) +
  theme_bw()


# Formatting growth data -------------------------------------------------------

buff01_sum <- buff01 %>% st_drop_geometry() %>% group_by( survives_tplus1 ) %>% 
                         summarise( count = length( survives_tplus1 ) )
buff01_sum$survives_tplus1 <- as.character( buff01_sum$survives_tplus1 )
buff01_sum[4,1] <- "misfits"
buff01_sum[4,2] <- mf01_sum[7,2]
buff01_sum[2,2] <- buff01_sum[2,2] - buff01_sum[4,2]
buff01_sum$buffer <- 1

buff02_sum <- buff02 %>% st_drop_geometry() %>% group_by( survives_tplus1 ) %>% 
  summarise( count = length( survives_tplus1 ) )
buff02_sum$survives_tplus1 <- as.character( buff02_sum$survives_tplus1 )
buff02_sum[4,1] <- "misfits"
buff02_sum[4,2] <- mf02_sum[6,2]
buff02_sum[2,2] <- buff02_sum[2,2] - buff02_sum[4,2]
buff02_sum$buffer <- 2

buff03_sum <- buff03 %>% st_drop_geometry() %>% group_by( survives_tplus1 ) %>% 
  summarise( count = length( survives_tplus1 ) )
buff03_sum$survives_tplus1 <- as.character( buff03_sum$survives_tplus1 )
buff03_sum[4,1] <- "misfits"
buff03_sum[4,2] <- mf03_sum[5,2]
buff03_sum[2,2] <- buff03_sum[2,2] - buff03_sum[4,2]
buff03_sum$buffer <- 3

buff04_sum <- buff04 %>% st_drop_geometry() %>% group_by( survives_tplus1 ) %>% 
  summarise( count = length( survives_tplus1 ) )
buff04_sum$survives_tplus1 <- as.character( buff04_sum$survives_tplus1 )
buff04_sum[4,1] <- "misfits"
buff04_sum[4,2] <- mf04_sum[5,2]
buff04_sum[2,2] <- buff04_sum[2,2] - buff04_sum[4,2]
buff04_sum$buffer <- 4

buff05_sum <- buff05 %>% st_drop_geometry() %>% group_by( survives_tplus1 ) %>% 
  summarise( count = length( survives_tplus1 ) )
buff05_sum$survives_tplus1 <- as.character( buff05_sum$survives_tplus1 )
buff05_sum[4,1] <- "misfits"
buff05_sum[4,2] <- mf05_sum[7,2]
buff05_sum[2,2] <- buff05_sum[2,2] - buff05_sum[4,2]
buff05_sum$buffer <- 5

buff_sum <- buff01_sum %>% bind_rows( buff02_sum ) %>% bind_rows( buff03_sum ) %>% bind_rows( buff04_sum ) %>% bind_rows( buff05_sum )


# Plotting survival (number of polygons which survive) -------------------------

buff_sum %>% ggplot( aes( x = buffer,
                          y = count,
                          fill = as.character( survives_tplus1 ) ) ) +
  geom_bar( position = "stack", stat = "identity" ) +
  scale_fill_manual( values = c( "firebrick", "lightgreen", "plum", "gray" ), 
                     labels = c( "Dies", "Survives", "Misfits", "NA" ) ) +
  labs( x = "Buffer width",
        y = "Individuals",
        fill = "Survival" ) +
  theme_bw()


# Quick plots of number of individuals by buffer size

ind_buff <- buff_sum %>% group_by( buffer ) %>% summarise( total = sum( count ) )

mf_buff  <- buff_sum %>% filter( survives_tplus1 == "misfits" )
mf_buff$prop <- mf_buff$count / ind_buff$total

ind_buff %>% ggplot( aes( x = buffer,
                          y = total ) ) +
  geom_point() +
  geom_line() +
  geom_point( data = mf_buff, aes( x = buffer, y = count ), color = "plum" ) +
  geom_line( data = mf_buff, aes( x = buffer, y = count ), color = "plum" ) +
  labs( x = "Buffer size (cm)",
        y = "Number of individual transitions" ) +
  theme_bw()



# Testing vital rate models ----------------------------------------------------

# Formatting

format_data <- function( df ){
  
  out <- df %>% select( -c( suspect, near_edge, site ) ) %>%
                mutate( across( c( quad ), as.factor ) ) %>%
                rename( size_t0  = basal_area_genet,
                        size_t1  = size_tplus1,
                        survives = survives_tplus1,
                        track_id = track_id ) %>%
                mutate( logsize_t0   = log( size_t0 ),
                        logsize_t1   = log( size_t1 ),
                        logsize_t0_2 = logsize_t0^2,
                        logsize_t0_3 = logsize_t0^3 ) %>%
                filter( logsize_t0 > -10.7 | is.na( logsize_t0 ) ) %>% 
                filter( logsize_t1 > -10.7 | is.na( logsize_t1 ) )
  
  out$id_quad_year <- paste( out$track_id, out$quad, out$year, sep = "_" )
  
  return( out )
}

buff01_f <- format_data( buff01 )
buff02_f <- format_data( buff02 )
buff03_f <- format_data( buff03 )
buff04_f <- format_data( buff04 )
buff05_f <- format_data( buff05 )


# First, drop all misfits

buff01_nmf <- buff01_f[-which( buff01_f$id_quad_year %in% mf01$id_quad_year ),]
buff02_nmf <- buff02_f[-which( buff02_f$id_quad_year %in% mf02$id_quad_year ),]
buff03_nmf <- buff03_f[-which( buff03_f$id_quad_year %in% mf03$id_quad_year ),]
buff04_nmf <- buff04_f[-which( buff04_f$id_quad_year %in% mf04$id_quad_year ),]
buff05_nmf <- buff05_f[-which( buff05_f$id_quad_year %in% mf05$id_quad_year ),]

# Now, only drop misfits that aren't A or AM

mfd01 <- mf01[-which(mf01$group %in% c( "A", "AM" ) ),]
mfd02 <- mf02[-which(mf02$group %in% c( "A", "AM" ) ),]
mfd03 <- mf03[-which(mf03$group %in% c( "A", "AM" ) ),]
mfd04 <- mf04[-which(mf04$group %in% c( "A", "AM" ) ),]
mfd05 <- mf05[-which(mf05$group %in% c( "A", "AM" ) ),]

buff01_smf <- buff01_f[-which( buff01_f$id_quad_year %in% mfd01$id_quad_year ),]
buff02_smf <- buff02_f[-which( buff02_f$id_quad_year %in% mfd02$id_quad_year ),]
buff03_smf <- buff03_f[-which( buff03_f$id_quad_year %in% mfd03$id_quad_year ),]
buff04_smf <- buff04_f[-which( buff04_f$id_quad_year %in% mfd04$id_quad_year ),]
buff05_smf <- buff05_f[-which( buff05_f$id_quad_year %in% mfd05$id_quad_year ),]


# Format data for fitting survival and growth models

format_fit <- function( df ){
  
  surv_df <- subset(df, !is.na(survives)) %>%
    subset(size_t0 != 0) %>%
    select(quad, track_id, year, size_t0, survives, size_t1, 
           logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)
  
  grow_df <- df %>% 
    subset(size_t0 != 0) %>%
    subset(size_t1 != 0) %>% 
    select(quad, track_id, year, size_t0, survives, size_t1, 
           logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)
  
  quad_df <- df %>%  
    st_drop_geometry() %>%
    group_by (quad, year) %>% 
    summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
    ungroup
  
  group_df <- quad_df %>%
    group_by (year) %>% 
    summarise(g_cov = mean(tot_p_area)) %>% 
    ungroup
  
  cover_df <- left_join(quad_df, group_df) %>%
    mutate(year = year + 1) %>%
    mutate(year = as.integer(year)) %>% 
    drop_na()
  
  recr_df <- df %>% 
    st_drop_geometry() %>%
    group_by (year, quad) %>%
    summarise(nr_rec = sum(recruit, na.rm = T)) %>% 
    ungroup
  
  recr_df <- full_join(cover_df, recr_df, by = c('quad', 'year'))
  
  out <- list( surv_df, grow_df, recr_df )
  
  return( out )
  
}

buff01_f_all <- format_fit( buff01_f )
buff02_f_all <- format_fit( buff02_f )
buff03_f_all <- format_fit( buff03_f )
buff04_f_all <- format_fit( buff04_f )
buff05_f_all <- format_fit( buff05_f )

buff01_nmf_f <- format_fit( buff01_nmf )
buff02_nmf_f <- format_fit( buff02_nmf )
buff03_nmf_f <- format_fit( buff03_nmf )
buff04_nmf_f <- format_fit( buff04_nmf )
buff05_nmf_f <- format_fit( buff05_nmf )

buff01_smf_f <- format_fit( buff01_smf )
buff02_smf_f <- format_fit( buff02_smf )
buff03_smf_f <- format_fit( buff03_smf )
buff04_smf_f <- format_fit( buff04_smf )
buff05_smf_f <- format_fit( buff05_smf )

# Fit models

fit_models <- function( list_in ){
  
  surv_mod <- glm( survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                   data = list_in[[1]], family = 'binomial')
  
  grow_mod <- lm( logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
                  data = list_in[[2]] )
  
  # recruitment
  recr_nona_nr_rec <- list_in[[3]] %>% filter( !is.na( nr_rec ) )
  rec_mod_mean <- MASS::glm.nb( nr_rec ~ 1, data = recr_nona_nr_rec )
  recr_nona_nr_rec <- recr_nona_nr_rec %>% 
    mutate(pred_mod_mean = predict( rec_mod_mean, type = 'response' ) ) 
  rec_sums_df_m <- recr_nona_nr_rec %>%
    summarize( nr_rec = sum( nr_rec ),
              pred_mod_mean = sum( pred_mod_mean ) )
  indiv_m <- list_in[[1]] %>%
    summarize( n_adults = n( ) )
  repr_pc_m <- indiv_m %>%
    bind_cols( rec_sums_df_m ) %>%
    mutate( repr_pc_mean = pred_mod_mean / n_adults ) %>%
    mutate( repr_pc_obs = nr_rec / n_adults ) %>%
    drop_na 
  
  out <- data.frame(
    coefficient = c( "surv_b0", "surv_b1", "surv_b2", "surv_b3",
                     "grow_b0", "grow_b1", "grow_b2", "grow_b3",
                     "fecu_b0" ),
    value = c( coef( surv_mod )[[1]], coef( surv_mod )[[2]], 
               coef( surv_mod )[[3]], coef( surv_mod )[[4]],
               coef( grow_mod )[[1]], coef( grow_mod )[[2]], 
               coef( grow_mod )[[3]], coef( grow_mod )[[4]],
               repr_pc_m$repr_pc_mean )
  )
  
  return( out )
  
}

buff01_f_coef <- fit_models( buff01_f_all )
buff02_f_coef <- fit_models( buff02_f_all )
buff03_f_coef <- fit_models( buff03_f_all )
buff04_f_coef <- fit_models( buff04_f_all )
buff05_f_coef <- fit_models( buff05_f_all )

buff01_nmf_coef <- fit_models( buff01_nmf_f )
buff02_nmf_coef <- fit_models( buff02_nmf_f )
buff03_nmf_coef <- fit_models( buff03_nmf_f )
buff04_nmf_coef <- fit_models( buff04_nmf_f )
buff05_nmf_coef <- fit_models( buff05_nmf_f )

buff01_smf_coef <- fit_models( buff01_smf_f )
buff02_smf_coef <- fit_models( buff02_smf_f )
buff03_smf_coef <- fit_models( buff03_smf_f )
buff04_smf_coef <- fit_models( buff04_smf_f )
buff05_smf_coef <- fit_models( buff05_smf_f )


# Plotting ---------------------------------------------------------------------

source('helper_functions/predictor_fun.R')
source('helper_functions/plot_binned_prop.R')

# Survival

plot_surv <- function( i, dfs, coefs ){
  
  dat <- dfs[[i]] %>% st_drop_geometry( )
  coef <- as.list( pivot_wider( coefs[[i]], names_from = "coefficient",
                                values_from = "value" ) )
  
  # Generate predictions for survival across a range of sizes
  surv_x <- seq( min( dat$logsize_t0, na.rm = T ), 
                 max( dat$logsize_t0, na.rm = T ), length.out = 100 )
  
  # Prepare data for survival plot
  surv_pred_df <- predictor_fun( surv_x, as.numeric( coef ) ) %>% 
    # Inverse logit for predictions
    boot::inv.logit() %>% 
    data.frame( logsize_t0 = surv_x, survives = . )
  
  # Plot
  out <- ggplot() +
    geom_jitter( data = dat, aes( x = logsize_t0,
                                  y = survives ),
                 alpha = 0.25, width = 0, height = 0.25 ) + 
    geom_line( data = surv_pred_df, aes( x = logsize_t0, 
                                         y = survives ),
              color = "red", 
              lwd   =  2 ) +  
    geom_point( data =  plot_binned_prop( dat, 10, logsize_t0, survives ),
                aes( x = logsize_t0,
                     y = survives ), color = "red" ) +
    geom_errorbar( data =  plot_binned_prop( dat, 10, logsize_t0, survives ),
                   aes( x = logsize_t0,
                        ymin = lwr,
                        ymax = upr) ) +
    xlim( -11, -2.5 ) +
    theme_bw()
  
  return( out )
  
}


raw_data <- list( buff01_f_all[[1]], buff01_smf_f[[1]], buff01_nmf_f[[1]],
                  buff02_f_all[[1]], buff02_smf_f[[1]], buff02_nmf_f[[1]],
                  buff03_f_all[[1]], buff03_smf_f[[1]], buff03_nmf_f[[1]],
                  buff04_f_all[[1]], buff04_smf_f[[1]], buff04_nmf_f[[1]],
                  buff05_f_all[[1]], buff05_smf_f[[1]], buff05_nmf_f[[1]] )

model_coefs <- list( buff01_f_coef[1:4,], buff01_smf_coef[1:4,], buff01_nmf_coef[1:4,],
                     buff02_f_coef[1:4,], buff02_smf_coef[1:4,], buff02_nmf_coef[1:4,],
                     buff03_f_coef[1:4,], buff03_smf_coef[1:4,], buff03_nmf_coef[1:4,],
                     buff04_f_coef[1:4,], buff04_smf_coef[1:4,], buff04_nmf_coef[1:4,],
                     buff05_f_coef[1:4,], buff05_smf_coef[1:4,], buff05_nmf_coef[1:4,] )

surv_plots <- lapply( 1:15, plot_surv, dfs = raw_data, coefs = model_coefs )

layout <- "
ABC
DEF
GHI
JKL
MNO
"

plots_text <- wrap_elements( grid::textGrob( 'Buffer size' ) ) +
  wrap_elements( grid::textGrob( 'All points' ) ) +
  wrap_elements( grid::textGrob( 'Problematic misfits removed' ) ) +
  wrap_elements( grid::textGrob( 'All misfits removed' ) ) +
  plot_layout( ncol = 4, nrow = 1, widths = c( 1, 2, 2, 2 ) )

buffer_text <- wrap_elements( grid::textGrob( '1 cm' ) ) +
  wrap_elements( grid::textGrob( '2 cm' ) ) +
  wrap_elements( grid::textGrob( '3 cm' ) ) +
  wrap_elements( grid::textGrob( '4 cm' ) ) +
  wrap_elements( grid::textGrob( '5 cm' ) ) +
  plot_layout( ncol = 1, nrow = 5 )

surv_plots_plotted <- wrap_plots( surv_plots ) + plot_layout( design = layout,
                                                              axis_titles = "collect" )

surv_plots_2 <- wrap_elements( buffer_text ) + wrap_elements( surv_plots_plotted ) +
  plot_layout( ncol = 2, widths = c( 1, 6 ) )

surv_plots_all <- wrap_elements( plots_text ) + wrap_elements( surv_plots_2 ) +
  plot_layout( nrow = 2, heights = c( 1, 5 ) )

surv_plots_all


# Growth

plot_grow <- function( i, dfs, coefs ){
  
  dat <- dfs[[i]] %>% st_drop_geometry( )
  coef <- as.list( pivot_wider( coefs[[i]], names_from = "coefficient",
                                values_from = "value" ) )
  
  # Generate predictions for growth across a range of sizes
  grow_x <- seq( min( dat$logsize_t0, na.rm = T ), 
                 max( dat$logsize_t0, na.rm = T ), length.out = 100 )
  
  # Prepare data for survival plot
  grow_pred_df <- predictor_fun( grow_x, as.numeric( coef ) ) %>% 
    data.frame( logsize_t0 = grow_x, logsize_t1 = . )
  
  # Plot
  out <- ggplot() +
    geom_point( data = dat, aes( x = logsize_t0,
                                  y = logsize_t1 ),
                alpha = 0.25 ) + 
    geom_line( data = grow_pred_df, aes( x = logsize_t0, 
                                         y = logsize_t1 ),
               color = "green", 
               lwd   =  2 ) +  
    xlim( -11, -2.5 ) +
    ylim( -11, -2.5 ) +
    theme_bw()
  
  return( out )
  
}

raw_dat2 <- list( buff01_f_all[[2]], buff01_smf_f[[2]], buff01_nmf_f[[2]],
                  buff02_f_all[[2]], buff02_smf_f[[2]], buff02_nmf_f[[2]],
                  buff03_f_all[[2]], buff03_smf_f[[2]], buff03_nmf_f[[2]],
                  buff04_f_all[[2]], buff04_smf_f[[2]], buff04_nmf_f[[2]],
                  buff05_f_all[[2]], buff05_smf_f[[2]], buff05_nmf_f[[2]] )

model_coef2 <- list( buff01_f_coef[5:8,], buff01_smf_coef[5:8,], buff01_nmf_coef[5:8,],
                     buff02_f_coef[5:8,], buff02_smf_coef[5:8,], buff02_nmf_coef[5:8,],
                     buff03_f_coef[5:8,], buff03_smf_coef[5:8,], buff03_nmf_coef[5:8,],
                     buff04_f_coef[5:8,], buff04_smf_coef[5:8,], buff04_nmf_coef[5:8,],
                     buff05_f_coef[5:8,], buff05_smf_coef[5:8,], buff05_nmf_coef[5:8,] )

grow_plots <- lapply( 1:15, plot_grow, dfs = raw_dat2, coefs = model_coef2 )

grow_plots_plotted <- wrap_plots( grow_plots ) + plot_layout( design = layout,
                                                              axis_titles = "collect" )

grow_plots_2 <- wrap_elements( buffer_text ) + wrap_elements( grow_plots_plotted ) +
  plot_layout( ncol = 2, widths = c( 1, 6 ) )

grow_plots_all <- wrap_elements( plots_text ) + wrap_elements( grow_plots_2 ) +
  plot_layout( nrow = 2, heights = c( 1, 5 ) )

grow_plots_all


recr_coefs <- list( buff01_f_coef[9,], buff01_smf_coef[9,], buff01_nmf_coef[9,],
                     buff02_f_coef[9,], buff02_smf_coef[9,], buff02_nmf_coef[9,],
                     buff03_f_coef[9,], buff03_smf_coef[9,], buff03_nmf_coef[9,],
                     buff04_f_coef[9,], buff04_smf_coef[9,], buff04_nmf_coef[9,],
                     buff05_f_coef[9,], buff05_smf_coef[9,], buff05_nmf_coef[9,] )
