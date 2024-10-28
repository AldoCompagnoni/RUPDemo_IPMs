# Example of how to construct an MPM from chart quadrat data

# Author: Aldo Compagnoni
# Email: aldo.compagnoni@gmail.com
# Date: 2024.10.23

# Original models from https://doi.org/10.1111/j.1365-2745.2009.01585.x
#   check the appendices for details on model fitting!

# Packages ---------------------------------------------------------------------

# read and read 
source( 'helper_functions/load_packages.R' )
load_packages( tidyverse, patchwork, skimr, lme4, ggthemes )


# Data -------------------------------------------------------------------------

# Define the species variable
species <- 'Psoralea tenuiflora'
sp_abb  <- 'pste'

# Species data frame
df <- read.csv( 'adler_2007_ks/data/pste/ks_pste.csv' ) %>% 
        filter(Species == species) %>%
        select( -c(Suspect, nearEdge, Site) ) %>% # geometry, 
        mutate(across(c(Quad), as.factor)) %>% 
        rename(species  = Species, 
               survives = survives_tplus1,
               year     = Year,
               quad     = Quad,
               track_id = trackID)

# first look at the data
skim(df)


# Data exploration -------------------------------------------------------------

# Quadrat inventory

# grouping by quard, year and counting the individuals
inv_plot_per_year <- df %>% 
  group_by(quad, year) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>% 
  # ploting year against quadrat
  ggplot() + 
  geom_point(aes(x = year, 
                 y = quad)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 5))

inv_plot_per_year

ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, 
              "/overall_inventory_quadrat_per_year.png"),  
       plot = inv_plot_per_year,
       width = 6, height = 4, dpi = 150)


# Getting the dfs ready --------------------------------------------------------

# Survival data that retains all ages
surv_out_df <- df %>% 
                subset( !is.na(survives) ) %>%
                select(quad, track_id, year, age, survives ) 

# Total number of individuals time t0
indiv_t0  <- df %>%  
              count(quad, year) %>% 
              mutate( year = as.integer(year) ) %>% 
              mutate( year = year + 1 ) %>% 
              rename( parents = n )

# cases for quadrats and year
cases_df  <- df %>% 
              dplyr::select( quad, year ) %>% 
              mutate( year = as.integer(year) ) %>% 
              unique

# Total number of recruits at time t1
indiv_t1 <- df %>%  
              subset( recruit == 1 ) %>% 
              count(quad, year) %>% 
              mutate( year = as.integer(year) ) %>% 
              rename( recruits = n ) %>% 
              # include all available quad X year combinations
              right_join( cases_df ) %>% 
              # if is.na(recruits), then recruits == 0
              mutate( recruits = replace( recruits,
                                          is.na(recruits),
                                          0) )

# recruits data frame
recr_df  <- left_join( indiv_t0, indiv_t1 ) %>% 
              # change this, so that it is year_t0 instead of year_t1
              mutate( year = year - 1 ) %>% 
              mutate( year = as.character(year) ) %>% 
              # calculate per capita recruitment. 
              #   Support: 0 to infinity.
              mutate( pcr = recruits / parents ) 

# Create the data folder for the species if it does not exits already
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))
}

# Save all data
write.csv(df,          paste0("adler_2007_ks/data/", sp_abb, "/data_df.csv"))
write.csv(surv_out_df, paste0("adler_2007_ks/data/", sp_abb, "/survival_df.csv"))
write.csv(recr_df,     paste0("adler_2007_ks/data/", sp_abb, "/recruitment_df.csv"))


# Survival ---------------------------------------------------------------------

# survival data. For this model, only two age classes
surv_df <- surv_out_df %>% 
                # consider all individuals past year 0 as if they had age 1
                mutate( age = dplyr::case_when(
                                age == 0 ~ 0,
                                age > 0  ~ 1 ) 
                      )

# Survival model: Random effect of year only
mod_surv   <- glmer( survives ~ age + (1 | year) + (1 | quad),
                     data = mutate(surv_df, age = as.factor(age)),
                     family = 'binomial' )

# calculate year-specific mean predictions
shat0       <- data.frame( age  = 0,
                           year = coef(mod_surv)$year %>% rownames,
                           stringsAsFactors = F ) %>%
                  mutate( shat = coef(mod_surv)$year[,1] %>% boot::inv.logit() )
shat1       <- data.frame( age  = 1,
                           year = coef(mod_surv)$year %>% rownames,
                           stringsAsFactors = F ) %>%
                  mutate( shat = coef(mod_surv)$year[,1:2] %>% 
                            rowSums %>% 
                            boot::inv.logit() )

# Hmmm, I realize now this is a sub-optimal object name...
shat_df     <- bind_rows( shat0, shat1 )

# plot raw data against model predictions
p_surv <- surv_df %>%
            subset( !is.na(age) ) %>%
            group_by( year, age ) %>%
            summarise( surv_p = sum(survives) / n() ) %>%
            ungroup %>%
            mutate( year = as.character(year) ) %>% 
            left_join( shat_df ) %>% 
            mutate( age = as.factor(age) ) %>% 
            pivot_longer( surv_p:shat, 
                          names_to  = 'Measure',
                          values_to = 'Survival' ) %>% 
            mutate( Measure = dplyr::case_when(
                          Measure == 'shat'   ~ 'Predicted',
                          Measure == 'surv_p' ~ 'Observed' ) 
                    ) %>% 
            ggplot() +
            geom_point( aes(x = age, 
                            y = Survival,
                            color = Measure),
                        alpha = 0.75 ) +
            theme_bw() +
            facet_wrap( ~ year, ncol = 8 ) +
            lims( y = c(0,1) ) +
            labs( x = 'Age',
                  y = 'Survival' ) +
            theme( legend.position = 'top' ) +
            scale_color_colorblind()

# store the plot
ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, 
              "/years_survival.png"),  
       plot = p_surv,
       width = 6, height = 9, dpi = 150)


# Recruitment two-stage model --------------------------------------------------
#   1. Model recruitment being there or not
#   2. Model per-capita recruitment (pcr) *conditional* on recruitment happening

# 1. yes or no recruits?
yesno_df <- recr_df %>% 
              mutate( recruits_yesno = as.numeric(recruits > 0) ) %>% 
              dplyr::select( quad, year, parents, recruits_yesno ) %>% 
              drop_na

# Yes-or-no recruits? (NOTE: parent number matters)
yesno_mod    <- glmer( recruits_yesno ~ (1 | year),
                       data = yesno_df,
                       family = binomial )

# predict zeroes
yesno_pred_df <- expand.grid( year    = coef(yesno_mod)$year %>% rownames,
                              stringsAsFactors = F ) %>% 
                    mutate( yesno_prob = predict( yesno_mod, 
                                                  newdata = .,
                                                  type = 'response' ) 
                    ) 


# 2. per capita recruitment **conditional** on recruitment happening
pcr_df       <- recr_df %>% subset( !(recruits == 0) )

# Per-capita recruitment model
pcr_mod      <- glmer( pcr ~ (1 | year),
                       data = pcr_df,
                       family = Gamma(link = "log") )

# predict per capita recruitment
pcr_pred_df  <- expand.grid( year    = coef(pcr_mod)$year %>% rownames,
                             stringsAsFactors = F ) %>% 
                  # conditional (on recruitment == 1) per capita recruitment
                  mutate( cond_pcr_hat = predict( pcr_mod, 
                                                  newdata = .,
                                                  type = 'response' ) 
                  )

# hurdle model predictions
recr_pred_df <- left_join( yesno_pred_df, 
                           pcr_pred_df ) %>% 
                  mutate( pcr_hat = cond_pcr_hat * yesno_prob ) 

# plot on recruitment
recr_p <- recr_df %>% 
            left_join( recr_pred_df ) %>% 
            mutate( recr_hat = pcr_hat * parents ) %>% 
            subset( year > 38 & year < 72 ) %>% 
            ggplot() + 
            geom_point( aes( x = parents,
                             y = recruits) ) +
            geom_point( aes( x = parents,
                             y = recr_hat),
                        color = 'red' ) +
            geom_abline( intercept = 0,
                         slope     = 1 ) +
            facet_wrap( ~ year, ncol = 4 ) +
            theme_bw()

# store the plot
ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, 
              "/years_recruitment.png"),  
       plot = recr_p,
       width = 6, height = 9, dpi = 150)


# compare "naive" per capita recruitment, with the model's
recr_df %>% 
  group_by( year ) %>% 
  summarise( pcr_naive = mean(pcr,na.rm=T) ) %>% 
  ungroup %>% 
  left_join(recr_pred_df ) %>% 
  ggplot() +
  geom_point( aes(  x = pcr_naive,
                    y = pcr_hat) ) +
  geom_abline( intercept = 0,
               slope     = 1 )


# Matrix population model ------------------------------------------------------

# MPM parameters
mpm_df <- shat_df %>% 
              pivot_wider( names_from  = 'age',
                           values_from = 'shat' ) %>% 
              left_join( dplyr::select( recr_pred_df,
                                        year, pcr_hat) ) %>% 
              drop_na %>% 
              ungroup

# lambda function
lambda_mpm <- function(x) x %>% eigen %>% .$value %>% Re %>% max

# year specific lambdas. 
#   Meant to work using rows of a data frame
year_lam <- function( year_x, mat_df ){
  
  tmp_pars <- mat_df %>% subset( year == year_x ) 
  
  mat      <- matrix(NA,2,2)
  
  # per capita recruitment values
  mat[1,]  <- tmp_pars$pcr_hat
  
  # survival values
  mat[2,1] <- tmp_pars$`0`
  mat[2,2] <- tmp_pars$`1`
  
  mat %>% lambda_mpm
  
}

# Make matrix. 
#   based on parameters in rows of a data frame
make_mat <- function( year_x, mat_df ){
  
  tmp_pars <- mat_df %>% subset( year == year_x ) 
  
  mat      <- matrix(NA,2,2)
  
  # per capita recruitment values
  mat[1,]  <- tmp_pars$pcr_hat
  
  # survival values
  mat[2,1] <- tmp_pars$`0`
  mat[2,2] <- tmp_pars$`1`
  
  # output the matrix
  mat
  
}

# project based on stage distribution
stage_counts <- df %>%
                  mutate( year = as.character(year) ) %>% 
                  mutate( age = dplyr::case_when(
                          age == 0 ~ 0,
                          age > 0  ~ 1 ) 
                          ) %>% 
                  group_by(year, age ) %>%
                  summarize(n_t0 = n()) %>% 
                  ungroup %>% 
                  # only retain year_t0 for which we have data
                  right_join( dplyr::select(mpm_df,year) ) %>% 
                  drop_na %>% 
                  pivot_wider( names_from  = 'age',
                               values_from = 'n_t0',
                               values_fill = 0 ) %>% 
                  drop_na %>% 
                  mutate( year = as.character(year) )
  
# projected lambda (takes into account observed stage distribution)
proj_lam <- function( year_x, mat_df ){
  
  # year-specific stage counts
  tmp_stage <- stage_counts %>% 
                subset( year == year_x ) %>% 
                .[,c('0','1')]
  
  # year-specific projection matrix
  tmp_mat   <- make_mat( year_x, mat_df )
  
  nt0 <- sum(tmp_stage)
  nt1 <- ( tmp_mat %*% as.numeric(tmp_stage) ) %>% sum

  nt1 / nt0
  
}

# model lambdas
lam_df        <- mpm_df %>% 
                  mutate( lam      = sapply( year, year_lam, mpm_df),
                          proj_lam = sapply( year, proj_lam, mpm_df) )
                
# Population counts at time t0
pop_counts_t0 <- df %>%
                  group_by(year, quad) %>%
                  summarize(n_t0 = n()) %>% 
                  ungroup %>% 
                  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
                  group_by(year, quad ) %>%
                  summarize(n_t1 = n()) %>% 
                  ungroup 

# Calculate observed population growth rates, 
#   accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
                # by dropping NAs, we remove gaps in sampling!
                drop_na %>% 
                group_by(year) %>% 
                summarise( n_t0 = sum(n_t0),
                           n_t1 = sum(n_t1) ) %>% 
                ungroup %>% 
                mutate( obs_pgr = n_t1 / n_t0) %>% 
                mutate( year = as.character(year) ) %>% 
                left_join( mpm_df ) %>% 
                mutate( year = as.character(year) ) %>% 
                left_join( lam_df )

# mostly over-predicted population growth
#   not as bad as the published data though (see below)!
ggplot( pop_counts ) +
  geom_point( aes( obs_pgr, lam ) ) +
  geom_point( aes( obs_pgr, proj_lam ),
              col = 'red' ) +
  geom_smooth( aes( obs_pgr, proj_lam ), method = 'lm') +
  geom_abline( intercept = 0,
               slope     = 1 )


# for now, store everything as an R object that resembles 
mpm_df %>% 
  mutate( year = paste0('19',year) ) %>% 
  mutate( SpeciesAuthor = 'Psoralea tenuiflora' ) %>% 
  setNames( c('year', 
              'surv_age0', 'surv_age1', 
              'per_capita_recruitment',
              'SpeciesAuthor') ) %>% 
  write.csv( paste0("adler_2007_ks/results/", 
                     sp_abb, 
                     "/mpm_pste.csv"),  
             row.names = F )


# compare with COMPADRE's estimates --------------------------------------------

# load latest version of COMPADRE 
#   download here: https://compadre-db.org/Data/Compadre
load( 'adler_2007_ks/data/COMPADRE_v.6.23.5.0.RData' )

# identify the matrices across COMPADRE 
id <- compadre$metadata %>% 
        subset( grepl('Psoralea', SpeciesAuthor ) ) %>% 
        subset( MatrixComposite == 'Individual' ) %>% 
        subset( MatrixStartYear %in% paste0('19', pop_counts$year) ) %>% 
        row.names %>% 
        as.numeric

# years we can directly compare
yrs_mat   <- compadre$metadata$MatrixStartYear[id] %>% 
              gsub('19','',.)

# compadre's matrices
comp_mats <- compadre$mat[id] %>% setNames( yrs_mat )

# published matrices
pub_mat_df   <- data.frame( year = yrs_mat ) %>% 
  # populate data frame with the three matrix elements
  mutate( `0`     = sapply( yrs_mat, function(x) pluck(comp_mats, x)$matA[2,1] ),
          `1`     = sapply( yrs_mat, function(x) pluck(comp_mats, x)$matA[2,2] ),
          pcr_hat = sapply( yrs_mat, function(x) pluck(comp_mats, x)$matA[1,1] )
          ) %>% 
  # calculate lambda and projected lambda
  mutate( lam      = sapply( year, year_lam, .),
          proj_lam = sapply( year, proj_lam, .) ) %>% 
  # re-set names to show these are from the publication
  setNames( c('year',    'pub_0',   'pub_1', 
              'pub_pcr_hat', 'pub_lam', 'pub_proj_lam') ) %>% 
  drop_na

# compare population growth rates across estimations
compare_df <- pop_counts %>% 
                left_join( pub_mat_df ) %>% 
                drop_na

# Recruitment: our estimation vs. publication
ggplot( compare_df ) +
  geom_point( aes( pcr_hat, pub_pcr_hat) ) +
  geom_abline( intercept = 0,
               slope     = 1 )

# Survival age 0: our estimation vs. publication
ggplot( compare_df ) +
  geom_point( aes( `0`, pub_0) ) +
  geom_abline( intercept = 0,
               slope     = 1 )

# Survival age 1: our estimation vs. publication
ggplot( compare_df ) +
  geom_point( aes( `1`, pub_1) ) +
  geom_abline( intercept = 0,
               slope     = 1 )

# Compare population growth rates
ggplot( compare_df ) +
  geom_point( aes( lam, pub_lam) ) +
  geom_smooth( aes( lam, pub_lam), method = 'lm' ) +
  geom_abline( intercept = 0,
               slope     = 1 )


# compare population growth rates
ggplot( compare_df ) +
  # geom_point( aes(pcr_hat,pub_pcr_hat) ) +
  geom_point( aes( obs_pgr, pub_lam ), col = 'green' ) +
  geom_point( aes( obs_pgr, pub_proj_lam ), col = 'blue' ) +
  geom_point( aes( obs_pgr, lam ), col = 'darkred' ) +
  geom_point( aes( obs_pgr, proj_lam ), col = 'red' ) +
  geom_abline( intercept = 0,
               slope     = 1 ) +
  geom_smooth( aes( obs_pgr, lam), method = 'lm' )


# calculate room mean squared error
rmse <- function( target, yhat ) (target - yhat)^2 %>% mean %>% sqrt  

# data frame to ease computation
compute_df <- compare_df %>% 
                dplyr::select( obs_pgr, pub_lam, pub_proj_lam, 
                               lam, proj_lam ) %>% 
                as.data.frame %>% 
                drop_na

# dramatically better fit than published lambdas
rmse( compute_df$obs_pgr,
      compute_df$pub_lam ) 
rmse( compute_df$obs_pgr,
      compute_df$lam ) 
rmse( compute_df$obs_pgr,
      compute_df$pub_proj_lam ) 
rmse( compute_df$obs_pgr,
      compute_df$proj_lam ) 

