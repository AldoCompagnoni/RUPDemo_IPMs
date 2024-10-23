# Example of how to construct an MPM from chart quadrat data

# Author: Aldo Compagnoni
# Email: aldo.compagnoni@gmail.com
# Date: 2024.10.23


# Packages ---------------------------------------------------------------------

# read and read 
source( 'helper_functions/load_packages.R' )
load_packages( tidyverse, patchwork, skimr )


# Data -------------------------------------------------------------------------

# Define the species variable
species <- 'Psoralea tenuiflora'
sp_abb  <- 'pste'

## Species data frame
df <- read.csv( 'adler_2007_ks/data/pste/ks_pste.csv' ) %>% 
        filter(Species == species) %>%
        select( -c(Suspect, nearEdge, Site) ) %>% # geometry, 
        mutate(across(c(Quad), as.factor)) %>% 
        rename(species  = Species, 
               survives = survives_tplus1,
               year     = Year,
               quad     = Quad,
               track_id = trackID)

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

# Survival
surv_df <- df %>% 
  subset( !is.na(survives) ) %>%
  select(quad, track_id, year, age, survives ) %>% 
  # consider all individuals past year 0 as if they had age 1
  mutate( age = dplyr::case_when(
                  age == 0 ~ 0,
                  age > 0  ~ 1
                  ) 
          )

# Total number of individuals time t0
indiv_t0 <- df %>%  
  count(quad, year) %>% 
  mutate( year = year + 1 ) %>% 
  mutate( year = as.integer(year) ) %>% 
  rename( parents = n )

# Total number of recruits at time t1
indiv_t1 <- df %>%  
  subset( recruit == 1 ) %>% 
  count(quad, year) %>% 
  mutate( year = as.integer(year) ) %>% 
  rename( recruits = n )

# recruits data frame
recr_df  <- left_join( indiv_t0, indiv_t1 ) %>% 
              # change this, so that it is year_t0 instead of year_t1
              mutate( year = year - 1 ) %>% 
              mutate( pcr = recruits / parents )

# Create the data folder for the species if it does not exits already
if (!dir.exists(paste0("adler_2007_ks/data/", sp_abb))) {
  dir.create(paste0("adler_2007_ks/data/", sp_abb))
}

# Save all data
write.csv(df,      paste0("adler_2007_ks/data/", sp_abb, "/data_df.csv"))
write.csv(surv_df, paste0("adler_2007_ks/data/", sp_abb, "/survival_df.csv"))
write.csv(recr_df, paste0("adler_2007_ks/data/", sp_abb, "/recruitment_df.csv"))


# vital rates ------------------------------------------------------------------

# survival. Random effect of year and quadrat
mod_surv <- glmer( survives ~ age + (1 | year) + (1 | quad),
                   data = mutate(surv_df, age = as.factor(age)),
                   family = 'binomial' )

# find quadrats
quad_v     <- setdiff( surv_df$quad %>% unique, 'e2qa-5' )

# find missing years
year_v     <- intersect( 35:70, coef(mod_surv)$year %>% rownames %>% as.numeric )

# predicted surv rates
shat_v     <- predict( mod_surv, 
                       newdata = expand.grid( age  = as.factor(0:1),
                                              year = year_v,
                                              quad = quad_v),
                       type = 'response') 

# mean predictions
shat_df    <- expand.grid( age  = 0:1,
                           year = year_v,
                           quad = quad_v) %>% 
                mutate( shat = shat_v ) %>% 
                # average across quadrats 
                group_by( year, age ) %>% 
                summarise( shat = mean(shat) ) %>% 
                ungroup
                
# plot raw data against model predictions
p_surv <- surv_df %>%
            subset( !is.na(age) ) %>%
            group_by( year, age ) %>%
            summarise( surv_p = sum(survives) / n() ) %>%
            ungroup %>%
            left_join( shat_df ) %>% 
            mutate( age = as.factor(age) ) %>% 
            ggplot() +
            geom_point( aes(x = age, surv_p),
                        alpha = 0.75 ) +
            geom_point( aes(x = age, shat),
                        color = 'red',
                        alpha = 0.75 ) +
            theme_bw() +
            facet_wrap( ~year, ncol = 3 ) +
            lims( y = c(0,1) ) +
            labs( x = 'Age',
                  y = 'Survival' )

# store the plot
ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, 
              "/years_survival.png"),  
       plot = p_surv,
       width = 6, height = 9, dpi = 150)



# recruitment
recr_mod <- glmer.nb( recruits ~ 0 + parents + (1 | year) + (1 | quad), 
                      data = recr_df )
pcr_mod  <- glmer( pcr ~ (1 | year) + (1 | quad), 
                   data = recr_df, family = Gamma(link = "log") )


# years in the recruitment model
year_v   <- coef(recr_mod)$year %>% rownames

# quadrats in the recruitment model
quad_v   <- coef(recr_mod)$quad %>% rownames

# predicted surv rates
rhat_v   <- predict( recr_mod, 
                     newdata = expand.grid( parents = 1,
                                            year    = year_v,
                                            quad    = quad_v,
                                            stringsAsFactors = F ),
                     type = 'response' ) 

# plot recruit data versus average predictions
p_recr   <- expand.grid( year    = as.double(year_v),
                         quad    = quad_v,
                         stringsAsFactors = F ) %>% 
              left_join( recr_df ) %>% 
              drop_na %>% 
              mutate( rhat = predict( recr_mod, 
                                      newdata = .,
                                      type = 'response' )  ) %>% 
              ggplot() +
              geom_point( aes( parents, recruits ) ) + 
              geom_point( aes( parents, rhat ),
                          col = 'red', alpha = 0.5) + 
              geom_abline( aes( intercept = 0 ,
                                slope    = 1) ) +
                          facet_wrap( ~ year, ncol = 4 ) +
              theme_bw()
  
# plot recruit data versus average predictions
p_recr   <- expand.grid( year    = as.double(year_v),
                         quad    = quad_v,
                         stringsAsFactors = F ) %>% 
  left_join( recr_df ) %>% 
  drop_na %>% 
  mutate( pcrhat = predict( pcr_mod, 
                          newdata = .,
                          type = 'response' ) ) %>% 
  mutate( recruits_hat = pcrhat*parents ) %>% 
  ggplot() +
  geom_point( aes( parents, recruits ) ) + 
  geom_point( aes( parents, recruits_hat ),
              col = 'red', alpha = 0.5) + 
  geom_abline( aes( intercept = 0 ,
                    slope    = 1) ) +
  facet_wrap( ~ year, ncol = 4 ) +
  theme_bw()


# store the plot
ggsave(paste0("adler_2007_ks/results/", 
              sp_abb, 
              "/years_recruitment.png"),  
       plot = p_recr,
       width = 6, height = 9, dpi = 150)

# PER CAPITA recruitment data frame
pcrhat_df <- expand.grid( year    = as.double(year_v),
                          quad    = quad_v,
                          stringsAsFactors = F ) %>% 
               mutate( pcr = predict( pcr_mod, 
                                      newdata = .,
                                      type = 'response' ) ) %>% 
               group_by( year ) %>% 
               summarise( pcr = mean(pcr) ) %>% 
               ungroup

# Matrix population model ------------------------------------------------------

# MPM parameters
mpm_pars <- shat_df %>% 
              pivot_wider( names_from  = 'age',
                           values_from = 'shat' ) %>% 
              left_join( pcrhat_df ) %>% 
              drop_na %>% 
              ungroup

# populate the matrices
populate_mats <- function( ii ){
  mat <- matrix(NA,2,2)
  
  # per capita recruitment values
  mat[1,]  <- mpm_pars[ii,]$pcr
  
  # survival values
  mat[2,1] <- mpm_pars[ii,]$`0`
  mat[2,2] <- mpm_pars[ii,]$`1`
  
  mat
  
}

# list of MPMs 
mpm_l     <- lapply( 1:nrow(mpm_pars), populate_mats )
  
# lambda function
lamda_mpm <- function(x) x %>% eigen %>% .$value %>% Re %>% max

# year-specific lambdas 
mpm_df    <-  mutate(mpm_pars, lam = sapply(mpm_l, lamda_mpm) )

# project based on stage distribution
stage_counts <- df %>%
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
               values_from = 'n_t0' ) %>% 
  drop_na
  
# projected lambda (takes into account stage distribution)
proj_lam <- function( ii ){
  
  nt0 <- sum(stage_counts[ii,c('0','1')])
  nt1 <- (mpm_l[[ii]] %*% as.numeric(stage_counts[ii,c('0','1')]) ) %>% 
          sum

  nt1 / nt0
  
}

# projected lambdas
proj_lam_df <- stage_counts %>% 
  dplyr::select( year ) %>% 
  mutate( proj_lam = sapply(1:27, proj_lam ) )

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
  summarise(n_t0 = sum(n_t0),
            n_t1 = sum(n_t1)) %>% 
  ungroup %>% 
  mutate(obs_pgr = n_t1 / n_t0) %>% 
  left_join( mpm_df ) %>% 
  left_join( proj_lam_df )

ggplot( pop_counts ) +
  geom_point( aes( obs_pgr, lam ) ) +
  geom_point( aes( obs_pgr, proj_lam ),
              col = 'red' ) +
  geom_abline( intercept = 0,
               slope     = 1 )

