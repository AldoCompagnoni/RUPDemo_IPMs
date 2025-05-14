# MPM year specific Adler 2007 Kansas for any forb species

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.30

# Original models from https://doi.org/10.1111/j.1365-2745.2009.01585.x
#   check the appendices for details on model fitting!


# Packages ---------------------------------------------------------------------

# read and read 
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork, skimr, lme4, ggthemes)


# Data -------------------------------------------------------------------------
# Type of population model
mod_type <- 'mpm'
# Species abbreviation
sp_abb  <- tolower(
  gsub(' ', '', paste(substr(unlist(strsplit(species, ' ')), 1, 2), 
                      collapse = '')))

# Directories 
data0_dir  <- 'adler_2007_ks/data/'
data_dir   <- paste0('adler_2007_ks/data/',    sp_abb, '_', mod_type)
result_dir <- paste0('adler_2007_ks/results/', sp_abb, '_', mod_type)
R_dir      <- paste0('adler_2007_ks/R/',       sp_abb, '_', mod_type)


# Plant tracker if its not already exists
if (!file.exists(paste0(data_dir, '/ks_', sp_abb, '.csv'))) {
  source(paste0(R_dir, '/', sp_abb, '_tracker.R'))
}

# Species data frame
df <- read.csv(paste0(data_dir, '/ks_', sp_abb, '.csv')) %>% 
  filter(Species == species) %>%
  select(-c(Suspect, nearEdge, Site)) %>% # geometry, 
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
  labs(title    = 'Sampling inventory',
       subtitle = species) +
  theme(axis.text.y = element_text(size = 5))

# Create output directory if it doesn't exist
if (!dir.exists(paste0(result_dir))) {
  dir.create(paste0(result_dir))
}

ggsave(paste0(result_dir, 
              '/0_inventory_quadrat_per_year.png'),  
       plot = inv_plot_per_year,
       width = 6, height = 4, dpi = 150)


# Getting the dfs ready --------------------------------------------------------

# Survival data that retains all ages
surv_out_df <- df %>% 
  subset(!is.na(survives)) %>%
  select(quad, track_id, year, age, survives) 

# Total number of individuals time t0
indiv_t0  <- df %>%  
  count(quad, year) %>% 
  mutate(year = as.integer(year)) %>% 
  mutate(year = year + 1) %>% 
  rename(parents = n)

# cases for quadrats and year
cases_df  <- df %>% 
  dplyr::select(quad, year) %>% 
  mutate(year = as.integer(year)) %>% 
  unique

# Total number of recruits at time t1
indiv_t1 <- df %>%  
  subset(recruit == 1) %>% 
  count(quad, year) %>% 
  mutate(year = as.integer(year)) %>% 
  rename(recruits = n) %>% 
  # include all available quad X year combinations
  right_join(cases_df) %>% 
  # if is.na(recruits), then recruits == 0
  mutate(recruits = replace(recruits,
                              is.na(recruits),
                              0))

# recruits data frame
recr_df  <- left_join(indiv_t0, indiv_t1) %>% 
  # change this, so that it is year_t0 instead of year_t1
  mutate(year = year - 1) %>% 
  mutate(year = as.character(year)) %>% 
  # calculate per capita recruitment. 
  #   Support: 0 to infinity.
  mutate(pcr = recruits / parents) 

# Create the data folder for the species if it does not exits already
if (!dir.exists(paste0(data_dir))) {
  dir.create(paste0(data_dir))
}

# Save all data
write.csv(df, row.names = F,          
          paste0(data_dir, '/data_df.csv'))
write.csv(surv_out_df, row.names = F, 
          paste0(data_dir, '/survival_df.csv'))
write.csv(recr_df, row.names = F,     
          paste0(data_dir, '/recruitment_df.csv'))


# Survival ---------------------------------------------------------------------

# survival data. For this model, only two age classes
surv_df <- surv_out_df %>% 
  # consider all individuals past year 0 as if they had age 1
  mutate(age = dplyr::case_when(
    age == 0 ~ 0,
    age > 0  ~ 1))

# Survival model: Random effect of year only
mod_surv   <- glmer(survives ~ age + (1 | year) + (1 | quad),
                     data = mutate(surv_df, age = as.factor(age)),
                     family = 'binomial')

# calculate year-specific mean predictions
surv_ysmpred_re_fun <- function(mod, re) {
  shat_df       <- expand.grid(age  = as.factor(c(0,1)),
                               year = coef(mod)$year %>% rownames,
                               quad = coef(mod)$quad %>% rownames,
                               stringsAsFactors = F) %>%
    mutate(shat = predict(mod,
                          newdata = .,
                          type    = 'response')) %>% 
    group_by( age, !!sym(re)) %>% 
    summarise(shat = mean(shat)) %>% 
    ungroup %>% 
    mutate(age = as.character(age) %>% as.numeric)
  return(shat_df)
}

shat_quad_df <- surv_ysmpred_re_fun(mod = mod_surv, re = 'quad')
shat_df <- surv_ysmpred_re_fun(mod = mod_surv, re = 'year')

# plot raw data against model predictions
su_fun <- function(df, re){
  p_surv <- surv_df %>% 
    subset(!is.na(age)) %>% 
    group_by(!!sym(re), age) %>%
    summarise(surv_p = sum(survives) / n()) %>%
    ungroup %>%
    mutate(!!sym(re) := as.character(!!sym(re))) %>% 
    left_join(df) %>% 
    mutate(age = as.factor(age)) %>% 
    pivot_longer(surv_p:shat, 
                 names_to  = 'Measure',
                 values_to = 'Survival') %>% 
    mutate(Measure = dplyr::case_when(
      Measure == 'shat'   ~ 'Predicted',
      Measure == 'surv_p' ~ 'Observed')
     ) %>% 
    ggplot() +
    geom_point(aes(x = age, 
                   y = Survival,
                   color = Measure),
               alpha = 0.75) +
    theme_bw() +
    facet_wrap(vars(!!sym(re)), ncol = 8) +
    lims(y = c(0,1)) +
    labs(x = 'Age',
         y = 'Survival', 
         title    = paste0(
           'Predicted vs observed survival probablity between' , re, 's'),
         subtitle = species) +
    theme(legend.position = 'top') +
    scale_color_colorblind()
  return(p_surv)
}

p_surv_quad <- su_fun(shat_quad_df, 'quad')
p_surv      <- su_fun(shat_df,      'year')

# store the plot
ggsave(paste0(result_dir, '/1_survival_years.png'),  
       plot = p_surv,      width = 6, height = 9, dpi = 150)
ggsave(paste0(result_dir, '/1_survival_quads.png'),  
       plot = p_surv_quad, width = 6, height = 9, dpi = 150)


# Recruitment two-stage model --------------------------------------------------
#   1. Model recruitment being there or not
#   2. Model per-capita recruitment (pcr) *conditional* on recruitment happening

# 1. yes or no recruits?
yesno_df <- recr_df %>% 
  mutate(recruits_yesno = as.numeric(recruits > 0)) %>% 
  dplyr::select(quad, year, parents, recruits_yesno) %>% 
  drop_na

# Yes-or-no recruits? (NOTE: parent number matters)
yesno_mod    <- glmer(recruits_yesno ~ (1 | year),
                       data = yesno_df,
                       family = binomial)

# predict zeroes
yesno_pred_df <- expand.grid(year    = coef(yesno_mod)$year %>% rownames,
                              stringsAsFactors = F) %>% 
  mutate(yesno_prob = predict(yesno_mod, 
                                newdata = .,
                                type = 'response') 
 ) 


# 2. per capita recruitment **conditional** on recruitment happening
pcr_df       <- recr_df %>% subset(!(recruits == 0))

# Per-capita recruitment model
pcr_mod      <- glmer(pcr ~ (1 | year),
                       data = pcr_df,
                       family = Gamma(link = 'log'))

# predict per capita recruitment
pcr_pred_df  <- expand.grid(year    = coef(pcr_mod)$year %>% rownames,
                             stringsAsFactors = F) %>% 
  # conditional (on recruitment == 1) per capita recruitment
  mutate(cond_pcr_hat = predict(pcr_mod, 
                                  newdata = .,
                                  type = 'response') 
 )

# 2.1 per capita recruitment unconditional 
pcr_uc_df       <- recr_df %>% 
  mutate(recruits_uc = replace(recruits,
                                recruits == 0,
                                0.1)) %>% 
  mutate(pcr_uc = recruits_uc / parents)


# Per-capita recruitment model
pcr_uc_mod      <- glmer(pcr_uc ~ (1 | year),
                          data = pcr_uc_df,
                          family = Gamma(link = 'log'))

# predict per capita recruitment
pcr_uc_pred_df  <- expand.grid(year    = coef(pcr_uc_mod)$year %>% rownames,
                                stringsAsFactors = F) %>% 
  # conditional (on recruitment == 1) per capita recruitment
  mutate(uc_pcr_hat = predict(pcr_uc_mod, 
                                newdata = .,
                                type = 'response') 
 )

# hurdle model predictions
recr_pred_df <- left_join(yesno_pred_df, pcr_uc_pred_df) %>% 
  left_join(pcr_pred_df) %>%
  mutate(pcr_uc_hat = uc_pcr_hat, 
          pcr_hat    = cond_pcr_hat * yesno_prob)

pcr_comp_mod <-
  recr_pred_df %>% 
  ggplot() +
  geom_point(aes(x = pcr_uc_hat, y = pcr_hat)) +
  geom_abline(intercept = 0, slope = 1, color = 'black', pch = 1) +
  labs(x = 'Unconditinal pcr hat',
       y = 'Conditional pcr hat', 
       title    = paste0(
         'Conditional on uncoditional per capita recruitment rate'),
       subtitle = species)

ggsave(paste0(result_dir, '/2_pcr_comp_mod.png'),  
       plot = pcr_comp_mod,
       width = 6, height = 9, dpi = 150)


# plot on recruitment
recr_p <- recr_df %>% 
  left_join(recr_pred_df) %>% 
  mutate(recr_hat = pcr_hat * parents) %>% 
  subset(year > 38 & year < 72) %>% 
  ggplot() + 
  geom_point(aes(x = parents,
                   y = recruits)) +
  geom_point(aes(x = parents,
                   y = recr_hat),
              color = 'red') +
  geom_abline(intercept = 0,
               slope     = 1) +
  facet_wrap(~ year, ncol = 4) +
  theme_bw() +
  labs(x = 'Recruits count',
       y = 'Parents count', 
       title    = paste0(
         'Observed vs modeled (red) recruitments'),
       subtitle = species)

# store the plot
ggsave(paste0(result_dir, '/2_recruitment_years.png'),  
       plot = recr_p,
       width = 6, height = 9, dpi = 150)


# compare 'naive' per capita recruitment, with the model's
pcr_comp_naive <- 
  recr_df %>% 
  group_by(year) %>% 
  summarise(pcr_naive = mean(pcr,na.rm=T)) %>% 
  ungroup %>% 
  left_join(recr_pred_df) %>% 
  ggplot() +
  geom_point(aes( x = pcr_naive,
                    y = pcr_hat)) +
  geom_point(aes( x = pcr_naive,
                    y = pcr_uc_hat), color = 'blue') +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'Native pcr',
       y = 'Modeled pcr', 
       title    = paste0(
         'Conditional vs uncoditional (blue)- on native per capita recruitment'),
       subtitle = species)

ggsave(paste0(result_dir, '/2_pcr_comp_naive.png'),  
       plot = pcr_comp_naive,
       width = 6, height = 9, dpi = 150)


# Matrix population model ------------------------------------------------------

# MPM parameters
mpm_df <- shat_df %>% 
  pivot_wider(names_from  = 'age',
               values_from = 'shat') %>% 
  left_join(dplyr::select(recr_pred_df,
                            year, pcr_hat)) %>% 
  drop_na %>% 
  ungroup

# lambda function
lambda_mpm <- function(x) x %>% eigen %>% .$value %>% Re %>% max

# year specific lambdas. 
#   Meant to work using rows of a data frame
year_lam <- function(year_x, mat_df){
  
  tmp_pars <- mat_df %>% subset(year == year_x) 
  
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
make_mat <- function(year_x, mat_df){
  
  tmp_pars <- mat_df %>% subset(year == year_x) 
  
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
  mutate(year = as.character(year)) %>% 
  mutate(age = dplyr::case_when(
    age == 0 ~ 0,
    age > 0  ~ 1) 
 ) %>% 
  group_by(year, age) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  # only retain year_t0 for which we have data
  right_join(dplyr::select(mpm_df,year)) %>% 
  drop_na %>% 
  pivot_wider(names_from  = 'age',
               values_from = 'n_t0',
               values_fill = 0) %>% 
  drop_na %>% 
  mutate(year = as.character(year))

# projected lambda (takes into account observed stage distribution)
proj_lam <- function(year_x, mat_df){
  
  # year-specific stage counts
  tmp_stage <- stage_counts %>% 
    subset(year == year_x) %>% 
    .[,c('0','1')]
  
  # year-specific projection matrix
  tmp_mat   <- make_mat(year_x, mat_df)
  
  nt0 <- sum(tmp_stage)
  nt1 <- (tmp_mat %*% as.numeric(tmp_stage)) %>% sum
  
  nt1 / nt0
  
}

# model lambdas
lam_df        <- mpm_df %>% 
  mutate(lam      = sapply(year, year_lam, mpm_df),
          proj_lam = sapply(year, proj_lam, mpm_df))

# Population counts at time t0
pop_counts_t0 <- df %>%
  group_by(year, quad) %>%
  summarize(n_t0 = n()) %>% 
  ungroup %>% 
  mutate(year = year + 1)

# Population counts at time t1
pop_counts_t1 <- df %>%
  group_by(year, quad) %>%
  summarize(n_t1 = n()) %>% 
  ungroup 

# Calculate observed population growth rates, 
#   accounting for discontinued sampling!
pop_counts <- left_join(pop_counts_t0, 
                        pop_counts_t1) %>% 
  mutate(year = year - 1) %>% 
  # by dropping NAs, we remove gaps in sampling!
  drop_na %>% 
  group_by(year) %>% 
  summarise(n_t0 = sum(n_t0),
             n_t1 = sum(n_t1)) %>% 
  ungroup %>% 
  mutate(obs_pgr = n_t1 / n_t0) %>% 
  mutate(year = as.character(year)) %>% 
  left_join(mpm_df) %>% 
  mutate(year = as.character(year)) %>% 
  left_join(lam_df) %>% 
  drop_na

write.csv(pop_counts, row.names = F, paste0(data_dir, '/metrics_df.csv'))

# mostly over-predicted population growth
#   not as bad as the published data though (see below)!
comp_obpogr_lam <- 
  ggplot(pop_counts) +
  geom_point(aes(obs_pgr, lam)) +
  geom_point(aes(obs_pgr, proj_lam),
              col = 'red') +
  geom_smooth(aes(obs_pgr, proj_lam), method = 'lm') +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'Observed population growth rates',
       y = 'Lambda', 
       title    = paste0(
         'Projected vs regular lambda on observed pop. growth rates'),
       subtitle = species)

ggsave(paste0(result_dir, '/3_comp_obpogr_lam.png'),  
       plot = comp_obpogr_lam,
       width = 6, height = 9, dpi = 150)

# for now, store everything as an R object that resembles 
mpm_df %>% 
  mutate(year = paste0('19',year)) %>% 
  mutate(SpeciesAuthor = species) %>% 
  setNames(c('year', 
              'surv_age0', 'surv_age1', 
              'per_capita_recruitment',
              'SpeciesAuthor')) %>% 
  write.csv(paste0(data_dir, '/mpm_pste.csv'), row.names = F)


# compare with COMPADRE's estimates --------------------------------------------

# load latest version of COMPADRE 
#   download here: https://compadre-db.org/Data/Compadre
load('adler_2007_ks/data/COMPADRE_v.6.23.5.0.RData')

# identify the matrices across COMPADRE 
id <- compadre$metadata %>% 
  subset(grepl(gsub(' ', '_', species), SpeciesAuthor)) %>% 
  subset(MatrixComposite == 'Individual') %>% 
  subset(MatrixStartYear %in% paste0('19', pop_counts$year)) %>% 
  row.names %>% 
  as.numeric

# years we can directly compare
yrs_mat   <- compadre$metadata$MatrixStartYear[id] %>% 
  gsub('19','',.)

# compadre's matrices
comp_mats <- compadre$mat[id] %>% setNames(yrs_mat)

# published matrices
pub_mat_df   <- data.frame(year = yrs_mat) %>% 
  # populate data frame with the three matrix elements
  mutate(`0`     = sapply(yrs_mat, function(x) pluck(comp_mats, x)$matA[2,1]),
          `1`     = sapply(yrs_mat, function(x) pluck(comp_mats, x)$matA[2,2]),
          pcr_hat = sapply(yrs_mat, function(x) pluck(comp_mats, x)$matA[1,1])
 ) %>% 
  # calculate lambda and projected lambda
  mutate(lam      = sapply(year, year_lam, .),
          proj_lam = sapply(year, proj_lam, .)) %>% 
  # re-set names to show these are from the publication
  setNames(c('year',    'pub_0',   'pub_1', 
              'pub_pcr_hat', 'pub_lam', 'pub_proj_lam')) %>% 
  drop_na

# compare population growth rates across estimations
compare_df <- pop_counts %>% 
  left_join(pub_mat_df) %>% 
  drop_na

# Recruitment: our estimation vs. publication
comp_pub_pcr <-
  ggplot(compare_df) +
  geom_point(aes(pcr_hat, pub_pcr_hat)) +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'pcr hat',
       y = 'pcr hat compadre', 
       title    = paste0(
         'PCR hat - estimation vs. publication'),
       subtitle = species)

ggsave(paste0(result_dir, '/4_comp_pub_pcr.png'),  
       plot = comp_pub_pcr,
       width = 6, height = 9, dpi = 150)

# Survival age 0: our estimation vs. publication
comp_pub_0 <-
  ggplot(compare_df) +
  geom_point(aes(`0`, pub_0)) +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'Suv 0',
       y = 'Suv 0 compadre', 
       title    = paste0(
         'Survival age 0 - estimation vs. publication'),
       subtitle = species)

ggsave(paste0(result_dir, '/4_comp_pub_0.png'),  
       plot = comp_pub_0,
       width = 6, height = 9, dpi = 150)

# Survival age 1: our estimation vs. publication
comp_pub_1 <-
  ggplot(compare_df) +
  geom_point(aes(`1`, pub_1)) +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'Suv 1',
       y = 'Suv 1 compadre', 
       title    = paste0(
         'Survival age 1 - estimation vs. publication'),
       subtitle = species)

ggsave(paste0(result_dir, '/4_comp_pub_1.png'),  
       plot = comp_pub_1,
       width = 6, height = 9, dpi = 150)

# Compare population growth rates
comp_pub_pgr <-
  ggplot(compare_df) +
  geom_point(aes(lam, pub_lam)) +
  geom_smooth(aes(lam, pub_lam), method = 'lm') +
  geom_abline(intercept = 0,
               slope     = 1) +
  labs(x = 'PGR',
       y = 'PGR compadre', 
       title    = paste0(
         'Population growth rates - estimation vs. publication'),
       subtitle = species)

ggsave(paste0(result_dir, '/4_comp_pub_pgr.png'),  
       plot = comp_pub_pgr,
       width = 6, height = 9, dpi = 150)

# compare population growth rates
comp_pgr <- 
  ggplot(compare_df) +
  # geom_point(aes(pcr_hat,pub_pcr_hat)) +
  geom_point(aes(obs_pgr, pub_lam), col = 'green') +
  geom_point(aes(obs_pgr, pub_proj_lam), col = 'blue') +
  geom_point(aes(obs_pgr, lam), col = 'darkred') +
  geom_point(aes(obs_pgr, proj_lam), col = 'red') +
  geom_abline(intercept = 0,
               slope     = 1) +
  geom_smooth(aes(obs_pgr, lam), method = 'lm') +
  labs(x = 'Lambda',
       y = 'Lambda compadre', 
       title    = paste0(
         'Lambda - estimation vs. publication'),
       subtitle = species)

ggsave(paste0(result_dir, '/4_comp_pgr.png'),  
       plot = comp_pgr,
       width = 6, height = 9, dpi = 150)

# calculate room mean squared error
rmse <- function(target, yhat) (target - yhat)^2 %>% mean %>% sqrt  

# data frame to ease computation
compute_df <- compare_df %>% 
  dplyr::select(obs_pgr, pub_lam, pub_proj_lam, 
                 lam, proj_lam) %>% 
  as.data.frame %>% 
  drop_na
