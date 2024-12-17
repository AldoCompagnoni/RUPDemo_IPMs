# MPM year specific MEAN Adler 2007 Kansas Ambrosia psilostachya

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com &
#  aldo.compagnoni@gmail.com
# Date: 2024.10.28

# Publication: https://doi.org/10.1890/0012-9658(2007)88[2673:LMQFKP]2.0.CO;2

# Original models from https://doi.org/10.1111/j.1365-2745.2009.01585.x
#   check the appendices for details on model fitting!

# Packages ---------------------------------------------------------------------
source( 'helper_functions/load_packages.R' )
load_packages( tidyverse, patchwork, skimr, lme4, ggthemes )


# Data -------------------------------------------------------------------------

# Define the species variable
species <- 'Ambrosia psilostachya'
sp_abb  <- tolower(
  gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                      collapse = "")))

# Species data frame
df <- read.csv(paste0('adler_2007_ks/data/', 
                      sp_abb, '/ks_', sp_abb, '.csv')) %>% 
  filter(Species == species) %>%
  select( -c(Suspect, nearEdge, Site) ) %>% # geometry, 
  mutate(across(c(Quad), as.factor)) %>% 
  rename(species  = Species, 
         survives = survives_tplus1,
         year     = Year,
         quad     = Quad,
         track_id = trackID)

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


# Survival ---------------------------------------------------------------------

# survival data. For this model, only two age classes
surv_df <- surv_out_df %>% 
  # consider all individuals past year 0 as if they had age 1
  mutate( age = dplyr::case_when(
    age == 0 ~ 0,
    age > 0  ~ 1 ) 
  )

surv_mean_df <- surv_df %>%
  subset( !is.na(age) ) %>%
  group_by( year, age) %>%
  summarise( surv_p = sum(survives) / n() ) %>%
  ungroup


ggplot() +
  geom_point(data = surv_mean_df, 
             aes(y= surv_p, x = age, colour = year))
ggplot() +
  geom_jitter(data = surv_df, 
             aes(y= survives, x = age, colour = year))


# Survival model: Random effect of year only
mod_surv   <- glm( surv_p ~ age,
                     data = mutate(surv_mean_df, age = as.factor(age)),
                     family = 'quasibinomial' )






# calculate year-specific mean predictions
shat0       <- data.frame( age  = 0,
                           stringsAsFactors = F ) %>%
  mutate( shat = coef(mod_surv)[1] %>% boot::inv.logit() )
shat1       <- data.frame( age  = 1,
                           stringsAsFactors = F ) %>%
  mutate( shat = (coef(mod_surv)[1] + coef(mod_surv)[2]) %>% 
            boot::inv.logit() )

shat_df     <- bind_rows( shat0, shat1 )


p_surv <- surv_df %>%
  subset( !is.na(age) ) %>%
  group_by( year, age ) %>%
  summarise( surv_p = sum(survives) / n() ) %>%
  ungroup %>%
  mutate( year = as.character(year) ) 


pred <- data_frame(year = rep(as.numeric(unique(p_surv$year)), 2),
                   age  = as.factor(c(
                     rep(0, length(unique(p_surv$year))), 
                     rep(1, length(unique(p_surv$year))))))
pred$pred_surv = boot::inv.logit(predict(mod_surv, pred))

  

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