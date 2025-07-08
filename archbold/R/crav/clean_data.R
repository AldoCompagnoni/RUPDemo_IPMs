df_og <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data_v2.csv')) %>% 
  janitor::clean_names() %>%  
  mutate(
    plant_id = as.factor(paste(site, quad, plant, sep = '_')),
    quad_id  = as.factor(paste(site, quad, sep = '_')),
    year     = as.numeric(substr(date, 1, 4)),  
    month    = as.numeric(substr(date, 6, 7)),
    site     = as.factor(site),
    quad     = as.factor(quad),
    mp       = as.factor(mp),
    plant    = as.factor(plant),
    caged    = as.factor(caged),
    veg      = as.factor(veg)) %>%
  arrange(site, quad, quad_id, plant, plant_id, year, month)


df_meta <- data.frame(variable = colnames(df)) %>% 
  mutate(definition = c(
    'study site', 'quadrat number',	'macroplot number',	
    'plant number (within quad)', 'direction within circular quad',	
    'distance from quad center',	'was quad caged from 2012 onward',	
    'vegetation type', 'year quadrat was initiated', 
    'year-month of observation',	'fire severity Dec. 2014 or Jan. 2015',
    'fire severity Aug. 2005',	'fire severity May/June 2009',
    "fire severity B's ridge 2016", 'fire severity Feb. 2017',
    'fire severity Oct. 2017', 'survival code for month', 'number of stems',
    'number of branch tips', 'number of flowers (corolla showing)', 
    'number of developing fruits', 'number of mature fruits', 'herbivory code',
    'plant identification', 'quadrat identification', 'sample year', 
    'sample month'))



# Recruits ---------------------------------------------------------------------
# s = 3 means Recruit/ we know that it showed up for the first time 
#  since the first recruits are recorded for 2010 we dont have to worrie about dormancy
df_og %>% 
  group_by(site, quad_id, plant_id) %>% 
  filter(any(s == 3)) %>% 
  ungroup()

df_og %>% 
  filter(s == 3) %>% 
  group_by(year) %>% 
  count


# Survival ---------------------------------------------------------------------
# The data starts with survival
names(df_og)

df_og %>% 
  group_by(site, quad_id, plant_id) %>% 
  filter(any(is.na(s))) %>% 
  ungroup()

df_og_exp <- df %>%
  # Survival:
  
  #  is the plant not appearing but it is 3 years before latest record, then NA
  mutate(date = ymd(paste0(date, "-01"))) %>%
  group_by(plant_id) %>%
  mutate(
    # s = 6 and 5 mean death (I think)
    latest_alive_date      = max(date[(s > 0 & s < 4)], na.rm = TRUE),
    earliest_recorded_date = min(date[ s > 0]         , na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    survival = case_when(
      #  is it before the earliest recorded data, then NA
      date <  earliest_recorded_date ~ NA_real_,
      #  is the plant still alive some time after, then 1
      date <  latest_alive_date & !is.infinite(latest_alive_date) ~ 1,
      #  is the plant not alive some time after, then 0
      date == latest_alive_date ~ 0,
      # is it after the latest date that the plant is alive, then NA
      date >  latest_alive_date ~ NA_real_,
      latest_alive_date == 2017 ~ NA_real_
    )) %>%
  # Recruits
  mutate(
    recruit = case_when(
      s == 3 ~ 1,
      date == earliest_recorded_date & year(earliest_recorded_date) > 1999 + 4 ~ 1,
      TRUE ~ NA_real_
    ))



# Mean data frame --------------------------------------------------------------
df_og_mean <- df_og_exp %>%
  group_by(site, quad_id , plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(survival )), NA_real_, min(survival,  na.rm = T)),
    size_t0  = if_else(all(is.na(br)),        NA_real_, max(br,        na.rm = T)),
    flower   = if_else(all(is.na(fl)),        NA_real_, max(fl,        na.rm = T)),  
    fruit    = if_else(all(is.na(fr)),        NA_real_, max(fr,        na.rm = T)),
    recruit  = if_else(all(is.na(recruit)),   NA_real_, min(recruit,   na.rm = T)),
    fire_sev = if_else(
      all(is.na(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f))),
      NA_real_, 
      mean(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f), na.rm = TRUE)),
    .groups  = 'drop'
  ) %>% 
  ungroup()


df_og_mean %>% 
  mutate(survives = if_else(survives == 0 & year > 2017 - 4, NA, survives))%>%
    # Define dormancy
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0)   ~ 1,
      size_t0  >  0 & !is.na(survives) ~ 0, 
      TRUE ~ NA_real_ 
    ),
    # Generate a new column 'dormancy_count' that counts consecutive 1s
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 & 
        lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1 & lag(dormancy, 5) == 1 ~ 6,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1 & lag(dormancy, 4) == 1                         ~ 5,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 
      & lag(dormancy, 3) == 1                                                 ~ 4,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1           ~ 3,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 0           ~ 2,
      dormancy == 1                                                           ~ 1,
      TRUE ~ dormancy
    )
  )


# Dormancy ---------------------------------------------------------------------
# In the year 2000 plant individual 1_12_3 goes dormant for 6 years,
#  it emerges with 10 stems (6 in the year before dormancy) in 2015
#  which in fact is a year after a fire (severity 1 for that plot in particular)
# Also in 2001 plant individual 1_12_9 goes dormant for 5 years, 
#  just to emerge with 9 stems (3 stems in 2000) in 2015 after a fire sev 1
