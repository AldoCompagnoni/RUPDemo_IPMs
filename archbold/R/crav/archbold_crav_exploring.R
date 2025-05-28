# Exploring data - Archbold - Menges 2016 - Crotalaria avonensis

# Author: Niklas Neisse*
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.03.03


# Website    : https://portal.edirepository.org/nis/mapbrowse?packageid=edi.219.1
# Publication: https://bioone.org/journals/southeastern-naturalist/volume-15/issue-3/058.015.0318/Ecology-and-Conservation-of-the-Endangered-Legume-Crotalaria-avonensis-in/10.1656/058.015.0318.short

# rm(list = ls())


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)

# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor, lme4)


# Specification ----------------------------------------------------------------
# Define head-directory 
v_head <- c('archbold')
# Define species
v_species <- c('Crotalaria avonensis')
# Customized delimiter for `read_delim` function, comma is predefined
custom_delimiter <- c()


# Create a unique species abbreviation for file naming
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Define script prefix
v_script_prefix <- str_c(v_head)

# Plot subtitle
v_ggp_suffix    <- paste(
  tools::toTitleCase(v_head), '-', v_species)


# Models
v_mod_set_gr <- c()
v_mod_set_su <- c()
v_mod_set_fl <- c()

# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_head))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

if (!dir.exists(dir_R     )) {dir.create(dir_R     )}
if (!dir.exists(dir_data  )) {dir.create(dir_data  )}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Functions --------------------------------------------------------------------
# function to plot your survival data 'binned' (instead of 'jittered')
source('helper_functions/plot_binned_prop.R')


# Data -------------------------------------------------------------------------
df <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data_v2.csv')) %>% 
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
  

# df <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data.csv')) %>% 
#   janitor::clean_names() %>%
#   mutate(across(c(5, 6, 11:length(.)), ~ na_if(., 9999))) %>%
#   mutate(
#     plant_id = as.factor(paste(site, quad, plant, sep = '_')),
#     quad_id  = as.factor(paste(site, quad, sep = '_')),
#     year     = as.numeric(substr(date, 1, 4)),  
#     month    = as.numeric(substr(date, 6, 7)),
#     site     = as.factor(site),
#     quad     = as.factor(quad),
#     mp       = as.factor(mp),
#     plant    = as.factor(plant),
#     caged    = as.factor(caged),
#     veg      = as.factor(veg)) %>%
#   arrange(site, quad, plant, year, month)

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

skimr::skim(df)
  

# Data exploration -------------------------------------------------------------
df %>% 
  group_by(site, date) %>% 
  summarise(nr_ind = length(.[2])) %>% 
  ungroup() %>% 
  # pivot_wider(names_from = year, values_from = nr_ind) %>%
  # Create a scatter plot of quadrat counts over the years
  ggplot() + 
  geom_point(aes(x = date, y = site)) +
  theme_bw() +
  labs(title    = 'Sampling inventory') +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.subtitle = element_text(size = 8))


# Prepare data frames for analysis ---------------------------------------------
odd_su_56 <- df %>% 
  filter(s > 4) %>%
  arrange(site, quad, plant, year, month)

odd_su_56_uniq <- odd_su_56 %>%
  distinct(site, quad, plant)

df_filtered <- df %>%
  semi_join(odd_su_56_uniq, by = c('site', 'quad', 'plant'))


#### What do the survival indices stand for?

## 3 = birth

## 6 
# 1,1,10 : 2005-06 -> death, no growth the same year, fire the year after 
# 4,18,2 : 2007-06 -> lives, no growth this month, growth the same year and the years after, fire two years before (s0)
# 4,19,2 : 2007-06 -> lives, no growth this month, growth the same year and the years after, fire two years before (s0)
# 4,19,4 : 2007-06 -> lives, no growth this month, growth the same year and the years after, fire two years before (s0)
## -> no data for this month?
### Exclude entirely from the dataset


## 5 
# 1,3,9  : 2006-02 -> birth, fire the year before (s2), dead the year after
# 1,7,10 : 2010-03 -> birth, fire the year before (s0), dead the same year
# 1,13,1 : 2006-06 -> birth, fire the year before (s1), lives many years after
# 3,48,8 : 2006-06 -> birth, fire the year before (s0), lives many years after
## -> birth after a fire in the year before
### Set the survival to 1 if they survive and 0 if they dont
###  and make a new column with age 1


# Original year mean dataframe -------------------------------------------------
df_mean_og <- df %>%
  filter(s < 6 | is.na(s)) %>%
  group_by(site, quad, quad_id, plant, plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(s )), NA_real_, max (s,  na.rm = TRUE)),
    size_t0  = if_else(all(is.na(br)), NA_real_, max (br, na.rm = TRUE)),
    fruit    = if_else(all(is.na(fr)), NA_real_, max (fr, na.rm = TRUE)),   
    flower   = if_else(all(is.na(fl)), NA_real_, max (fl, na.rm = TRUE)),  
    fire_sev = if_else(
      all(is.na(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f))),
      NA_real_, 
      mean(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f), na.rm = TRUE)),
    .groups  = 'drop'
  ) %>% 
  ungroup()


# Base mean dataframe ----------------------------------------------------------
df_mean <- df_mean_og %>% 
  group_by(site, quad, quad_id, plant, plant_id) %>% 
  # Handle survival based on previous dormancy status
  # If the current status is dead or NA
  # Check if survived in the previous 1, 2, or 3 years
  # Check if survives in the next 1, 2, or 3 years
  mutate(
    survives = if_else(
      (survives == 0 | is.na(survives)) & 
        (lag(survives, 1) %in% c(1, 3, 5) | 
           lag(survives, 2) %in% c(1, 3, 5) | 
           lag(survives, 3) %in% c(1, 3, 5)) & 
        (lead(survives, 1) %in% c(1, 3, 5) | 
           lead(survives, 2) %in% c(1, 3, 5) | 
           lead(survives, 3) %in% c(1, 3, 5)),
      1,  
      survives
    ),
    
    # Set survival to 0 if plant dies after 3 consecutive years of dormancy
    survives = if_else(
      survives == 1 & 
        lead(survives, 1) == 0 & 
        lead(survives, 2) == 0 & 
        lead(survives, 3) == 0, 
      0,  
      survives
    )
  ) %>%
  
  # Define dormancy
  mutate(
    dormancy = case_when(
      survives == 1 & is.na(size_t0) ~ 1,
      size_t0  >  0                  ~ 0, 
      TRUE ~ NA_real_ 
    ),
    # Generate a new column 'dormancy_count' that counts consecutive 1s
    dormancy_count = case_when(
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 0 ~ 2,
      dormancy == 1 & lag(dormancy, 1) == 1 & lag(dormancy, 2) == 1 ~ 3,
      dormancy == 1 ~ 1,
      TRUE ~ dormancy)
  ) %>% 
  
  # Define recruits
  mutate(
    recruit = case_when(
      (survives == 3 | survives == 5) ~ 1, 
      TRUE ~ NA_real_  
    )
  ) %>% 
  
  mutate(
    # Make sure that no survival value exceeds 1
    survives = if_else(survives > 1, 1, survives), 
    
    # Handle missing size_t0 for dead plants by setting NA in survival column
    survives = if_else(survives == 0 & is.na(size_t0), NA, survives)
  ) %>% 
  
  ungroup() %>% 
  arrange(site, quad, plant, year) %>%
  group_by(site, quad, plant, plant_id) %>%
  
  # Set size_t1 based on survival; propagate size_t0 if survives, otherwise set to NA
  mutate(
    size_t1 = case_when(
      survives == 1 ~ lead(size_t0), 
      TRUE ~ NA_real_  
    )
  ) %>%
  
  ungroup() %>% 
  
  # Compute log-transformed sizes and their powers for modeling
  mutate(
    logsize_t0   = log(size_t0),     
    logsize_t1   = log(size_t1),    
    logsize_t0_2 = logsize_t0^2,     
    logsize_t0_3 = logsize_t0^3      
  ) %>%
  group_by(site, quad, plant, plant_id) %>%
  mutate(age = NA_real_) %>%  # Initialize age as NA for all
  
  # Apply for loop to calculate age
  mutate(age = {
    # Initialize age vector
    age_vector <- rep(NA_real_, n())
    
    # Iterate over each row in the group (each site, quad, and plant combination)
    for (i in 1:n()) {
      if (!is.na(recruit[i]) && recruit[i] == 1) {
        age_vector[i] <- 1  # Set age to 1 for recruits
      } else if (!is.na(dormancy[i]) && dormancy[i] >= 0) {
        # If plant survives, increment the age based on previous value
        if (i > 1 && !is.na(age_vector[i - 1])) {
          age_vector[i] <- age_vector[i - 1] + 1
        }
      }
    }
    
    # Return the computed age vector
    age_vector
  }) %>%
  ungroup() %>%
  arrange(site, quad, plant, plant_id, year) %>%
  group_by(site, quad, plant, plant_id) %>%
  mutate(
    fire_event = ifelse(!is.na(fire_sev) & fire_sev > 0, 1, 0),
    fire_gap = {
      gap <- numeric(n())
      counter <- NA_real_
      for (i in seq_along(fire_event)) {
        if (is.na(fire_event[i])) {
          # Keep NA if we haven't seen a fire yet
          gap[i] <- counter
        } else if (fire_event[i] == 1) {
          counter <- 0
          gap[i] <- counter
        } else {
          counter <- ifelse(is.na(counter), NA_real_, counter + 1)
          gap[i] <- counter
        }
      }
      gap
    }
  ) %>%
  ungroup()


# BUGs to fix ------------------------------------------------------------------
df_issue <- df_mean[c(932:950,
                      1397:1406,
                      1418:1425,
                      1749:1759,
                      1768:1778,
                      2461:2470,
                      2902:2907,
                      3003:3019,
                      3432:3438,
                      3906:3912),]


# Histogram of recruit size (at age 1)
ggplot(df_mean %>% filter(recruit == 1), aes(x = size_t0)) +
  geom_histogram(binwidth = 1, fill = 'lightgray', color = 'black') +
  theme_minimal()

# Size over age of individuals with tracked recruit 
ggplot(df_mean, aes(x = age, y = logsize_t0)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +  
  theme_minimal() + 
  labs(x = 'Age', y = 'Log Size', title = 'Log Size vs Age of Recruits') + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, by = 1))


# Histogram of size after leaving dormancy
ggplot(df_mean %>% filter(dormancy == 1), aes(x = size_t1)) +
  geom_histogram(binwidth = 1, fill = 'lightgray', color = 'black') +
  theme_minimal()


# data frame for individuals that are flagged with potential recruits or dormancy
df_mean_f <- df_mean %>%
  group_by(site, quad, plant, plant_id) %>%
  mutate(
    flag_1 = if_else(
      any(year == 1999 & is.na(size_t0)) |
        (any(year %in% c(1999, 2000)) & all(is.na(size_t0[year %in% c(1999, 2000)]))) |
        (any(year %in% c(1999, 2000, 2001)) & all(is.na(size_t0[year %in% c(1999, 2000, 2001)]))),
      1, 0
    ),
    flag_2 = if_else(any(year %in% c(1999, 2000, 2001) & !is.na(size_t1)), 1, 0)
  ) %>%
  ungroup() %>% 
  filter(flag_1 == 1 & flag_2 == 1) 

# Initialize the age column
df_mean_f$age <- NA

# Loop through each unique combination of site, quad, and plant
for (i in 1:nrow(df_mean_f)) {
  
  current_site <- df_mean_f$site[i]
  current_quad <- df_mean_f$quad[i]
  current_plant <- df_mean_f$plant[i]
  
  # Subset the data for the current combination of site, quad, and plant
  group_data <- df_mean_f %>%
    filter(site == current_site, quad == current_quad, plant == current_plant) %>%
    arrange(year)  # Ensure the rows are ordered by year
  
  # Initialize the current age starting at NA
  current_age <- 1  # We will start at 1 when the first non-NA size_t0 is found
  
  # Loop through the rows within this group
  for (j in 1:nrow(group_data)) {
    if (!is.na(group_data$size_t0[j])) {
      # If the size_t0 is not NA, assign the current age
      df_mean_f$age[df_mean_f$site == current_site & 
                      df_mean_f$quad == current_quad & 
                      df_mean_f$plant == current_plant & 
                      df_mean_f$year == group_data$year[j]] <- current_age
      # Increment the age for the next year
      current_age <- current_age + 1
    } else if (j > 1) {
      # If size_t0 is NA, carry the age from the previous row
      previous_year_age <- df_mean_f$age[df_mean_f$site == current_site &
                                           df_mean_f$quad == current_quad &
                                           df_mean_f$plant == current_plant &
                                           df_mean_f$year == group_data$year[j - 1]]
      df_mean_f$age[df_mean_f$site == current_site & 
                      df_mean_f$quad == current_quad & 
                      df_mean_f$plant == current_plant & 
                      df_mean_f$year == group_data$year[j]] <- previous_year_age
    }
  }
}


# Histogram of recruit size (at 'first encounter')
ggplot(df_mean_f %>% filter(age == 1), aes(x = size_t0)) +
  geom_histogram(binwidth = 1, fill = 'lightgray', color = 'black') +
  theme_minimal()

# 
ggplot() +
  geom_jitter(
    data = df_mean_f, 
    aes(x = age, y = logsize_t0)) +
  theme_minimal()



# Fit the first model: logsize_t0 ~ age + (1 | plant_id) (linear model)
mod_flag_linear <- lmer(age ~ logsize_t0 + (1 | plant_id), data = df_mean_f)
summary(mod_flag_linear)

# Predict values using the linear model
df_mean_f$predicted_linear <- predict(mod_flag_linear, newdata = df_mean_f, re.form = NULL)

# Fit the second model: logsize_t0 ~ age + I(age^2) + (1 | plant_id) (quadratic model)
mod_flag_quad <- lmer(age ~ logsize_t0 + I(logsize_t0^2) + (1 | plant_id), data = df_mean_f)
summary(mod_flag_quad)

# Predict values using the quadratic model
df_mean_f$predicted_quad <- predict(mod_flag_quad, newdata = df_mean_f, re.form = NULL)

# Fit the third model: logsize_t0 ~ age + (1 | plant_id) (reversed model)
mod_flag_logsize_age <- lmer(logsize_t0 ~ age + (1 | plant_id), data = df_mean_f)
summary(mod_flag_logsize_age)

# Predict values using the logsize_t0 ~ age model
df_mean_f$predicted_logsize_age <- predict(mod_flag_logsize_age, newdata = df_mean_f, re.form = NULL)

# Fit the fourth model: logsize_t0 ~ age + I(age^2) + (1 | plant_id) (quadratic term for age)
mod_flag_quad_age <- lmer(logsize_t0 ~ age + I(age^2) + (1 | plant_id), data = df_mean_f)
summary(mod_flag_quad_age)

# Predict values using the logsize_t0 ~ age + I(age^2) model
df_mean_f$predicted_quad_age <- predict(mod_flag_quad_age, newdata = df_mean_f, re.form = NULL)

# Plot the first model (linear relationship)
ggplot(df_mean_f, aes(x = logsize_t0, y = age)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_linear, color = plant_id), alpha = 0.5) +  
  theme_minimal() +
  labs(title = 'Linear Model with Random Intercepts: logsize_t0 vs Age',
       y = 'Age',
       x = 'Log Size at t0') +
  theme(legend.position = 'none')

# Plot the second model (quadratic relationship for logsize_t0 ~ logsize_t0^2)
ggplot(df_mean_f, aes(x = logsize_t0, y = age)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_quad, color = plant_id), alpha = 0.5) +  
  theme_minimal() +
  labs(title = 'Quadratic Model with Random Intercepts: logsize_t0 vs Age',
       y = 'Age',
       x = 'Log Size at t0') +
  theme(legend.position = 'none')

# Plot the third model (logsize_t0 ~ age relationship)
ggplot(df_mean_f, aes(x = age, y = logsize_t0)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_logsize_age, color = plant_id), alpha = 0.5) +  
  theme_minimal() +
  labs(title = 'Log Size Model with Random Intercepts: Age vs Log Size at t0',
       y = 'Log Size at t0',
       x = 'Age') +
  theme(legend.position = 'none')

# Plot the fourth model (logsize_t0 ~ age + I(age^2) relationship)
ggplot(df_mean_f, aes(x = age, y = logsize_t0)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_quad_age, color = plant_id), alpha = 0.5) + 
  theme_minimal() +
  labs(title = 'Quadratic Model with Random Intercepts: Age vs Log Size at t0',
       y = 'Log Size at t0',
       x = 'Age') +
  theme(legend.position = 'none')



# Individuals to check out because of certains flags ---------------------------
#   1, 1, 23, 2017
# size t0 after recruit
#   1, 1, 23, 2012, 1, 14
#   1, 1, 37, 2014, 1, 13
#   1, 4, 17, 2012, 1, NA
#   1, 2, 22, 2010, 1, 19, 1, 8



# Investigation: Fire ----------------------------------------------------------
names(df_mean)
df_fire1 <- df_mean %>% 
  group_by(quad_id, year) %>% 
  summarise(fire_sev = mean(fire_sev, na.rm = T))

ggplot(df_fire1, aes(x = quad_id, y = year, size = fire_sev, color = fire_sev)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = 'Severity') +
  scale_color_continuous(name = '') +
  labs(x = 'Quadrat ID', y = 'Year', title = 'Fire Severity by Year and Plot') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Fire on Growth timing --------------------------------------------------------
names(df)
df_fi_hit <- df %>% 
  select(
    year, month, burn_a, burn_b, burn_c, burn_d, burn_e, burn_f, st, fl) %>% 
  mutate(burn_a = as.numeric(burn_a),
         burn_d = as.numeric(burn_d),
         burn_e = as.numeric(burn_e)) %>% 
  ungroup()

summary(df_fi_hit)

{
  # library(dplyr)
  # library(tidyr)
  # library(ggplot2)
  # library(patchwork)
  
  # Step 1: Extract months with fire (non-zero burn)
  fire_months <- df_fi_hit %>%
    select(month, starts_with('burn_')) %>%
    pivot_longer(cols = starts_with('burn_'), names_to = 'fire_type', values_to = 'value') %>%
    filter(!is.na(value), value > 0) %>%
    distinct(month) %>%
    pull(month)
  
  # Step 2: Filter NA data
  df_st <- df_fi_hit %>% filter(!is.na(st), !is.na(month))
  df_fl <- df_fi_hit %>% filter(!is.na(fl), !is.na(month))
  
  # Step 3: Plot ST (with fire lines)
  plot_st <- ggplot(df_st, aes(x = month, y = st)) +
    geom_jitter(height = 0.1, width = 0.1, alpha = 0.5, color = 'darkgreen') +
    geom_vline(xintercept = fire_months, color = 'red', linetype = 'dashed') +
    labs(y = 'stem count') +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # Step 4: Plot FL (inverted y + fire lines)
  plot_fl <- ggplot(df_fl, aes(x = month, y = fl)) +
    geom_jitter(height = 0.1, width = 0.1, alpha = 0.5, color = 'darkblue') +
    geom_vline(xintercept = fire_months, color = 'red', linetype = 'dashed') +
    scale_y_reverse() +
    labs(x = 'month', y = 'flower count') +
    theme_bw()
  
  # Step 5: Combine plots
  combined <- (plot_st / plot_fl) +
    plot_annotation(title = 'Timing of growth and flowering against fire events')
  
  # Save final plot to global environment
  assign('fig_fi_timing_overall', combined, envir = .GlobalEnv)
  }
fig_fi_timing_overall



fig_fi_timing <- {
  # library(dplyr)
  # library(tidyr)
  # library(ggplot2)
  # library(purrr)
  # library(patchwork)
  
  # Step 1: Identify all fire months and years (where any burn_* > 0)
  fires_by_year <- df_fi_hit %>%
    select(year, month, starts_with('burn_')) %>%
    pivot_longer(cols = starts_with('burn_'), names_to = 'fire_type', values_to = 'value') %>%
    filter(!is.na(value), value > 0) %>%
    distinct(year, month)
  
  fire_years <- unique(fires_by_year$year)
  
  # Step 2: Filter the dataset to only rows with non-NA st (growth counts)
  # but keep all rows for fire years, including those with NA st (for proper plot x axis range)
  df_st_fire_years <- df_fi_hit %>%
    filter(year %in% fire_years)
  
  df_st_nonfire_years <- df_fi_hit %>%
    filter(!year %in% fire_years, !is.na(st))
  
  # Step 3: Create plots per fire year (including months with NA st as gaps)
  plots_fire <- df_st_fire_years %>%
    split(.$year) %>%
    map(function(df_year) {
      yr <- unique(df_year$year)
      fire_months <- fires_by_year %>%
        filter(year == yr) %>%
        pull(month)
      
      ggplot(df_year, aes(x = month, y = st)) +
        geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = 'darkgreen', na.rm = TRUE) +
        geom_smooth(method = 'glm', method.args = list(family = 'poisson'), color = 'black', na.rm = TRUE) +
        geom_vline(xintercept = fire_months, color = 'red', linetype = 'dashed') +
        scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
        labs(title = paste('Year:', yr), y = 'Stem count', x = 'Month') +
        theme_bw()
    })
  
  # Step 4: Plot for non-fire years
  plot_nonfire <- ggplot(df_st_nonfire_years, aes(x = month, y = st)) +
    geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = 'gray40') +
    geom_smooth(method = 'glm', method.args = list(family = 'poisson'), color = 'black') +
    scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
    labs(title = 'All Non-Fire Years Combined', y = 'Stem count', x = 'Month') +
    theme_bw()
  
  # Step 5: Combine all plots in a 2-column layout
  wrap_plots(c(plots_fire, list(plot_nonfire)), ncol = 2) +
    plot_annotation(title = 'Stem Count per Month: Fire vs Non-Fire Years with GLM Fit')
}
fig_fi_timing



fig_fi_fl_timing <- {
  # Step 1: Identify fire months and years (same as before)
  fires_by_year <- df_fi_hit %>%
    select(year, month, starts_with('burn_')) %>%
    pivot_longer(cols = starts_with('burn_'), names_to = 'fire_type', values_to = 'value') %>%
    filter(!is.na(value), value > 0) %>%
    distinct(year, month)
  
  fire_years <- unique(fires_by_year$year)
  
  # Step 2: Filter dataset for flower count (fl)
  # Keep all rows for fire years (including NAs for full x-axis coverage)
  df_fl_fire_years <- df_fi_hit %>%
    filter(year %in% fire_years)
  
  # Non-fire years only with non-NA fl
  df_fl_nonfire_years <- df_fi_hit %>%
    filter(!year %in% fire_years, !is.na(fl))
  
  # Step 3: Create plots for each fire year
  plots_fire_fl <- df_fl_fire_years %>%
    split(.$year) %>%
    map(function(df_year) {
      yr <- unique(df_year$year)
      fire_months <- fires_by_year %>%
        filter(year == yr) %>%
        pull(month)
      
      ggplot(df_year, aes(x = month, y = fl)) +
        geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = 'darkblue', na.rm = TRUE) +
        geom_smooth(method = 'glm', method.args = list(family = 'poisson'), color = 'black', na.rm = TRUE) +
        geom_vline(xintercept = fire_months, color = 'red', linetype = 'dashed') +
        scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
        labs(title = paste('Year:', yr), y = 'Flower count', x = 'Month') +
        theme_bw()
    })
  
  # Step 4: Plot for non-fire years
  plot_nonfire_fl <- ggplot(df_fl_nonfire_years, aes(x = month, y = fl)) +
    geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = 'gray40') +
    geom_smooth(method = 'glm', method.args = list(family = 'poisson'), color = 'black') +
    scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
    labs(title = 'All Non-Fire Years Combined', y = 'Flower count', x = 'Month') +
    theme_bw()
  
  # Step 5: Combine plots in 2-column layout
  wrap_plots(c(plots_fire_fl, list(plot_nonfire_fl)), ncol = 2) +
    plot_annotation(title = 'Flower Count per Month: Fire vs Non-Fire Years with GLM Fit')
}
fig_fi_fl_timing


# Fire on Growth ---------------------------------------------------------------
df_fire_grow <- df_mean %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(plant_id, year, size_t0, size_t1, age,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, 
         fire_event, fire_sev) %>% 
  mutate(fire_event = as.factor(fire_event),
         fire_sev   = if_else(is.na(fire_sev), 0, fire_sev))

df_fire_grow %>% count(fire_event == 1)

ggplot(data = df_fire_grow, aes(x = logsize_t0, y = logsize_t1)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.3, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.8, size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', aes(colour = fire_event)) +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red')) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = 'Fire on Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1]))



# Model fire on growth with event ----------------------------------------------
mod_fi_gr_0.0 <- lm(logsize_t1 ~ 1, 
                  data = df_fire_grow)

# Linear model
mod_fi_gr_1.0 <- lm(logsize_t1 ~ logsize_t0, 
                  data = df_fire_grow)

# Linear model
mod_fi_gr_1.1 <- lm(logsize_t1 ~ logsize_t0 + fire_event, 
                  data = df_fire_grow)

# Linear model
mod_fi_gr_1 <- lm(logsize_t1 ~ logsize_t0 +
                    fire_event + fire_event:logsize_t0, 
                  data = df_fire_grow)

# Quadratic model
mod_fi_gr_2.0 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + fire_event, 
                  data = df_fire_grow)  

# Quadratic model
mod_fi_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                    fire_event + fire_event:logsize_t0 + 
                    fire_event:logsize_t0_2, 
                  data = df_fire_grow)  

# Cubic model
mod_fi_gr_3.0 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + fire_event, 
                  data = df_fire_grow)


# Intercept model
mod_fi_gr_0 <- lm(logsize_t1 ~ fire_event, 
                  data = df_fire_grow)
# Linear model
mod_fi_gr_1 <- lm(logsize_t1 ~ logsize_t0 +
                    fire_event + fire_event:logsize_t0, 
                  data = df_fire_grow)
# Quadratic model
mod_fi_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                    fire_event + fire_event:logsize_t0 + 
                    fire_event:logsize_t0_2, 
                  data = df_fire_grow)  
# Cubic model
mod_fi_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + 
                    fire_event + fire_event:logsize_t0 + 
                    fire_event:logsize_t0_2 + fire_event:logsize_t0_3, 
                  data = df_fire_grow)

mods_fi_gr      <- list(mod_fi_gr_0,   mod_fi_gr_1,   mod_fi_gr_2,   mod_fi_gr_3,
                        mod_fi_gr_0.0, mod_fi_gr_1.0, mod_fi_gr_2.0, mod_fi_gr_3.0)
mods_fi_gr_dAIC <- AICtab(mods_fi_gr, weights = T, sort = F)$dAIC


# Get the sorted indices of dAIC values
mods_fi_gr_sorted <- order(mods_fi_gr_dAIC)

# Establish the index of model complexity
v_mod_fi_set_gr <- c()
if (length(v_mod_fi_set_gr) == 0) {
  mod_fi_gr_index_bestfit <- mods_fi_gr_sorted[1]
  v_mod_fi_gr_index       <- mod_fi_gr_index_bestfit - 1 
} else {
  mod_fi_gr_index_bestfit <- v_mod_fi_set_gr +1
  v_mod_fi_gr_index       <- v_mod_fi_set_gr
}

mod_fi_gr_bestfit         <- mods_fi_gr[[mod_fi_gr_index_bestfit]]
mod_fi_gr_ranef           <- coef(mod_fi_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_fire_grow$pred <- predict(mod_fi_gr_bestfit, type = 'response')

source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

g_fi_grow_line <- ggplot(
  df_fire_grow, aes(x = logsize_t0, y = logsize_t1, colour = fire_event)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +
  geom_function(fun = function(x) predictor_fun(x, mod_fi_gr_ranef),
                color = line_color_pred_fun(mod_fi_gr_ranef),
                lwd = 2) +
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red')) +
  labs(title    = 'Growth prediction with fire event',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_fi_grow_pred <- ggplot(
  df_fire_grow, aes(x = pred, y = logsize_t1, colour = fire_event)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1),
              color = 'red', lwd = 2) +
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red'))

g_fi_grow_overall_pred <- g_fi_grow_line + g_fi_grow_pred + plot_layout()
g_fi_grow_overall_pred


# Model fire on growth with event & severity -----------------------------------
# Intercept model
mod_fi2_gr_0 <- lm(logsize_t1 ~ fire_sev + fire_event + fire_sev:fire_event, 
                   data = df_fire_grow)
# Linear model
mod_fi2_gr_1 <- lm(logsize_t1 ~ logsize_t0 +
                     fire_sev:fire_event + 
                     fire_sev   + fire_sev:logsize_t0 +
                     fire_event + fire_event:logsize_t0, 
                   data = df_fire_grow)
# Quadratic model
mod_fi2_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 +
                     fire_sev:fire_event + 
                     fire_sev   + fire_sev:logsize_t0 + 
                     fire_sev   + fire_sev:logsize_t0_2 +
                     fire_event + fire_event:logsize_t0 +
                     fire_event + fire_event:logsize_t0_2, 
                   data = df_fire_grow)  
# Cubic model
mod_fi2_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +
                     fire_sev:fire_event + 
                     fire_sev   + fire_sev:logsize_t0 + 
                     fire_sev   + fire_sev:logsize_t0_2 + 
                     fire_sev   + fire_sev:logsize_t0_3 +
                     fire_event + fire_event:logsize_t0 +
                     fire_event + fire_event:logsize_t0_2 +
                     fire_event + fire_event:logsize_t0_3, 
                   data = df_fire_grow)

mods_fi2_gr      <- list(mod_fi2_gr_0, mod_fi2_gr_1, mod_fi2_gr_2, mod_fi2_gr_3)
mods_fi2_gr_dAIC <- AICtab(mods_fi2_gr, weights = T, sort = F)$dAIC


# Get the sorted indices of dAIC values
mods_fi2_gr_sorted <- order(mods_fi2_gr_dAIC)

# Establish the index of model complexity
v_mod_fi2_set_gr <- c()
if (length(v_mod_fi2_set_gr) == 0) {
  mod_fi2_gr_index_bestfit <- mods_fi2_gr_sorted[1]
  v_mod_fi2_gr_index       <- mod_fi2_gr_index_bestfit - 1 
} else {
  mod_fi2_gr_index_bestfit <- v_mod_fi2_set_gr +1
  v_mod_fi2_gr_index       <- v_mod_fi2_set_gr
}

mod_fi2_gr_bestfit         <- mods_fi2_gr[[mod_fi2_gr_index_bestfit]]
mod_fi2_gr_ranef           <- coef(mod_fi2_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_fire_grow$pred2 <- predict(mod_fi2_gr_bestfit, type = 'response')

source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

g_fi2_grow_line <- ggplot(
  df_fire_grow, aes(x = logsize_t0, y = logsize_t1, colour = fire_event)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +
  geom_function(fun = function(x) predictor_fun(x, mod_fi2_gr_ranef), 
                color = line_color_pred_fun(mod_fi2_gr_ranef), 
                lwd = 2) +
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red')) + 
  labs(title    = 'Growth prediction with fire event',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_fi2_grow_pred <- ggplot(
  df_fire_grow, aes(x = pred2, y = logsize_t1)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red'))

g_fi2_grow_overall_pred <- g_fi2_grow_line + g_fi2_grow_pred + plot_layout() 
g_fi2_grow_overall_pred


# Model fire on growth with severity -------------------------------------------
mod_fi3_gr_0 <- lm(logsize_t1 ~ fire_sev, 
                   data = df_fire_grow)
# Linear model
mod_fi3_gr_1 <- lm(logsize_t1 ~ logsize_t0 +
                     fire_sev + fire_sev:logsize_t0, 
                   data = df_fire_grow)
# Quadratic model
mod_fi3_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                     fire_sev + fire_sev:logsize_t0 + 
                     fire_sev:logsize_t0_2, 
                   data = df_fire_grow)  
# Cubic model
mod_fi3_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + 
                     fire_sev + fire_sev:logsize_t0 + 
                     fire_sev:logsize_t0_2 + fire_sev:logsize_t0_3, 
                   data = df_fire_grow)

mods_fi3_gr      <- list(mod_fi3_gr_0, mod_fi3_gr_1, mod_fi3_gr_2, mod_fi3_gr_3)
mods_fi3_gr_dAIC <- AICtab(mods_fi3_gr, weights = T, sort = F)$dAIC


# Get the sorted indices of dAIC values
mods_fi3_gr_sorted <- order(mods_fi3_gr_dAIC)

# Establish the index of model complexity
v_mod_fi3_set_gr <- c()
if (length(v_mod_fi3_set_gr) == 0) {
  mod_fi3_gr_index_bestfit <- mods_fi3_gr_sorted[1]
  v_mod_fi3_gr_index       <- mod_fi3_gr_index_bestfit - 1 
} else {
  mod_fi3_gr_index_bestfit <- v_mod_fi3_set_gr +1
  v_mod_fi3_gr_index       <- v_mod_fi3_set_gr
}

mod_fi3_gr_bestfit         <- mods_fi3_gr[[mod_fi3_gr_index_bestfit]]
mod_fi3_gr_ranef           <- coef(mod_fi3_gr_bestfit)

summary(mod_fi3_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_fire_grow$pred3 <- predict(mod_fi3_gr_bestfit, type = 'response')

source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

g_fi3_grow_line <- ggplot(
  df_fire_grow, aes(x = logsize_t0, y = logsize_t1, colour = fire_event)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +
  geom_function(fun = function(x) predictor_fun(x, mod_fi3_gr_ranef), 
                color = line_color_pred_fun(mod_fi3_gr_ranef), 
                lwd = 2) +
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red')) + 
  labs(title    = 'Growth prediction with fire event',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_fi3_grow_pred <- ggplot(
  df_fire_grow, aes(x = pred3, y = logsize_t1)) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 0),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.2, size = 1.5) +
  geom_jitter(data = subset(df_fire_grow, fire_event == 1),
              aes(colour = fire_event, shape = fire_event),
              alpha = 0.6, size = 1.5) +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw() +
  scale_colour_manual(values = c('0' = 'black', '1' = 'red'))

g_fi3_grow_overall_pred <- g_fi3_grow_line + g_fi3_grow_pred + plot_layout() 
g_fi3_grow_overall_pred


# Fire on growth models comparison ---------------------------------------------

mods_fi_gr_all <- list(
  mod_fi_gr_0,  mod_fi_gr_1,  mod_fi_gr_2,  mod_fi_gr_3,
  mod_fi2_gr_0, mod_fi2_gr_1, mod_fi2_gr_2, mod_fi2_gr_3,
  mod_fi3_gr_0, mod_fi3_gr_1, mod_fi3_gr_2, mod_fi3_gr_3
)

mods_fi_gr_all_dAIC <- AICtab(mods_fi_gr_all, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fi_gr_all_sorted <- order(mods_fi_gr_all_dAIC)

# Establish the index of model complexity
mods_fi_gr_all_index_bestfit <- mods_fi_gr_all_sorted[1]

mods_fi_gr_all_bestfit         <- mods_fi_gr_all[[mods_fi_gr_all_index_bestfit]]
mods_fi_gr_all_ranef           <- coef(mods_fi_gr_all_bestfit)

summary(mods_fi_gr_all_bestfit)



# Fire on Survival -------------------------------------------------------------
df_fire_surv_0 <- df_mean %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire_event) %>% 
  filter(fire_event == 0)

df_fire_surv_1 <- df_mean %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire_event) %>% 
  filter(fire_event == 1)

p_fire_suv_0 <- ggplot(
  data = plot_binned_prop(df_fire_surv_0, 10, logsize_t0, survives)) +
  geom_point(aes(x = logsize_t0,
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Fire on Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Survival to time t1')) +
  theme(plot.subtitle = element_text(size = 8))

p_fire_suv_1 <- ggplot(
  data = plot_binned_prop(df_fire_surv_1, 10, logsize_t0, survives)) +
  geom_point(aes(x = logsize_t0,
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Fire on Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Survival to time t1')) +
  theme(plot.subtitle = element_text(size = 8))

p_fire_suv <- p_fire_suv_0 + p_fire_suv_1 + plot_layout()
p_fire_suv


# Fire and Survival data
df_fi_su <- df_mean %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3,
         fire_event)

# Intercept model
mod_fi_su_0   <- glm(survives ~ 1,
                     data = df_fi_su, family = 'binomial')
# Logistic regression
mod_fi_su_1   <- glm(survives ~ logsize_t0,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_1.1 <- glm(survives ~ logsize_t0 + fire_event,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_1.2 <- glm(survives ~ logsize_t0 + fire_event +
                       logsize_t0 : fire_event,
                     data = df_fi_su, family = 'binomial')
# Quadratic logistic model
mod_fi_su_2   <- glm(survives ~ logsize_t0 + logsize_t0_2,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_2.1 <- glm(survives ~ logsize_t0 + logsize_t0_2 +  fire_event,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_2.2 <- glm(survives ~ logsize_t0 + logsize_t0_2 +  fire_event +
                       logsize_t0   : fire_event +
                       logsize_t0_2 : fire_event +
                       logsize_t0   : logsize_t0_2,
                     data = df_fi_su, family = 'binomial')
# Cubic logistic model
mod_fi_su_3   <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_3.1 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +  fire_event,
                     data = df_fi_su, family = 'binomial')
mod_fi_su_3.2 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 +  fire_event +
                       logsize_t0   : fire_event   +
                       logsize_t0_2 : fire_event   +
                       logsize_t0_3 : fire_event   +
                       logsize_t0   : logsize_t0_2 +
                       logsize_t0   : logsize_t0_3 +
                       logsize_t0_2 : logsize_t0_3,
                     data = df_fi_su, family = 'binomial')


# Compare models using AIC
mods_fi_su      <- list(mod_fi_su_0,
                        mod_fi_su_1, mod_fi_su_1.1, mod_fi_su_1.2,
                        mod_fi_su_2, mod_fi_su_2.1, mod_fi_su_2.2,
                        mod_fi_su_3, mod_fi_su_3.1, mod_fi_su_3.2)

mods_fi_su_dAIC <- AICtab(mods_fi_su, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fi_su_sorted <- order(mods_fi_su_dAIC)

v_mod_set_fi_su <- c()

# Establish the index of model complexity
if (length(v_mod_set_fi_su) == 0) {
  mod_fi_su_index_bestfit <- mods_fi_su_sorted[1]
  v_mod_fi_su_index       <- mod_fi_su_index_bestfit - 1 
} else {
  mod_fi_su_index_bestfit <- v_mod_set_fi_su +1
  v_mod_fi_su_index       <- v_mod_set_fi_su
}


mod_fi_su_bestfit   <- mods_fi_su[[mod_fi_su_index_bestfit]]
mod_fi_su_ranef     <- coef(mod_fi_su_bestfit)



# Fire on Flowers --------------------------------------------------------------
df_fire_fl <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(
    flower, fire_sev, fire_event, fire_gap, 
    logsize_t0, logsize_t0_2, logsize_t0_3) %>% 
  group_by(site, quad_id, year) %>% 
  summarise(fl_quad = mean(flower, na.rm = T), 
            fire_gap = max(fire_gap, na.rm = T))


# Plot with counts
ggplot(data = df_fire_fl, aes(x = as.factor(fire_gap), y = fl_quad)) + 
  geom_boxplot() +
  geom_text(data = df_fire_fl %>%
              group_by(fire_gap) %>%
              summarise(n = n()), 
            aes(
    x = as.factor(fire_gap), 
    y = max(df_fire_fl$fl_quad, na.rm = TRUE) + 1.5, label = n),
    inherit.aes = FALSE, size = 3) + 
  theme_minimal() +
  labs(y = expression('mean nr flowers per quadrat'),
       x = expression('year gap after fire'))


# Fire on Flower ---------------------------------------------------------------
df_fire_fl_0 <- df_mean %>% 
  filter(!is.na(flower)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, flower, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire_event) %>% 
  mutate(flower = if_else(flower > 0, 1, flower)) %>% 
  filter(fire_event == 0)

df_fire_fl_1 <- df_mean %>% 
  filter(!is.na(flower)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, flower, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire_event) %>% 
  mutate(flower = if_else(flower > 0, 1, flower)) %>% 
  filter(fire_event == 1)

p_fire_fl_0 <- ggplot(
  data = plot_binned_prop(df_fire_fl_0, 10, logsize_t0, flower)) +
  geom_point(aes(x = logsize_t0,
                 y = flower),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Fire on Flowering without fire',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Flowering at time t0')) +
  theme(plot.subtitle = element_text(size = 8))

p_fire_fl_1 <- ggplot(
  data = plot_binned_prop(df_fire_fl_1, 10, logsize_t0, flower)) +
  geom_point(aes(x = logsize_t0,
                 y = flower),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Fire on Flowering with fire',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Flowering at time t0')) +
  theme(plot.subtitle = element_text(size = 8))

p_fire_fl <- p_fire_fl_0 + p_fire_fl_1 + plot_layout()
p_fire_fl


# Fire on Flowering one year after ---------------------------------------------
df_fr_fl_t0 <- df_mean %>% 
  group_by(site, quad_id, plant_id) %>% 
  select(year, flower, fire_event) %>% 
  rename(fl_t0 = flower)

df_fr_fl_t1 <- df_mean %>% 
  group_by(site, quad_id, plant_id) %>% 
  select(year, flower) %>% 
  mutate(year = year - 1) %>% 
  rename(fl_t1 = flower)
  
df_fr_fl_t <- df_fr_fl_t0 %>% 
  left_join(df_fr_fl_t1, by = c('site', 'quad_id', 'plant_id', 'year')) %>% 
  mutate(fire_event = as.factor(fire_event))

ggplot(data = df_fr_fl_t, aes(x = fl_t0, y = fl_t1)) +
  geom_jitter(data = subset(df_fr_fl_t, fire_event == 0),
              aes(colour = factor(fire_event)),
              alpha = 0.2) +
  geom_jitter(data = subset(df_fr_fl_t, fire_event == 1),
              aes(colour = factor(fire_event)),
              alpha = 0.6) +
  geom_smooth(aes(colour = factor(fire_event)), method = 'lm') +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(
    values = c('0' = 'black', '1' = 'red'),
    name = 'Fire Event',
    labels = c('No Fire', 'Fire')
  ) +
  theme_minimal() +
  labs(subtitle = v_ggp_suffix,
       x        = expression('Flowers at time t0'),
       y        = expression('Flowers at time t1')) +
  theme(plot.subtitle = element_text(size = 8))



# Fire on Flowering on plant id level data -------------------------------------
df_fi_fl_p <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(flower, fire_event, logsize_t0, logsize_t0_2, logsize_t0_3) %>% 
  rename(flower_t0 = flower) %>% 
  filter(flower_t0 >= 0) %>% 
  mutate(flower_t0 = if_else(flower_t0 > 0, 1, flower_t0)) %>% 
  left_join(df_mean %>% 
              group_by(site, quad_id, plant_id, year) %>% 
              select(fire_event) %>% 
              mutate(year = year + 1) %>% 
              rename(fire_event_t_1 = fire_event),
            by = c('site', 'quad_id', 'plant_id', 'year'))


df_fi_fl_p %>% 
  summary()


# Fire on flowering ------------------------------------------------------------

ggplot(data = df_fi_fl_p, aes(y = flower_t0, x = as.factor(fire_event), colour = logsize_t0)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.6, size = 1.5) +
  scale_color_viridis_c(option = 'D', name = 'Log Size (t0)') +
  theme_minimal(base_size = 12) +
  labs(
    x = 'Fire Event',
    y = 'Flowering (binary)',
    title    = 'Flowering by Fire Event',
    subtitle = v_ggp_suffix) +
  theme(
    plot.title = element_text(face = 'bold'),
    axis.title = element_text(face = 'bold'),
    legend.position = 'right',
    plot.subtitle = element_text(size = 8)) 

  
# Fire on flowering models -----------------------------------------------------
# Control models without fire --------------------------------------------------
mod_fi_fl_control0   <- glm(flower_t0 ~ 1,
                            data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_control1   <- glm(flower_t0 ~ logsize_t0,
                            data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_control2.1 <- glm(flower_t0 ~ logsize_t0 + logsize_t0_2,
                            data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_control2.2 <- glm(flower_t0 ~ logsize_t0 * logsize_t0_2,
                            data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_control3.1 <- glm(flower_t0 ~ logsize_t0 + logsize_t0_2 + 
                              logsize_t0_3,
                            data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_control3.2 <- glm(flower_t0 ~ logsize_t0 * logsize_t0_2 * 
                              logsize_t0_3,
                            data = df_fi_fl_p, family = 'binomial')

# Fire the same year -----------------------------------------------------------
# Intercept model
mod_fi_fl_0.0 <- glm(flower_t0 ~ fire_event,
                   data = df_fi_fl_p, family = 'binomial') 
# Logistic regression
mod_fi_fl_1.0 <- glm(flower_t0 ~ fire_event + 
                     logsize_t0 +
                     fire_event:logsize_t0,
                   data = df_fi_fl_p, family = 'binomial') 
mod_fi_fl_1.01 <- glm(flower_t0 ~ fire_event + 
                       logsize_t0,
                     data = df_fi_fl_p, family = 'binomial') 
# Quadratic logistic model
mod_fi_fl_2.0 <- glm(flower_t0 ~ fire_event + 
                     logsize_t0 + logsize_t0_2 +
                     fire_event:logsize_t0 + 
                     fire_event:logsize_t0_2,
                   data = df_fi_fl_p, family = 'binomial')  
# Quadratic logistic model
mod_fi_fl_2.01 <- glm(flower_t0 ~ fire_event + 
                     logsize_t0 + logsize_t0_2 +
                     fire_event:logsize_t0_2,
                   data = df_fi_fl_p, family = 'binomial') 
# Quadratic logistic model
mod_fi_fl_2.02 <- glm(flower_t0 ~ fire_event + 
                       logsize_t0 + logsize_t0_2,
                     data = df_fi_fl_p, family = 'binomial')  
# Cubic logistic model
mod_fi_fl_3.0 <- glm(flower_t0 ~ fire_event + 
                     logsize_t0 + logsize_t0_2 + logsize_t0_3 +
                     fire_event:logsize_t0 + 
                     fire_event:logsize_t0_2 + 
                     fire_event:logsize_t0_3,
                   data = df_fi_fl_p, family = 'binomial')  


# Fire the year before on flowering models -------------------------------------
# Intercept model
mod_fi_fl_0.1 <- glm(flower_t0 ~ fire_event_t_1,
                     data = df_fi_fl_p, family = 'binomial') 
# Logistic regression
mod_fi_fl_1.1 <- glm(flower_t0 ~ fire_event_t_1 + 
                       logsize_t0 +
                       fire_event_t_1:logsize_t0,
                     data = df_fi_fl_p, family = 'binomial') 

mod_fi_fl_1.11 <- glm(flower_t0 ~ fire_event_t_1 + 
                       logsize_t0,
                     data = df_fi_fl_p, family = 'binomial') 
# Quadratic logistic model
mod_fi_fl_2.1 <- glm(flower_t0 ~ fire_event_t_1 + 
                       logsize_t0 + logsize_t0_2 +
                       fire_event_t_1:logsize_t0 + 
                       fire_event_t_1:logsize_t0_2,
                     data = df_fi_fl_p, family = 'binomial')  
# Quadratic logistic model
mod_fi_fl_2.11 <- glm(flower_t0 ~ fire_event_t_1 + 
                        logsize_t0 + logsize_t0_2 +
                        fire_event_t_1:logsize_t0_2,
                      data = df_fi_fl_p, family = 'binomial') 
# Quadratic logistic model
mod_fi_fl_2.12 <- glm(flower_t0 ~ fire_event_t_1 + 
                        logsize_t0 + logsize_t0_2,
                      data = df_fi_fl_p, family = 'binomial')  
# Cubic logistic model
mod_fi_fl_3.1 <- glm(flower_t0 ~ fire_event_t_1 + 
                       logsize_t0 + logsize_t0_2 + logsize_t0_3 +
                       fire_event_t_1:logsize_t0 + 
                       fire_event_t_1:logsize_t0_2 + 
                       fire_event_t_1:logsize_t0_3,
                     data = df_fi_fl_p, family = 'binomial')  


# Fire and fire the year before on flowering -----------------------------------
# Quadratic model
mod_fi_fl_4.2 <- glm(flower_t0 ~ fire_event + fire_event_t_1 +
                     logsize_t0 + logsize_t0_2,
                   data = df_fi_fl_p, family = 'binomial')

mod_fi_fl_4.21 <- glm(flower_t0 ~ fire_event * fire_event_t_1 +
                       logsize_t0 + logsize_t0_2 +
                       fire_event:logsize_t0   + fire_event_t_1:logsize_t0 +
                       fire_event:logsize_t0_2 + fire_event_t_1:logsize_t0_2,
                     data = df_fi_fl_p, family = 'binomial')

# Compare models using AIC
mods_fi_fl      <- list(
  mod_fi_fl_control0, mod_fi_fl_control1, 
  mod_fi_fl_control2.1, mod_fi_fl_control2.2,
  mod_fi_fl_control3.1, mod_fi_fl_control3.2,  
  mod_fi_fl_0.0, mod_fi_fl_1.0, mod_fi_fl_1.01, mod_fi_fl_2.0, mod_fi_fl_3.0,
  mod_fi_fl_2.01, mod_fi_fl_2.02,
  mod_fi_fl_0.1, mod_fi_fl_1.1, mod_fi_fl_1.11, mod_fi_fl_2.1, mod_fi_fl_3.1,
  mod_fi_fl_2.11, mod_fi_fl_2.12,
  mod_fi_fl_4.2, mod_fi_fl_4.21)

mods_fi_fl_dAIC <- bbmle::AICtab(mods_fi_fl, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fi_fl_sorted <- order(mods_fi_fl_dAIC)

v_mod_sets_fi_fl <- c()

# Establish the index of model complexity
if (length(v_mod_sets_fi_fl) == 0) {
  mod_fi_fl_index_bestfit <- mods_fi_fl_sorted[1]
  v_mod_fi_fl_index       <- mod_fi_fl_index_bestfit - 1 
} else {
  mod_fi_fl_index_bestfit <- v_mod_sets_fi_fl +1
  v_mod_fi_fl_index       <- v_mod_sets_fi_fl
}


mod_fi_fl_bestfit   <- mods_fi_fl[[mod_fi_fl_index_bestfit]]
mod_fi_fl_ranef     <- coef(mod_fi_fl_bestfit)

summary(mod_fi_fl_bestfit)

# anova(glm(flower_t0 ~ logsize_t0 + logsize_t0_2, 
#           family = 'binomial', 
#           data = df_fi_fl_p %>% drop_na()),
#       glm(flower_t0 ~ fire_event + 
#             logsize_t0 + logsize_t0_2, 
#           family = 'binomial',
#           data = df_fi_fl_p %>% drop_na()),
#       mod_fi_fl_bestfit, test = 'Chisq')
# 
# anova(glm(flower_t0 ~ fire_event + 
#             logsize_t0 + logsize_t0_2, 
#           family = 'binomial',
#           data = df_fi_fl_p %>% drop_na()),
#       mod_fi_fl_bestfit, test = 'Chisq')


# Generate predictions for survival across a range of sizes
mod_fi_fl_x <- seq(
  min(df_fi_fl_p$logsize_t0, na.rm = T),
  max(df_fi_fl_p$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_fi_flow_pred <- predictor_fun(mod_fi_fl_x, mod_fi_fl_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_fi_fl_x, flower_t0 = .)

g_fi_flow_line <- ggplot() +
  geom_jitter(data = df_fi_fl_p, 
              aes(x = logsize_t0, y = flower_t0, colour = as.factor(fire_event)),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_fi_flow_pred, 
            aes(x = logsize_t0, y = flower_t0),
            color = line_color_pred_fun(mod_fi_fl_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Flowering prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_fi_flow_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_fi_fl_p, 10, logsize_t0, flower_t0), 
    aes(x = logsize_t0, y = flower_t0) ) +
  geom_errorbar(
    data = plot_binned_prop(df_fi_fl_p, 10, logsize_t0, flower_t0), 
    aes(x = logsize_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_fi_flow_pred, 
            aes(x = logsize_t0, y = flower_t0),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
g_fi_flow_overall_pred <- g_fi_flow_line + g_fi_flow_bin + plot_layout()
g_fi_flow_overall_pred



df_fi_fl_ls_pred <- expand.grid(
  logsize_t0 = seq(
    min(df_fi_fl_p$logsize_t0, na.rm = TRUE),
    max(df_fi_fl_p$logsize_t0, na.rm = TRUE), length.out = 100),
  fire_event = c(0, 1),
  fire_event_t_1 = mean(df_fi_fl_p$fire_event_t_1, na.rm = TRUE)) %>% 
  mutate(logsize_t0_2 = rep(seq(
    min(df_fi_fl_p$logsize_t0_2, na.rm = TRUE),
    max(df_fi_fl_p$logsize_t0_2, na.rm = TRUE), length.out = 100),2))

# Predict with standard errors
pred <- predict(mod_fi_fl_bestfit, df_fi_fl_ls_pred, type = 'link', se.fit = TRUE)

# Add predicted values and CIs (back-transform from logit scale)
df_fi_fl_ls_pred$predicted_prob <- plogis(pred$fit)
df_fi_fl_ls_pred$lower_ci <- plogis(pred$fit - 1.96 * pred$se.fit)
df_fi_fl_ls_pred$upper_ci <- plogis(pred$fit + 1.96 * pred$se.fit)

# Plot
ggplot() +
  # Raw data
  geom_jitter(data = df_fi_fl_p, 
              aes(x = logsize_t0, y = flower_t0), 
              alpha = 0.2, width = 0, height = 0.05, color = 'black') +
  
  # Confidence ribbon
  geom_ribbon(data = df_fi_fl_ls_pred, 
              aes(x = logsize_t0, ymin = lower_ci, ymax = upper_ci),
              fill = 'blue', alpha = 0.2) +
  
  # Prediction line
  geom_line(data = df_fi_fl_ls_pred, 
            aes(x = logsize_t0, y = predicted_prob), 
            color = 'blue', size = 1.2) +
  
  facet_wrap(~ fire_event, labeller = labeller(fire_event = c('0' = 'No Fire', '1' = 'Fire'))) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = 'Predicted Flowering by Size and Fire Event',
    subtitle = v_ggp_suffix,
    x = 'Log(Size at t0)',
    y = 'Flowering') +
  theme(
    plot.title = element_text(face = 'bold', size = 16),
    strip.text = element_text(face = 'bold'),
    plot.subtitle = element_text(size = 8))



# ------------------------------------------------------------------------------
# Prepare prediction grid
logsize_seq <- seq(min(df_fi_fl_p$logsize_t0, na.rm = TRUE), 
                   max(df_fi_fl_p$logsize_t0, na.rm = TRUE), 
                   length.out = 100)

newdata <- expand.grid(
  fire_event = c(0, 1),
  fire_event_t_1 = mean(df_fi_fl_p$fire_event_t_1, na.rm = TRUE),
  logsize_t0 = logsize_seq
)
newdata$logsize_t0_2 <- newdata$logsize_t0^2

# Get predicted probabilities
newdata$predicted_prob <- predict(mod_fi_fl_bestfit, newdata, type = "response")

# Plot with raw data and custom colors and alpha
ggplot() +
  # Raw data points (jittered) with different alpha values based on fire_event
  geom_jitter(data = df_fi_fl_p, 
              aes(x = logsize_t0, y = flower_t0, color = factor(fire_event),
                  alpha = factor(fire_event)),  # map alpha to fire_event
              width = 0, height = 0.1, size = 1.5) +
  
  # Model prediction lines with custom colors
  geom_line(data = newdata, 
            aes(x = logsize_t0, y = predicted_prob, color = factor(fire_event)),
            size = 1.2) +
  
  # Custom color scale for fire_event
  scale_color_manual(values = c("0" = "black", "1" = "red")) +
  
  # Custom alpha scale: no fire (0) with alpha 0.3, fire (1) with alpha 0.9
  scale_alpha_manual(values = c("0" = 0.1, "1" = 1)) +
  
  labs(x = "log(size at t0)",
       y = "Flowering Probability",
       color = "Fire Event (t0)") +
  theme_minimal()
# ------------------------------------------------------------------------------


# Fire on Fruits ---------------------------------------------------------------
df_fire_fr <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(fruit, fire_sev, fire_event, fire_gap) %>% 
  group_by(site, quad_id, year) %>% 
  summarise(fr_quad = mean(fruit, na.rm = T), 
            fire_gap = as.factor(max(fire_gap, na.rm = T)))

# Plot
ggplot(data = df_fire_fr, aes(x = fire_gap, y = fr_quad)) + 
  geom_boxplot() +
  geom_text(data = df_fire_fr %>%
              group_by(fire_gap) %>%
              summarise(n = n(),
                        y_pos = max(fr_quad, na.rm = TRUE) + 0.2), 
            aes(x = fire_gap, y = y_pos, label = n),
            inherit.aes = FALSE, size = 3) +
  theme_minimal() +
  labs(y = expression('mean nr fruits per quadrat'),
       x = expression('year gap after fire')) +
  coord_cartesian(ylim = c(0, max(df_fire_fr$fr_quad, na.rm = TRUE) + 1))


# Effect of fire on fruiting
df_fi_quad <- df_mean %>%  
  group_by (year, site, quad) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

df_fi_group <- df_fi_quad %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

df_fi_cover <- left_join(df_fi_quad, df_fi_group) %>%
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()


df_fi_fruit <- df_mean %>%
  group_by (year, site, quad) %>% 
  summarise(nr_quad    = mean(fruit, na.rm = T),
            fire_event = max( fire_event, na.rm = T )) %>% 
  ungroup

df_fi_fruit <- left_join(df_fi_cover, df_fi_fruit) %>% 
  mutate(fire_event = as.factor(if_else(is.na(fire_event), 0, fire_event)))

summary(df_fi_fruit$fire_event)

ggplot(
  df_fi_fruit, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1) +  
  theme_bw() + 
  labs(title    = 'Fruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of fruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8)) +
  facet_wrap('fire_event')


# Fire on Fruiting one year after ----------------------------------------------
df_fi_fr_t0 <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(fruit, fire_event) %>% 
  rename(fr_t0 = fruit)

df_fi_fr_t1 <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(fruit) %>% 
  mutate(year = year - 1) %>% 
  rename(fr_t1 = fruit)

df_fi_fr_t <- df_fi_fr_t0 %>% 
  left_join(df_fi_fr_t1) %>% 
  mutate(fire_event = if_else(is.na(fire_event), 0, fire_event),
         fire_event = as.factor(fire_event))

ggplot(data = df_fi_fr_t, aes(x = fr_t0, y = fr_t1, color = fire_event)) +
  geom_jitter(alpha = 0.4, width = 0.1, height = 0.1, size = 2) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(
    values = c('0' = 'black', '1' = 'red'),
    name = 'Fire Event',
    labels = c('No Fire', 'Fire')) +
  theme_minimal(base_size = 14) +
  labs(
    title = 'Fruit Production: Current vs. Previous Year',
    x = 'Fruit in Year t',
    y = 'Fruit in Year t+1') +
  theme(
    plot.title = element_text(face = 'bold', size = 16),
    legend.title = element_text(face = 'bold'),
    legend.position = 'top')
  


# Fire on Flower to Fruit transition -------------------------------------------
df_fi_ftf_og <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(logsize_t0, flower, fruit, fire_event) %>% 
  mutate(ftf = fruit/flower) %>% 
  ungroup()

df_fi_ftf_og %>%
  filter(!is.na(fruit), !is.na(flower)) %>%               
  mutate(ex_fr = if_else(fruit > flower, 1, 0)) %>%       
  count(ex_fr) %>%                                        
  mutate(prop = n / sum(n))

df_fi_ftf_og %>%
  filter(!is.na(fruit), !is.na(flower)) %>% 
  arrange(flower)

df_fi_ftf_1 <- df_fi_ftf_og %>% 
  filter(!(ftf == 0 & flower == 0))

df_fi_ftf_1 %>%
  filter(!is.na(fruit), !is.na(flower)) %>%               
  mutate(ex_fr = if_else(fruit > flower, 1, 0)) %>%       
  count(ex_fr) %>%                                        
  mutate(prop = n / sum(n))

df_fi_ftf_1 %>%
  filter(!is.na(fruit), !is.na(flower)) %>%               
  mutate(ex_fr = if_else(fruit > flower, 1, 0)) %>%
  filter(ftf > 0) %>% 
  count(ex_fr) %>%                                        
  mutate(prop = n / sum(n))


df_fi_ftf <- df_fi_ftf_1 %>% 
  mutate(fruit = if_else(fruit > flower, flower, fruit),
         ftf   = fruit/flower)

df_fi_ftf %>% 
  arrange(desc(fruit))

df_fi_ftf_og %>% 
  filter(flower > 0) %>% 
  arrange(desc(ftf))

df_fi_ftf_og %>% 
  filter(flower == 0) %>% 
  arrange(desc(ftf)) %>% 
  head(5)

df_fi_ftf_og %>% 
  filter(flower > 0) %>% 
  arrange(desc(ftf)) %>% 
  head(5)

ggplot(data = df_fi_ftf_og) +
  geom_histogram(aes(x = ftf))

ggplot(data = df_fi_ftf) +
  geom_histogram(aes(x = ftf))

ggplot(data = df_fi_ftf_og) +
  geom_histogram(aes(x = ftf)) + 
  facet_wrap('fire_event')


ggplot(data = df_fi_ftf) +
  geom_histogram(aes(x = ftf)) + 
  facet_wrap('fire_event')


ggplot(data = df_fi_ftf_og %>% filter(flower > 0)) +
  geom_jitter(aes(y = ftf, x = logsize_t0)) +
  facet_wrap('fire_event')

ggplot(data = df_fi_ftf %>% filter(flower > 0)) +
  geom_jitter(aes(x = logsize_t0, y = ftf), 
              width = 0.1, height = 0.05, alpha = 0.4, color = 'black') +
  facet_wrap(~ fire_event, labeller = labeller(fire_event = c('0' = 'No Fire', '1' = 'Fire'))) +
  theme_minimal(base_size = 14) +
  labs(
    title = 'Fruit-to-Flower Ratio by Size and Fire Event',
    x = 'Log(Size at t0)',
    y = 'Fruit-to-Flower Ratio') +
  theme(
    plot.title = element_text(face = 'bold', size = 16),
    strip.text = element_text(face = 'bold'))


# Fire on FtF-R models ---------------------------------------------------------

# Binomial model with and without fire even
mod0 <- glm( cbind(fruit,flower-fruit) ~ 1, 
             data = subset( df_fi_ftf, !is.nan(ftf) ), 
             family = 'binomial') 
mod1 <- glm( cbind(fruit,flower-fruit) ~ fire_event, 
             data = subset( df_fi_ftf, !is.nan(ftf) ), 
             family = 'binomial') 
mod2 <- glm( cbind(fruit,flower-fruit) ~ logsize_t0, 
             data = subset( df_fi_ftf, !is.nan(ftf) ), 
             family = 'binomial') 
mod3 <- glm( cbind(fruit,flower-fruit) ~ logsize_t0 + fire_event, 
             data = subset( df_fi_ftf, !is.nan(ftf) ), 
             family = 'binomial') 
mod4 <- glm( cbind(fruit,flower-fruit) ~ logsize_t0*fire_event, 
             data = subset( df_fi_ftf, !is.nan(ftf) ), 
             family = 'binomial') 


# Compare the two models
bbmle::AICctab(mod0,mod1,mod2,mod3,mod4, weights = T)

# confidence intervals are enormous. This is a mess.
confint(mod3)

# Vast majority of data points are 0!!!
#   I don't think we have enough data to confidently model fire effect
subset( df_fi_ftf, !is.nan(ftf) ) %>% 
  count( fire_event, ftf)


# Fire on Recruits the same year -----------------------------------------------
df_quad <- df_mean %>%  
  group_by (year, site, quad) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

df_group <- df_quad %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

df_cover <- left_join(df_quad, df_group) %>%
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

df_fire_recr <- df_mean %>%
  group_by (year, site, quad) %>% 
  summarise(nr_quad = sum(recruit, na.rm = T), 
            fire_event = max(fire_event, na.rm = T)) %>% 
  ungroup

df_fire_recr <- left_join(df_cover, df_fire_recr) %>% 
  filter(!is.na(fire_event))

ggplot(
  df_fire_recr, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Fire on Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8)) + 
  facet_wrap(~ fire_event)




df_fi_re_yn <- df_fire_recr %>% 
  mutate(re_yn = if_else(nr_quad > 0, 1, nr_quad))


df_mod_fi_re_yn <- df_fi_re_yn %>% filter(!is.na(re_yn))
# Fit a negative binomial model for recruitment
mod_fi_re_0 <- glm(re_yn ~ 1, data = df_mod_fi_re_yn, family = 'binomial')
mod_fi_re_1 <- glm(re_yn ~ fire_event, data = df_mod_fi_re_yn, family = 'binomial')

bbmle::AICtab(mod_fi_re_0, mod_fi_re_1, weights = T)


# Compare models using AIC
mods_fi_re      <- list(mod_fi_re_0, mod_fi_re_1)

mods_fi_re_dAIC <- AICtab(mods_fi_re, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fi_re_sorted <- order(mods_fi_re_dAIC)

v_mod_set_fi_re <- c()

# Establish the index of model complexity
if (length(v_mod_set_fi_re) == 0) {
  mod_fi_re_index_bestfit <- mods_fi_re_sorted[1]
  v_mod_fi_re_index       <- mod_fi_re_index_bestfit - 1 
} else {
  mod_fi_re_index_bestfit <- v_mod_set_fi_re +1
  v_mod_fi_re_index       <- v_mod_set_fi_re
}


mod_fi_re_bestfit   <- mods_fi_re[[mod_fi_re_index_bestfit]]
mod_fi_re_ranef     <- coef(mod_fi_re_bestfit)


# Generate predictions for recruitment
df_fi_re_mod <- df_fi_re_yn %>% 
  mutate(mod_pred = predict(mod_fi_re_bestfit, type = 'response')) 

# Summarize total number of recruits and predictions
df_rec_sums_m <- df_fi_re_mod %>%
  summarize(nr_quad  = sum(nr_quad),
            mod_pred = sum(mod_pred))

# Count number of adult individuals
indiv_m <- df_surv %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
repr_pc_m <- indiv_m %>%
  bind_cols(df_rec_sums_m) %>%
  mutate(repr_pc_mean = mod_pred / n_adults) %>%
  mutate(repr_pc_obs = nr_quad / n_adults) %>%
  drop_na 



# Fire on Recruits the next year -----------------------------------------------
df_fire_re <- df_mean %>% 
  group_by(site, quad_id, plant_id, year) %>% 
  select(recruit, fire_sev, fire_event, fire_gap) %>% 
  group_by(site, quad_id, year) %>% 
  summarise(re_quad = sum(recruit, na.rm = T), 
            fire_gap = max(fire_gap, na.rm = T))

# Plot with counts
ggplot(data = df_fire_re, aes(x = as.factor(fire_gap), y = re_quad)) + 
  geom_boxplot() +
  geom_text(data = df_fire_re %>%
              group_by(fire_gap) %>%
              summarise(n = n()), 
            aes(
              x = as.factor(fire_gap), 
              y = max(df_fire_re$re_quad, na.rm = TRUE) + 1.5, label = n),
            inherit.aes = FALSE, size = 3) + 
  theme_minimal() +
  labs(y = expression('sum recuits per quadrat'),
       x = expression('year gap after fire'))

ggplot(data = df_fire_re %>%
         filter(fire_gap >= 0,
                fire_gap <= 1), 
       aes(x = as.factor(fire_gap), y = re_quad)) + 
  geom_boxplot() +
  geom_jitter(width = 0.05, height = 0.2) +
  theme_minimal()


summary(lm(data = df_fire_re %>% filter(fire_gap >= 0, fire_gap <= 1),
           formula = re_quad ~ fire_gap))



# Fire on Per-capita Recruit ---------------------------------------------------
df_fire_rec_pc <- df_mean %>% 
  group_by(site, quad_id, year) %>% 
  summarise(rec_nr_t1 = sum(recruit, na.rm = T), 
            fire_gap = max(fire_gap, na.rm = T)) %>% 
  left_join(df_mean %>%
              group_by(site, quad_id, plant_id, year) %>%
              summarise(rep_nr_t0 = as.integer(any(fruit > 0)), .groups = 'drop') %>% 
              ungroup() %>% 
              group_by(site, quad_id, year) %>% 
              summarise(rep_nr_t0 = sum(rep_nr_t0, na.rm = T)) %>% 
              mutate(year = year + 1)
            , by = c('site', 'quad_id', 'year')) %>% 
  mutate(rec_pc = rec_nr_t1 / rep_nr_t0) %>%
  mutate(rec_pc = ifelse(is.nan(rec_pc), 0, rec_pc))
  

ggplot(df_fire_rec_pc <- df_fire_rec_pc %>%
         mutate(fire_gap = as.factor(fire_gap)), aes(x = fire_gap, y = rec_pc)) +
  geom_boxplot() +
  geom_text(data = df_fire_rec_pc %>%
              group_by(fire_gap) %>%
              summarise(n = n(),
                        y_pos = min(max(rec_pc, na.rm = TRUE) + 0.2, 5)),
            aes(x = fire_gap, y = y_pos, label = n),
            inherit.aes = FALSE, size = 3) +
  theme_minimal() +
  labs(y = expression('per capita recruits'),
       x = expression('year gap after fire')) +
  coord_cartesian(ylim = c(0, 5))


# Data Processes ---------------------------------------------------------------
names(df_mean)


# Growth data ------------------------------------------------------------------
df_grow <- df_mean %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(plant_id, year, size_t0, size_t1, age,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)


ggplot(
  data  = df_grow, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))



# Survival data ----------------------------------------------------------------
df_surv <- df_mean %>% 
  filter(!is.na(survives)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

g_surv_overall <- ggplot(
  data = plot_binned_prop(df_mean, 10, logsize_t0, survives)) +
  geom_point(aes(x = logsize_t0,
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Survival to time t1')) +
  theme(plot.subtitle = element_text(size = 8))
g_surv_overall

ggplot(data = df_surv) +
  geom_jitter(aes(x = logsize_t0, y = survives), 
              position = position_jitter(width = 0.1, height = 0.3)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Survival',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('Survival to time t1')) +
  theme(plot.subtitle = element_text(size = 8))


# Recruitment data -------------------------------------------------------------
df_quad <- df_mean %>%  
  group_by (year, site, quad) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

df_group <- df_quad %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

df_cover <- left_join(df_quad, df_group) %>%
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

df_recr <- df_mean %>%
  group_by (year, site, quad) %>% 
  summarise(nr_quad = sum(recruit, na.rm = T)) %>% 
  ungroup

df_recr <- left_join(df_cover, df_recr)


ggplot(
  df_recr, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))


df_re_qd <- df_mean %>% 
  group_by(site, quad_id, year) %>%
  select(recruit) %>% 
  summarise(rec_qd_t1 = sum(recruit, na.rm = T)) %>%
  left_join(df_mean %>% 
              group_by(site, quad_id, year) %>% 
              summarise(nr_ind = sum(!is.na(size_t0))) %>% 
              mutate(year = year - 1),
            by = c('site', 'quad_id', 'year'))

ggplot(data = df_re_qd) + 
  geom_jitter(aes(y = rec_qd_t1, x = nr_ind)) + 
  geom_smooth(aes(y = rec_qd_t1, x = nr_ind), method = 'lm') + 
  theme_minimal() + 
  labs(title    = 'Recruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of recruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))


# Fruiting data ----------------------------------------------------------------
df_quad <- df_mean %>%  
  group_by (year, site, quad) %>% 
  summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
  ungroup

df_group <- df_quad %>% 
  group_by (year) %>% 
  summarise(g_cov = mean(tot_p_area)) %>% 
  ungroup

df_cover <- left_join(df_quad, df_group) %>%
  mutate(year = year + 1) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na()

hist(df$fl)

df_fruit <- df_mean %>%
  group_by (year, site, quad) %>% 
  summarise(nr_quad = sum(fruit, na.rm = T)) %>% 
  ungroup

df_fruit <- left_join(df_cover, df_fruit)

ggplot(
  df_fruit, aes(x = tot_p_area, y = nr_quad)) + 
  geom_point(alpha = 0.5, pch = 16, size = 1, color = 'red') +  
  theme_bw() + 
  labs(title    = 'Fruitment',
       subtitle = v_ggp_suffix,
       x        = expression('Total parent plant area '[t0]),   
       y        = expression('Number of fruits '     [t1])) +
  theme(plot.subtitle = element_text(size = 8))


# Fecundity data ---------------------------------------------------------------

df_fecu <- df_fruit %>%
  rename(nr_fruit = nr_quad) %>% 
  left_join(
    df_recr %>%
      mutate(year = year - 1) %>%
      rename(nr_rec_t1 = nr_quad) %>% 
      select('year', 'site', 'quad', 'nr_rec_t1'),
    by = c('year', 'site', 'quad')
  )


ggplot(data = df_fruit %>%
         rename(nr_fruit = nr_quad) %>% 
         left_join(
           df_recr %>%
             rename(nr_rec_t0 = nr_quad) %>% 
             select('year', 'site', 'quad', 'nr_rec_t0'),
           by = c('year', 'site', 'quad')
         )) +
  geom_point(aes(x = nr_fruit, y = nr_rec_t0))


ggplot(data = df_fecu) +
  geom_point(aes(x = nr_fruit, y = nr_rec_t1))



ggplot(data = df_fruit %>%
         rename(nr_fruit = nr_quad) %>% 
         left_join(
           df_recr %>%
             mutate(year = year - 2) %>% 
             select('year', 'site', 'quad', 'nr_quad'),
           by = c('year', 'site', 'quad')
         ) %>% 
       rename(nr_rec_t2 = nr_quad)
       ) +
  geom_point(aes(x = nr_fruit, y = nr_rec_t2))



df_fecu_mod <- df_recr %>%
  rename(nr_recr = nr_quad) %>% 
  left_join(
    df_fruit %>%
      mutate(year = year + 1) %>% 
      select('year', 'site', 'quad', 'nr_quad'),
    by = c('year', 'site', 'quad')
  ) %>% 
  rename(nr_fru_t1 = nr_quad) %>% 
  left_join(
    df_fruit %>%
      mutate(year = year + 2) %>% 
      select('year', 'site', 'quad', 'nr_quad'),
    by = c('year', 'site', 'quad')
  ) %>% 
  rename(nr_fru_t2 = nr_quad) %>% 
  left_join(
    df_fruit %>%
      mutate(year = year + 3) %>% 
      select('year', 'site', 'quad', 'nr_quad'),
    by = c('year', 'site', 'quad')
  ) %>% 
  rename(nr_fru_t3 = nr_quad) %>% 
  left_join(
    df_fruit %>%
      mutate(year = year + 4) %>% 
      select('year', 'site', 'quad', 'nr_quad'),
    by = c('year', 'site', 'quad')
  ) %>% 
  rename(nr_fru_t4 = nr_quad)


ggplot(df_fecu_mod %>%
         pivot_longer(
           cols = starts_with('nr_fru_t'),
           names_to = 'fru_type',
           values_to = 'fru_count'
         ), 
       aes(x = fru_count, y = nr_recr, color = fru_type)) +
  geom_smooth(alpha = 0.2, method = 'lm') +
  geom_jitter(aes(shape = fru_type)) +
  labs(
    x = 'Number of Fruits',
    y = 'Number of Recruits',
    color = 'Fruit Type',
    title = 'Relationship between Fruit Counts and Recruits'
  ) +
  theme_minimal()


names(df_fecu_mod)

mod_fec <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
              data = df_fecu_mod)

summary(mod_fec)



# Fecundity year specific ------------------------------------------------------
#  fires are highlighted

df_mean %>%
  group_by(year) %>%
  summarize(total_fire = sum(fire_sev, na.rm = TRUE)) %>%
  filter(total_fire > 0)


df_fecu_mod %>%
  mutate(
    quad_id = paste(site, quad, sep = '_'),
    year_fire = ifelse(year %in% c(2005, 2009, 2014, 2016, 2017), 'fire', 'normal')
  ) %>%
  ggplot(aes(x = nr_fru_t1, y = nr_recr, color = year_fire)) +
  geom_jitter(alpha = 0.7) +
  facet_wrap(~year) +
  scale_color_manual(values = c('fire' = 'red', 'normal' = 'gray')) +
  theme_minimal()

# Model
mod_fec <- lmer(nr_recr ~ nr_fru_t1 + (1|year), data = df_fecu_mod)
summary(mod_fec)



# Fec - site Level -------------------------------------------------------------
df_fecu_mod_site <- df_fecu_mod %>%
  group_by(year, site) %>% 
  summarise(
    tot_p_area = sum(tot_p_area, na.rm = TRUE),
    g_cov      = sum(g_cov, na.rm = TRUE),
    nr_recr    = sum(nr_recr, na.rm = TRUE),
    nr_fru_t1  = sum(nr_fru_t1, na.rm = TRUE),
    nr_fru_t2  = sum(nr_fru_t2, na.rm = TRUE),
    nr_fru_t3  = sum(nr_fru_t3, na.rm = TRUE),
    nr_fru_t4  = sum(nr_fru_t4, na.rm = TRUE))


ggplot(df_fecu_mod_site %>%
         pivot_longer(
           cols = starts_with('nr_fru_t'),
           names_to = 'fruit_type',
           values_to = 'fruit_count'
         ), aes(x = fruit_count, y = nr_recr, color = fruit_type)) +
  geom_smooth(se = FALSE, method = 'lm') +
  geom_jitter(aes(shape = fruit_type)) +
  labs(
    title = 'Relationship between Fruit Counts and Recruits',
    x = 'Fruit Count',
    y = 'Number of Recruits',
    color = 'Fruit Type'
  ) +
  theme_minimal()


mod_fec_s <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
                data = df_fecu_mod_site)
summary(mod_fec_s)


ggplot(df_fecu_mod_site, aes(y = nr_recr, x = nr_fru_t2)) +
  geom_jitter(width = .25, height = .5) +
  geom_smooth(method = 'lm')


# Fec- year Level --------------------------------------------------------------
df_fecu_mod_year <- df_fecu_mod %>%
  group_by(year) %>% 
  summarise(
    tot_p_area = sum(tot_p_area, na.rm = TRUE),
    g_cov      = sum(g_cov, na.rm = TRUE),
    nr_recr    = sum(nr_recr, na.rm = TRUE),
    nr_fru_t1  = sum(nr_fru_t1, na.rm = TRUE),
    nr_fru_t2  = sum(nr_fru_t2, na.rm = TRUE),
    nr_fru_t3  = sum(nr_fru_t3, na.rm = TRUE),
    nr_fru_t4  = sum(nr_fru_t4, na.rm = TRUE))


ggplot(df_fecu_mod_year %>%
         pivot_longer(
           cols = starts_with('nr_fru_t'),
           names_to = 'fruit_type',
           values_to = 'fruit_count'
         ), aes(x = fruit_count, y = nr_recr, color = fruit_type)) +
  geom_smooth(se = FALSE, method = 'lm') +
  geom_jitter(aes(shape = fruit_type)) +
  labs(
    title = 'Relationship between Fruit Counts and Recruits',
    x = 'Fruit Count',
    y = 'Number of Recruits',
    color = 'Fruit Type'
  ) +
  theme_minimal()


mod_fec_s <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
                data = df_fecu_mod_year)
summary(mod_fec_s)


ggplot(df_fecu_mod_year, aes(y = nr_recr, x = nr_fru_t2)) +
  geom_jitter(width = .25, height = .5) +
  geom_smooth(method = 'lm')



# Per capita recruits ----------------------------------------------------------
#  nr recruits t1 / nr reproductive indv t0
#  maybe on the single quadrat

df_rec_pc <- df_mean %>% 
  group_by(site, quad_id, year) %>% 
  summarise(rec_nr_t1 = sum(recruit, na.rm = T)) %>% 
  left_join(df_mean %>%
              group_by(site, quad_id, plant_id, year) %>%
              summarise(rep_nr_t0 = as.integer(any(fruit > 0)), 
                        .groups = 'drop') %>% 
              ungroup() %>% 
              group_by(site, quad_id, year) %>% 
              summarise(rep_nr_t0 = sum(rep_nr_t0, na.rm = T)) %>% 
              mutate(year = year + 1)
  , by = c('site', 'quad_id', 'year')) %>% 
  mutate(
    rec_pc = rec_nr_t1 / rep_nr_t0,
    fire_year = ifelse(year %in% c(2005, 2009, 2014, 2016, 2017), 'fire', 'normal'),
    rec_pc = ifelse(is.nan(rec_pc) | is.infinite(rec_pc), 0, rec_pc),  
    fire_year = ifelse(year %in% c(2006, 2010, 2015, 2017, 2018), 'fire_t1', fire_year),
    fire_year = ifelse(year %in% c(2017), 'FIRE', fire_year)
  )
  

ggplot(data = df_rec_pc, aes(y = rec_nr_t1, x = rep_nr_t0)) +
  geom_point() + 
  facet_wrap('quad_id', scales = 'free') +
  geom_smooth(method = 'lm') +
  theme_minimal()

ggplot(data = df_rec_pc, aes(y = rec_nr_t1, x = rep_nr_t0)) +
  geom_jitter(width = 0.05) + 
  geom_smooth(method = 'lm') +
  facet_wrap('quad_id', scales = 'free') +
  theme_minimal()

ggplot(data = df_rec_pc, aes(y = rec_nr_t1, x = rep_nr_t0)) +
  geom_jitter(width = 0.05) + 
  geom_smooth(method = 'lm') + 
  facet_wrap('site') +
  theme_minimal()


df_rec_pc <- df_rec_pc %>%
  mutate(
    rec_pc = rec_nr_t1 / rep_nr_t0,
    fire_year = ifelse(year %in% c(2005, 2009, 2014, 2016, 2017), 'fire', 'normal'),
    rec_pc = ifelse(is.nan(rec_pc) | is.infinite(rec_pc), 0, rec_pc),  # <- key fix
    fire_year = ifelse(year %in% c(2006, 2010, 2015, 2017, 2018), 'fire_t1', fire_year),
    fire_year = ifelse(year %in% c(2017), 'FIRE', fire_year)
  )


ggplot(df_rec_pc, aes(x = year, y = rec_pc, color = fire_year)) +
  geom_point() +
  facet_wrap('quad_id') +
  scale_color_manual(values = c('fire' = 'red', 'normal' = 'gray', 'fire_t1' = 'pink', 'FIRE' = 'purple')) +
  theme_minimal() +
  labs(y = expression('Per capita recruits'),
       x = expression('Year'))


# Flowering data ---------------------------------------------------------------
df_flower <- df_mean %>% 
  filter(!is.na(flower)) %>%
  filter(size_t0 != 0) %>%
  select(plant_id, year, size_t0, flower, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>% 
  mutate(flower = if_else(flower > 0, 1, flower))

g_flower_overall <- ggplot(
  data = plot_binned_prop(df_flower, 10, logsize_t0, flower)) +
  geom_point(aes(x = logsize_t0,
                 y = flower),
             alpha = 1, pch = 16, color = 'red' ) +
  geom_errorbar(aes(x = logsize_t0, ymin = lwr, ymax = upr),
                size = 0.5, width = 0.5) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1.01) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('log(size)'[t0]),
       y        = expression('Flowering to time t1')) +
  theme(plot.subtitle = element_text(size = 8))
g_flower_overall

ggplot(data = df_flower) +
  geom_jitter(aes(x = logsize_t0, y = flower), 
              position = position_jitter(width = 0.1, height = 0.3)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Flowering',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('Flowering to time t1')) +
  theme(plot.subtitle = element_text(size = 8))



# Survival model ---------------------------------------------------------------
# Logistic regression
mod_su_0 <- glm(survives ~ 1,
                data = df_surv, family = 'binomial') 
# Logistic regression
mod_su_1 <- glm(survives ~ logsize_t0,
                data = df_surv, family = 'binomial') 
# Quadratic logistic model
mod_su_2 <- glm(survives ~ logsize_t0 + logsize_t0_2,
                data = df_surv, family = 'binomial')  
# Cubic logistic model
mod_su_3 <- glm(survives ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                data = df_surv, family = 'binomial')  


# Compare models using AIC
mods_su      <- list(mod_su_0, mod_su_1, mod_su_2, mod_su_3)
mods_su_dAIC <- AICtab(mods_su, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_su_sorted <- order(mods_su_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_su) == 0) {
  mod_su_index_bestfit <- mods_su_sorted[1]
  v_mod_su_index       <- mod_su_index_bestfit - 1 
} else {
  mod_su_index_bestfit <- v_mod_set_su +1
  v_mod_su_index       <- v_mod_set_su
}

mod_su_bestfit   <- mods_su[[mod_su_index_bestfit]]
mod_su_ranef         <- coef(mod_su_bestfit)

# Generate predictions for survival across a range of sizes
mod_su_x <- seq(
  min(df_surv$logsize_t0, na.rm = T),
  max(df_surv$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_surv_pred <- predictor_fun(mod_su_x, mod_su_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_su_x, survives = .)

g_surv_line <- ggplot() +
  geom_jitter(data = df_surv, aes(x = logsize_t0, 
                                  y = survives),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_surv_pred, aes(x = logsize_t0, 
                                     y = survives),
            color = line_color_pred_fun(mod_su_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_surv_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_mean, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        y = survives) ) +
  geom_errorbar(
    data = plot_binned_prop(df_mean, 10, logsize_t0, survives), 
    aes(x = logsize_t0, 
        ymin = lwr,
        ymax = upr) ) +
  geom_line(data = df_surv_pred, aes(x = logsize_t0, 
                                     y = survives),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
g_surv_overall_pred <- g_surv_line + g_surv_bin + plot_layout()
g_surv_overall_pred


# Growth model -----------------------------------------------------------------
mod_gr_0 <- lm(logsize_t1 ~ 1, 
               data = df_grow)
# Linear model
mod_gr_1   <- lm(logsize_t1 ~ logsize_t0, 
                 data = df_grow)
# Quadratic model
mod_gr_2 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2, 
               data = df_grow)  
# Cubic model
mod_gr_3 <- lm(logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3, 
               data = df_grow)

mods_gr      <- list(mod_gr_0, mod_gr_1, mod_gr_2, mod_gr_3)
mods_gr_dAIC <- AICtab(mods_gr, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_gr_sorted <- order(mods_gr_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_gr) == 0) {
  mod_gr_index_bestfit <- mods_gr_sorted[1]
  v_mod_gr_index       <- mod_gr_index_bestfit - 1 
} else {
  mod_gr_index_bestfit <- v_mod_set_gr +1
  v_mod_gr_index       <- v_mod_set_gr
}

mod_gr_bestfit         <- mods_gr[[mod_gr_index_bestfit]]
mod_gr_ranef           <- coef(mod_gr_bestfit)

# Predict size at time t1 using the mean growth model
df_grow$pred <- predict(mod_gr_bestfit, type = 'response')

source('helper_functions/line_color_pred_fun.R')
source('helper_functions/predictor_fun.R')

g_grow_line <- ggplot(
  df_grow, aes(x = logsize_t0, y = logsize_t1)) +
  # Plot observed data
  geom_point() +
  geom_function(fun = function(x) predictor_fun(x, mod_gr_ranef), 
                color = line_color_pred_fun(mod_gr_ranef), 
                lwd = 2) +
  theme_bw() + 
  labs(title    = 'Growth prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_grow_pred <- ggplot(
  df_grow, aes(x = pred, y = logsize_t1)) +
  geom_point() +  
  geom_abline(aes(intercept = 0, slope = 1),  
              color = 'red', lwd = 2) + 
  theme_bw()

g_grow_overall_pred <- g_grow_line + g_grow_pred + plot_layout() 
g_grow_overall_pred


# Fit a model to assess variance in growth
# Fitted values from growth model
mod_gr_x   <- fitted(mod_gr_bestfit)  
# Squared residuals
mod_gr_y   <- resid(mod_gr_bestfit)^2  
# Non-linear model for variance
mod_gr_var <- nls(
  mod_gr_y ~ a * exp(b * mod_gr_x), start = list(a = 1, b = 0),
  control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE) ) 


# Recruitment model ------------------------------------------------------------
df_recr_nona_nr_quad <- df_recr %>% filter(!is.na(nr_quad))
# Fit a negative binomial model for recruitment
mod_rec <- MASS::glm.nb(nr_quad ~ 1, data = df_recr_nona_nr_quad)

# Generate predictions for recruitment
df_recr_nona_nr_quad <- df_recr_nona_nr_quad %>% 
  mutate(mod_pred = predict(mod_rec, type = 'response')) 

# Summarize total number of recruits and predictions
df_rec_sums_m <- df_recr_nona_nr_quad %>%
  summarize(nr_quad = sum(nr_quad),
            mod_pred = sum(mod_pred))

# Count number of adult individuals
indiv_m <- df_surv %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
repr_pc_m <- indiv_m %>%
  bind_cols(df_rec_sums_m) %>%
  mutate(repr_pc_mean = mod_pred / n_adults) %>%
  mutate(repr_pc_obs = nr_quad / n_adults) %>%
  drop_na 


# Fruiting model ----------------------------------------------------------------
df_fruit_nona_nr_quad <- df_fruit %>% filter(!is.na(nr_quad))
# Fit a negative binomial model
mod_fru <- MASS::glm.nb(nr_quad ~ 1, data = df_fruit_nona_nr_quad)

# Generate predictions
df_fruit_nona_nr_quad <- df_fruit_nona_nr_quad %>% 
  mutate(mod_pred = predict(mod_fru, type = 'response')) 

# Summarize total number of fruits and predictions
df_fruit_sums_m <- df_fruit_nona_nr_quad %>%
  summarize(nr_quad = sum(nr_quad),
            mod_pred = sum(mod_pred))

# Count number of adult individuals
indiv_m <- df_surv %>%
  summarize(n_adults = n())

# Calculate reproduction per capita (both observed and predicted)
fruit_pc_m <- indiv_m %>%
  bind_cols(df_fruit_sums_m) %>%
  mutate(fruit_pc_mean = mod_pred / n_adults) %>%
  mutate(fruit_pc_obs  = nr_quad  / n_adults) %>%
  drop_na 


# Flowering model --------------------------------------------------------------
# Logistic regression
mod_fl_0 <- glm(flower ~ 1,
                data = df_flower, family = 'binomial') 
# Logistic regression
mod_fl_1 <- glm(flower ~ logsize_t0,
                data = df_flower, family = 'binomial') 
# Quadratic logistic model
mod_fl_2 <- glm(flower ~ logsize_t0 + logsize_t0_2,
                data = df_flower, family = 'binomial')  
# Cubic logistic model
mod_fl_3 <- glm(flower ~ logsize_t0 + logsize_t0_2 + logsize_t0_3,
                data = df_flower, family = 'binomial')  


# Compare models using AIC
mods_fl      <- list(mod_fl_0, mod_fl_1, mod_fl_2, mod_fl_3)
mods_fl_dAIC <- AICtab(mods_fl, weights = T, sort = F)$dAIC

# Get the sorted indices of dAIC values
mods_fl_sorted <- order(mods_fl_dAIC)

# Establish the index of model complexity
if (length(v_mod_set_fl) == 0) {
  mod_fl_index_bestfit <- mods_fl_sorted[1]
  v_mod_fl_index       <- mod_fl_index_bestfit - 1 
} else {
  mod_fl_index_bestfit <- v_mod_set_fl +1
  v_mod_fl_index       <- v_mod_set_fl
}


mod_fl_bestfit   <- mods_fl[[mod_fl_index_bestfit]]
mod_fl_ranef     <- coef(mod_fl_bestfit)


# Generate predictions for survival across a range of sizes
mod_fl_x <- seq(
  min(df_flower$logsize_t0, na.rm = T),
  max(df_flower$logsize_t0, na.rm = T), length.out = 100)

# Prepare data for survival plot
df_flow_pred <- predictor_fun(mod_fl_x, mod_fl_ranef) %>% 
  # Inverse logit for predictions
  boot::inv.logit() %>% 
  data.frame(logsize_t0 = mod_fl_x, flower = .)

g_flow_line <- ggplot() +
  geom_jitter(data = df_flower, 
              aes(x = logsize_t0, y = flower),
              alpha = 0.25, width = 0.08, height = 0.3) +
  geom_line(data = df_flow_pred, 
            aes(x = logsize_t0, y = flower),
            color = line_color_pred_fun(mod_fl_ranef), 
            lwd   = 2) +  
  theme_bw() + 
  labs(title    = 'Survival prediction',
       subtitle = v_ggp_suffix) +
  theme(plot.subtitle = element_text(size = 8))

g_flow_bin <- ggplot() +
  geom_point(data =  plot_binned_prop(
    df_flower, 10, logsize_t0, flower), 
    aes(x = logsize_t0, y = flower) ) +
  geom_errorbar(
    data = plot_binned_prop(df_flower, 10, logsize_t0, flower), 
    aes(x = logsize_t0, ymin = lwr, ymax = upr) ) +
  geom_line(data = df_flow_pred, 
            aes(x = logsize_t0, y = flower),
            color = 'red', lwd   = 2) + 
  theme_bw() +
  ylim(0, 1)

# Combine survival plots
g_flow_overall_pred <- g_flow_line + g_flow_bin + plot_layout()
g_flow_overall_pred



# Exporting parameter estimates ------------------------------------------------
# Growth
grow_fe  <- data.frame(coefficient = names(coef(mod_gr_bestfit)),
                       value       = coef(mod_gr_bestfit))
grow_var <- data.frame(coefficient = names(coef(mod_gr_var)),
                       value       = coef(mod_gr_var))

grow_out <- Reduce(function(...) rbind(...), list(grow_fe, grow_var)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# write.csv(grow_out, row.names = F, paste0(
#   dir_data, '/', v_script_prefix, '_', v_sp_abb, '_grow_pars_mean.csv'))


# Survival
surv_fe  <- data.frame(coefficient = names(coef(mod_su_bestfit)),
                       value       = coef(mod_su_bestfit))

surv_out <- Reduce(function(...) rbind(...), list(surv_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# write.csv(surv_out, row.names = F, paste0(
#   dir_data, '/', v_script_prefix, '_', v_sp_abb, '_surv_pars_mean.csv'))


# Flower
flwr_fe  <- data.frame(coefficient = names(coef(mod_fl_bestfit)),
                       value       = coef(mod_fl_bestfit))

flwr_out <- Reduce(function(...) rbind(...), list(flwr_fe)) %>%
  mutate(coefficient = as.character(coefficient)) %>%
  mutate(coefficient = replace(
    coefficient, grepl('Intercept', coefficient), 'b0'))

# write.csv(flwr_out, row.names = F, paste0(
#   dir_data, '/', v_script_prefix, '_', v_sp_abb, '_surv_pars_mean.csv'))

# Recruitment 
rec_size <- df_mean %>% subset(recruit == 1)
fru_size <- df_mean %>% subset(fruit   == 1)

others   <- data.frame(coefficient = c('rec_siz', 'rec_sd', 
                                       'fru_siz', 'fru_sd', 
                                       'max_siz', 'min_siz',
                                       'fecu_b0',
                                       'frui_b0'),
                       value       = c(mean(log(rec_size$size_t0), na.rm = T), 
                                       sd(  log(rec_size$size_t0), na.rm = T),
                                       mean(log(fru_size$size_t0), na.rm = T), 
                                       sd(  log(fru_size$size_t0), na.rm = T),
                                       df_grow$logsize_t0 %>% max, 
                                       df_grow$logsize_t0 %>% min,
                                       repr_pc_m$repr_pc_mean,
                                       fruit_pc_m$fruit_pc_mean))

# write.csv(others, row.names = F, paste0(
#   dir_data, '/', v_script_prefix, '_', v_sp_abb, '_other_pars_mean.csv'))



# Building the IPM from scratch ------------------------------------------------
extr_value <- function(x, field){
  subset(x, coefficient == field)$value
}

pars <- Filter(function(x) length(x) > 0, list(
  prefix  = v_script_prefix,
  species = v_species,
  grow_b0 = extr_value(grow_out, 'b0'),
  grow_b1 = extr_value(grow_out, 'logsize_t0'),
  grow_b2 = extr_value(grow_out, 'logsize_t0_2'),
  grow_b3 = extr_value(grow_out, 'logsize_t0_3'),
  a       = extr_value(grow_out, 'a'),
  b       = extr_value(grow_out, 'b'),
  surv_b0 = extr_value(surv_out, 'b0'),
  surv_b1 = extr_value(surv_out, 'logsize_t0'),
  surv_b2 = extr_value(surv_out, 'logsize_t0_2'),
  surv_b3 = extr_value(surv_out, 'logsize_t0_3'),
  flwr_b0 = extr_value(flwr_out, 'b0'),
  flwr_b1 = extr_value(flwr_out, 'logsize_t0'),
  flwr_b2 = extr_value(flwr_out, 'logsize_t0_2'),
  flwr_b3 = extr_value(flwr_out, 'logsize_t0_3'),
  fecu_b0 = extr_value(others, 'fecu_b0'),
  recr_sz = extr_value(others, 'rec_siz'),
  recr_sd = extr_value(others, 'rec_sd'),
  frui_b0 = extr_value(others, 'frui_b0'),
  frui_sz = extr_value(others, 'fru_siz'),
  frui_sd = extr_value(others, 'fru_sd'),
  L       = extr_value(others, 'min_siz'),
  U       = extr_value(others, 'max_siz'),
  mat_siz = 200,
  mod_gr_index = v_mod_gr_index,
  mod_su_index = v_mod_su_index
))

# write.csv(pars, row.names = F, paste0(
#   dir_data, '/', v_script_prefix, '_', v_sp_abb, '_pars.csv'))

# Function describing standard deviation of growth model
grow_sd <- function(x, pars) {
  pars$a * (exp(pars$b* x)) %>% sqrt 
}

# Growth from size x to size y
gxy <- function(x, y, pars, num_pars = v_mod_gr_index) {
  mean_value <- 0
  for (i in 0:num_pars) {
    param_name <- paste0('grow_b', i)
    if (!is.null(pars[[param_name]])) {
      mean_value <- mean_value + pars[[param_name]] * x^i
    }
  }
  sd_value <- grow_sd(x, pars)
  return(dnorm(y, mean = mean_value, sd = sd_value))
}


# Function describing the invert logit
inv_logit <- function(x) {exp(x) / (1 + exp(x))}


# Survival of x-sized individual to time t1
sx <- function(x, pars, num_pars = v_mod_su_index) {
  survival_value <- pars$surv_b0
  for (i in 1:num_pars) {
    param_name <- paste0('surv_b', i)
    if (!is.null(pars[[param_name]])) {
      survival_value <- survival_value + pars[[param_name]] * x^(i)
    }
  }
  return(inv_logit(survival_value))
}

# Function describing the transition kernel
pxy <- function(x, y, pars) {
  return(sx(x, pars) * gxy(x, y, pars))
}

# Function describing the recruitment 
fy <- function(y, pars, h){
  n_recr  <- pars$fecu_b0
  recr_y  <- dnorm(y, pars$recr_sz, pars$recr_sd) * h
  recr_y  <- recr_y / sum(recr_y)
  f       <- n_recr * recr_y
  return(f)
}


# Kernel
kernel <- function(pars) {
  
  # number of bins over which to integrate
  n   <- pars$mat_siz 
  # lower limit of integration
  L   <- pars$L  
  # upper limit of integration
  U   <- pars$U       
  # bin size
  h   <- (U - L) / n  
  # lower boundaries of bins
  b   <- L + c(0:n) * h             
  # midpoints of bins
  y   <- 0.5 * (b[1:n] + b[2:(n + 1)]) 
  
  # Fertility matrix
  Fmat        <- matrix(0, n, n)
  Fmat[]      <- matrix(fy(y, pars, h), n, n)
  
  # Survival vector
  Smat   <- c()
  Smat   <- sx(y, pars)
  
  # Growth matrix
  Gmat   <- matrix(0, n, n)
  Gmat[] <- t(outer(y, y, gxy, pars)) * h
  
  # Growth/survival transition matrix
  Tmat   <- matrix(0, n, n)
  
  # Correct for eviction of offspring
  for(i in 1:(n / 2)) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for(i in (n / 2 + 1):n) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum(Gmat[,i])
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrices
  k_yx <- Fmat + Tmat
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y))
}

lambda_ipm <- function(i) {
  return(Re(eigen(kernel(i)$k_yx)$value[1]))
}

# mean population growth rate
lam_mean <- lambda_ipm(pars)
lam_mean
