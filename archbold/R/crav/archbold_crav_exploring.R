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
load_packages(tidyverse, patchwork, skimr, ipmr, binom, bbmle)


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


# Data -------------------------------------------------------------------------
df <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data.csv')) %>% 
  janitor::clean_names() %>%
  mutate(across(c(5, 6, 11:length(.)), ~ na_if(., 9999))) %>%
  mutate(
    plant_id = paste(site, quad, plant, sep = "_"),
    year     = substr(date, 1, 4),  
    month    = substr(date, 6, 7)) %>%
  arrange(site, quad, plant, year, month)

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
    'plant identification', 'sample year', 'sample month'))

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
  semi_join(odd_su_56_uniq, by = c("site", "quad", "plant"))


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


# Original year mean dataframe
df_mean_og <- df %>%
  filter(s < 6 | is.na(s)) %>%
  group_by(site, quad, plant, plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(s )), NA_real_, max(s,  na.rm = TRUE)),
    size_t0 =  if_else(all(is.na(br)), NA_real_, max(br, na.rm = TRUE)),
    .groups = "drop"
  ) %>% 
  ungroup()



# Base mean dataframe
df_mean <- df_mean_og %>% 
  group_by(site, quad, plant, plant_id) %>% 
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
  
  ungroup()

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
  geom_histogram(binwidth = 1, fill = "lightgray", color = "black") +
  theme_minimal()

# Size over age of individuals with tracked recruit 
ggplot(df_mean, aes(x = age, y = logsize_t0)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +  
  theme_minimal() + 
  labs(x = "Age", y = "Log Size", title = "Log Size vs Age of Recruits") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, by = 1))


# Histogram of size after leaving dormancy
ggplot(df_mean %>% filter(dormancy == 1), aes(x = size_t1)) +
  geom_histogram(binwidth = 1, fill = "lightgray", color = "black") +
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


# Histogram of recruit size (at "first encounter")
ggplot(df_mean_f %>% filter(age == 1), aes(x = size_t0)) +
  geom_histogram(binwidth = 1, fill = "lightgray", color = "black") +
  theme_minimal()

# 
ggplot() +
  geom_jitter(
    data = df_mean_f, 
    aes(x = age, y = logsize_t0)) +
  theme_minimal()


# Load necessary libraries
library(lme4)
library(ggplot2)

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
  labs(title = "Linear Model with Random Intercepts: logsize_t0 vs Age",
       y = "Age",
       x = "Log Size at t0") +
  theme(legend.position = "none")

# Plot the second model (quadratic relationship for logsize_t0 ~ logsize_t0^2)
ggplot(df_mean_f, aes(x = logsize_t0, y = age)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_quad, color = plant_id), alpha = 0.5) +  
  theme_minimal() +
  labs(title = "Quadratic Model with Random Intercepts: logsize_t0 vs Age",
       y = "Age",
       x = "Log Size at t0") +
  theme(legend.position = "none")

# Plot the third model (logsize_t0 ~ age relationship)
ggplot(df_mean_f, aes(x = age, y = logsize_t0)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_logsize_age, color = plant_id), alpha = 0.5) +  
  theme_minimal() +
  labs(title = "Log Size Model with Random Intercepts: Age vs Log Size at t0",
       y = "Log Size at t0",
       x = "Age") +
  theme(legend.position = "none")

# Plot the fourth model (logsize_t0 ~ age + I(age^2) relationship)
ggplot(df_mean_f, aes(x = age, y = logsize_t0)) +
  geom_jitter(aes(color = plant_id), width = 0.1, height = 0.1, alpha = 0.6) +  
  geom_line(aes(y = predicted_quad_age, color = plant_id), alpha = 0.5) + 
  theme_minimal() +
  labs(title = "Quadratic Model with Random Intercepts: Age vs Log Size at t0",
       y = "Log Size at t0",
       x = "Age") +
  theme(legend.position = "none")










# Individuals to check out because of certains flags ---------------------------
#   1, 1, 23, 2017
# size t0 after recruit
#   1, 1, 23, 2012, 1, 14
#   1, 1, 37, 2014, 1, 13
#   1, 4, 17, 2012, 1, NA
#   1, 2, 22, 2010, 1, 19, 1, 8


# Model df's -------------------------------------------------------------------
# Survival data frame
surv_df <- df_mean_dor %>% 
  filter(!is.na(survives))
  # subset(df_mean_dor, !is.na(survives)) %>%
  # subset(size_t0 != 0) %>%
  # select(quad, track_id, year, size_t0, survives, size_t1,
  #        logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Growth data frame
grow_df <- df_mean_dor %>% 
  filter(!is.na(size_t0))

# # Total area data frame
# quad_df <- df %>%  
#   group_by (quad, year) %>% 
#   summarise(tot_p_area = sum(size_t0, na.rm = T)) %>% 
#   ungroup
# 
# group_df <- quad_df %>% 
#   group_by (year) %>% 
#   summarise(g_cov = mean(tot_p_area)) %>% 
#   ungroup
# 
# cover_df <- left_join(quad_df, group_df) %>%
#   mutate(year = year + 1) %>% 
#   mutate(year = as.integer(year)) %>% 
#   drop_na()
# 
# # Recruitment data frame
# recr_df <- df %>%
#   group_by (year, quad) %>% 
#   summarise(nr_quad = sum(recruit, na.rm = T)) %>% 
#   ungroup
# 
# recr_df <- left_join(cover_df, recr_df)




##########################################
# Size t0 and t1 Histogram
g_hist_t0 <- ggplot(
  df_mean_dor, aes(x = logsize_t0)) +
  geom_histogram(binwidth = 0.2, fill = 'grey', color = 'black') +
  # labs(title    = 'Histogram',
  #      subtitle = v_ggp_suffix, 
  #      x        = 'Size at time t0') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(size = 8))
g_hist_t1 <- ggplot(
  df_mean_dor, aes(x = logsize_t1)) +
  geom_histogram(binwidth = 0.2, fill = 'white', color = 'black') +
  labs(x = 'Size at time t1') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

g_hist_sizes_log <- g_hist_t0 + g_hist_t1


# Graphical exploration - Growth
g_gr_overall <- ggplot(
  data  = grow_df, aes( x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  # labs(title    = 'Growth',
  #      subtitle = v_ggp_suffix,
  #      x        = expression('log(size) ' [t0]),
  #      y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8))



# Survival analysis
# Load custom function for plotting binned proportions
source('helper_functions/plot_binned_prop.R')
# plot_binned_prop(df_mean_dor, 10, size_t0, survives)

# Generate a plot for overall survival based on size at t0
g_surv_overall <- ggplot(
  data = df_mean) +
  geom_point(aes(x = size_t0, 
                 y = survives),
             alpha = 1, pch = 16, color = 'red' ) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9)) +
  ylim(0, 1) +
  theme_bw()


