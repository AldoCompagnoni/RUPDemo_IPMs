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
# function to plot your survival data "binned" (instead of "jittered")
source('helper_functions/plot_binned_prop.R')


# Data -------------------------------------------------------------------------
df <- read_csv(file.path(dir_data, 'crotalaria_avonensis_data.csv')) %>% 
  janitor::clean_names() %>%
  mutate(across(c(5, 6, 11:length(.)), ~ na_if(., 9999))) %>%
  mutate(
    plant_id = as.factor(paste(site, quad, plant, sep = "_")),
    quad_id  = as.factor(paste(site, quad, sep = "_")),
    year     = as.numeric(substr(date, 1, 4)),  
    month    = as.numeric(substr(date, 6, 7)),
    site     = as.factor(site),
    quad     = as.factor(quad),
    mp       = as.factor(mp),
    plant    = as.factor(plant),
    caged    = as.factor(caged),
    veg      = as.factor(veg)) %>%
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


# Original year mean dataframe -------------------------------------------------
df_mean_og <- df %>%
  filter(s < 6 | is.na(s)) %>%
  group_by(site, quad, quad_id, plant, plant_id, year) %>%
  summarise(
    survives = if_else(all(is.na(s )), NA_real_, max (s,  na.rm = TRUE)),
    size_t0  = if_else(all(is.na(br)), NA_real_, max (br, na.rm = TRUE)),
    fruit    = if_else(all(is.na(fr)), NA_real_, mean(fr, na.rm = TRUE)),
    flower   = if_else(all(is.na(fl)), NA_real_, sum (fl, na.rm = TRUE)),
    fire_sev = if_else(
      all(is.na(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f))),
      NA_real_, 
      mean(c(burn_a, burn_b, burn_c, burn_d, burn_e, burn_f), na.rm = TRUE)),
    .groups  = "drop"
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



# Investigation: Fire ----------------------------------------------------------
names(df_mean)
df_fire1 <- df_mean %>% 
  group_by(quad_id, year) %>% 
  summarise(fire_sev = mean(fire_sev, na.rm = T))

ggplot(df_fire1, aes(x = quad_id, y = year, size = fire_sev, color = fire_sev)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Severity") +
  scale_color_continuous(name = "") +
  labs(x = "Quadrat ID", y = "Year", title = "Fire Severity by Year and Plot") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Fire on Growth
df_fire_grow <- df_mean %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(plant_id, year, size_t0, size_t1, age,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, fire_event)

ggplot(
  data  = df_fire_grow, aes(x = logsize_t0, y = logsize_t1)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  geom_abline(aes(slope = 1, intercept = 0)) + 
  geom_smooth(method = 'lm') +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        title     = element_text(size = 10)) +
  labs(title    = 'Fire on Growth',
       subtitle = v_ggp_suffix,
       x        = expression('log(size) ' [t0]),
       y        = expression('log(size)  '[t1])) +
  theme(plot.subtitle = element_text(size = 8)) +
  facet_wrap(~ fire_event)


# Fire on survival
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


# Fire on Recruits
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




# Fertility data ---------------------------------------------------------------

df_fert <- df_fruit %>%
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


ggplot(data = df_fert) +
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



df_fert_mod <- df_recr %>%
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


ggplot(df_fert_mod %>%
         pivot_longer(
           cols = starts_with("nr_fru_t"),
           names_to = "fru_type",
           values_to = "fru_count"
         ), 
       aes(x = fru_count, y = nr_recr, color = fru_type)) +
  geom_smooth(alpha = 0.2, method = 'lm') +
  geom_jitter(aes(shape = fru_type)) +
  labs(
    x = "Number of Fruits",
    y = "Number of Recruits",
    color = "Fruit Type",
    title = "Relationship between Fruit Counts and Recruits"
  ) +
  theme_minimal()


names(df_fert_mod)

mod_fer <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
              data = df_fert_mod)

summary(mod_fer)



# Site Level
df_fert_mod_site <- df_fert_mod %>%
  group_by(year, site) %>% 
  summarise(
    tot_p_area = sum(tot_p_area, na.rm = TRUE),
    g_cov      = sum(g_cov, na.rm = TRUE),
    nr_recr    = sum(nr_recr, na.rm = TRUE),
    nr_fru_t1  = sum(nr_fru_t1, na.rm = TRUE),
    nr_fru_t2  = sum(nr_fru_t2, na.rm = TRUE),
    nr_fru_t3  = sum(nr_fru_t3, na.rm = TRUE),
    nr_fru_t4  = sum(nr_fru_t4, na.rm = TRUE))


ggplot(df_fert_mod_site %>%
         pivot_longer(
           cols = starts_with("nr_fru_t"),
           names_to = "fruit_type",
           values_to = "fruit_count"
         ), aes(x = fruit_count, y = nr_recr, color = fruit_type)) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_jitter(aes(shape = fruit_type)) +
  labs(
    title = "Relationship between Fruit Counts and Recruits",
    x = "Fruit Count",
    y = "Number of Recruits",
    color = "Fruit Type"
  ) +
  theme_minimal()


mod_fer_s <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
                data = df_fert_mod_site)
summary(mod_fer_s)


ggplot(df_fert_mod_site, aes(y = nr_recr, x = nr_fru_t2)) +
  geom_jitter(width = .25, height = .5) +
  geom_smooth(method = 'lm')


# Year Level
df_fert_mod_year <- df_fert_mod %>%
  group_by(year) %>% 
  summarise(
    tot_p_area = sum(tot_p_area, na.rm = TRUE),
    g_cov      = sum(g_cov, na.rm = TRUE),
    nr_recr    = sum(nr_recr, na.rm = TRUE),
    nr_fru_t1  = sum(nr_fru_t1, na.rm = TRUE),
    nr_fru_t2  = sum(nr_fru_t2, na.rm = TRUE),
    nr_fru_t3  = sum(nr_fru_t3, na.rm = TRUE),
    nr_fru_t4  = sum(nr_fru_t4, na.rm = TRUE))


ggplot(df_fert_mod_year %>%
         pivot_longer(
           cols = starts_with("nr_fru_t"),
           names_to = "fruit_type",
           values_to = "fruit_count"
         ), aes(x = fruit_count, y = nr_recr, color = fruit_type)) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_jitter(aes(shape = fruit_type)) +
  labs(
    title = "Relationship between Fruit Counts and Recruits",
    x = "Fruit Count",
    y = "Number of Recruits",
    color = "Fruit Type"
  ) +
  theme_minimal()


mod_fer_s <- lm(nr_recr ~ nr_fru_t1 * nr_fru_t2 * nr_fru_t3, 
                data = df_fert_mod_year)
summary(mod_fer_s)


ggplot(df_fert_mod_year, aes(y = nr_recr, x = nr_fru_t2)) +
  geom_jitter(width = .25, height = .5) +
  geom_smooth(method = 'lm')



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



# Models -----------------------------------------------------------------------
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
