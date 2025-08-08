library(tidyverse)

# # Download all data ------------------------------------------------------------
# # Set the root folder
# {root_folder <- "archbold/R"
# 
# # Find all matching R scripts recursively
# script_files <- list.files(
#   path = root_folder,
#   pattern = "^archbold_[A-Za-z0-9]{4}_download\\.R$",
#   recursive = TRUE,
#   full.names = TRUE
# )
# 
# # Run each script using source()
# for (script in script_files) {
#   message("Running: ", script)
#   source(script)
# }}

# Read the data ----------------------------------------------------------------
csv_files <- list.files(path = "archbold/data", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
data_files <- grep("_data(_v2)?\\.csv$", csv_files, value = TRUE)

df_erlo <- read_csv(data_files[3]) %>% 
  janitor::clean_names() %>%
  separate(date, into = c("year", "month"), sep = "-", convert = TRUE) %>% 
  mutate(species = "erlo")

df_hycu <- read_csv(data_files[4]) %>% 
  janitor::clean_names() %>% 
  mutate(species = "hycu") %>% 
  rename(yr_encounter = year,
         year = time)

df_lioh <- read_csv(data_files[5]) %>% 
  janitor::clean_names() %>% 
  mutate(species = "lioh")

df_pole <- read_csv(data_files[6]) %>% 
  janitor::clean_names() %>%
  separate(date, into = c("year", "month"), sep = "-", convert = TRUE) %>% 
  mutate(species = "pole")

df_sood <- read_csv(data_files[7]) %>% 
  janitor::clean_names() %>% 
  mutate(species = "sood")

# ------
# Count number of observations per year per species
bind_rows(df_erlo %>% head, df_sood %>% head())
yearly_counts <- bind_rows(df_erlo, df_hycu, df_lioh, df_pole, df_sood) %>%
  group_by(species, year) %>%
  summarise(n = n(), .groups = "drop")

ggplot(yearly_counts, aes(x = as.numeric(year), y = n, color = species)) +
  geom_line(size = 1) +
  geom_point() +
  labs(
    title = "Number of Observations per Year per Species",
    x = "Year",
    y = "Number of Observations",
    color = "Species") +
  ylim(0,10000) +
  theme_minimal()

yearly_counts %>% filter(species == 'sood')




