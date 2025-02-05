library(sf)
library(dplyr)
library(tidyverse)


base_dir <- "C:/code/RUPDemo_IPMs/moore_2021_az/data/Species_Shapefile_Extractions/"

shapefiles <- list.files(base_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)


# Function to extract year and quadrat information from file names
extract_info <- function(file_path) {
  # Extract the year (assuming it's a 4-digit number in the filename)
  year <- sub(".*_(\\d{4})\\.shp$", "\\1", basename(file_path))
  
  # Extract quadrat identifier (assuming it's before the year in the filename)
  quadrat <- sub("^(.*)_(\\d{4})\\.shp$", "\\1", basename(file_path))
  
  return(c(quadrat = quadrat, year = year))  # Return as a vector
}

# Apply the function to all shapefiles and store results
file_info <- lapply(shapefiles, extract_info)

# Convert to a data frame for easier manipulation
file_info_df <- do.call(rbind, file_info)
colnames(file_info_df) <- c("quadrat", "year")
file_info_df <- as.data.frame(file_info_df)

plot(as.factor(file_info_df$quadrat), file_info_df$year)


file_info_df <- file_info_df %>% 
  mutate(year = str_sub(as.factor(year))) %>%
  rename(Year = year) %>% 
  mutate(Year = as.numeric(str_sub(Year, start = -2))) %>%
  mutate(quadrat = gsub("Quadrat_", "", quadrat),
         quadrat = gsub("_bar_", " / ", quadrat))

# Get a list of unique quadrats for each year
unique_quadrats_by_year <- file_info_df %>%
  group_by(Year) %>%
  summarise(unique_quadrats = list(unique(quadrat))) %>%
  arrange(Year)

# Get a list of unique year for each quadrats
unique_quadrats_by_plot <- file_info_df %>%
  group_by(quadrat) %>%
  summarise(Year = list(unique(Year))) %>%
  arrange(quadrat)



# Convert to a list
unique_quadrats_list <- setNames(unique_quadrats_by_year$unique_quadrats, unique_quadrats_by_year$Year)
unique_year_list     <- setNames(unique_quadrats_by_plot$Year           , unique_quadrats_by_plot$quadrat)
inv_sgs <- unique_year_list

# Print the result
print(unique_quadrats_list)




# Create a data frame with each year-quadrat combination
years <- rep(unique(file_info_df$year), each = length(unique(file_info_df$quadrat)))
quadrats <- rep(unique(file_info_df$quadrat), times = length(unique(file_info_df$year)))

# Create a data frame of quadrat-year combinations
plot_data <- data.frame(year = years, quadrat = quadrats)

# Create a logical column to indicate whether the quadrat exists for that year
plot_data$match <- ifelse(paste(plot_data$quadrat, plot_data$year) %in% paste(file_info_df$quadrat, file_info_df$year), 1, 0)

# View the first few rows of the prepared plot data
head(plot_data)




ggplot(plot_data, aes(x = year, y = quadrat, color = factor(match))) +
  geom_point(shape = 16, size = 3) +   # Shape 16 is a solid dot
  scale_color_manual(values = c("0" = "white", "1" = "blue")) +  # White for no match, blue for match
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x labels for readability
        axis.title = element_blank(),  # Remove axis titles
        legend.position = "none") +  # Remove the legend
  labs(title = "Presence of Quadrats by Year")




#### ----------








