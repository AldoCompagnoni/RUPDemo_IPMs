# Packages --------------------------------------------------------------------
library(tidyverse)
library(Rpadrino)

# Life history functions ------------------------------------------------------
source("uncertainty_life_history/life_history_metrics.R")
source("uncertainty_life_history/padrino_lh_test.R")

# IPM objects -----------------------------------------------------------------
ipm_artr <- readRDS("zachmann_2016_id/data/artr/zachmann16_artr_mean_ipmr.rds")


# Andropogon gerardii ---------------------------------------------------------
artr_record <- make_ipm_record(
  ipm = ipm_artr,
  ipm_id = "zachmann16_artr_mean",
  species = "Artemisia tripartita",
  P_kernels = "P",
  F_kernels = "F",
  timestep = 1,
  maturity_threshold = 0.1)


# Record storage --------------------------------------------------------------
overwrite <- FALSE
record <- artr_record

dir_records <- "uncertainty_life_history/data"
file_records <- file.path(dir_records, "ipm_records.rds")
file_metrics <- file.path(dir_records, "life_history_metrics.csv")

dir.create(dir_records, recursive = TRUE, showWarnings = FALSE)

# Load existing records -------------------------------------------------------
ipm_records <- if (file.exists(file_records)) {
  readRDS(file_records)
} else {
  list()
}

life_history_df <- if (file.exists(file_metrics)) {
  read.csv(file_metrics)
} else {
  tibble(ipm_id = character(), species = character())
}

# Update IPM records ----------------------------------------------------------
model_id <- record$ipm_id

if (overwrite || !(model_id %in% names(ipm_records))) {
  ipm_records[[model_id]] <- record
  saveRDS(ipm_records, file_records)
}

# Update life history metrics -------------------------------------------------
if (overwrite || !(model_id %in% life_history_df$ipm_id)) {
  new_metrics <- calculate_life_history_metrics(ipm_records[[model_id]])
  
  life_history_df <- life_history_df %>%
    filter(ipm_id != model_id) %>%
    bind_rows(new_metrics)
  
  write.csv(life_history_df, file_metrics, row.names = FALSE)
}

# Check stored models ---------------------------------------------------------
names(ipm_records)
life_history_df %>% select(ipm_id, species)


# Test stored records ---------------------------------------------------------
ipm_records <- readRDS(
  "uncertainty_life_history/data/ipm_records.rds")

life_history_df <- calculate_life_history_metrics_batch(ipm_records)

life_history_df
