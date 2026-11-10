# Packages --------------------------------------------------------------------
library(tidyverse)
library(Rpadrino)

# Life history functions ------------------------------------------------------
source("uncertainty_life_history/life_history_metrics.R")
source("uncertainty_life_history/padrino_lh_test.R")

# IPM objects -----------------------------------------------------------------
ipm_ange <- readRDS("adler_2007_ks/data/ange/adler07_ange_mean_ipmr.rds")


# Andropogon gerardii ---------------------------------------------------------
ange_record <- make_ipm_record(
  ipm = ipm_ange,
  ipm_id = "adler07_ange_mean",
  species = "Andropogon gerardii",
  P_kernels = "P",
  F_kernels = "F",
  timestep = 1,
  maturity_threshold = 0.1)


# Life history metrics --------------------------------------------------------
life_history_df <- calculate_life_history_metrics(ange_record)

life_history_df


# IPM records -----------------------------------------------------------------
ipm_records <- list(adler07_ange_mean = ange_record)

# Record directory ------------------------------------------------------------
dir_records <- "uncertainty_life_history/data"
file_records <- file.path(dir_records, "ipm_records.rds")

dir.create(dir_records, recursive = TRUE, showWarnings = FALSE)

# Save records ----------------------------------------------------------------
if (!file.exists(file_records)) {
  saveRDS(ipm_records, file_records)
}

# Load records ----------------------------------------------------------------
ipm_records <- readRDS(file_records)

names(ipm_records)


# Life history metrics --------------------------------------------------------
life_history_df <- calculate_life_history_metrics_batch(ipm_records)

life_history_df


# Life history results --------------------------------------------------------
file_metrics <- file.path(dir_records, "life_history_metrics.csv")

# Save results ----------------------------------------------------------------
if (!file.exists(file_metrics)) {
  write.csv(life_history_df, file_metrics, row.names = FALSE)
}

# Load results ----------------------------------------------------------------
life_history_df <- read.csv(file_metrics)

life_history_df
