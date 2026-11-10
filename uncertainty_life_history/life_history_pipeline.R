
# Packages --------------------------------------------------------------------
library(tidyverse)
library(Rpadrino)

# Life history functions ------------------------------------------------------
source("uncertainty_life_history/life_history_metrics.R")
source("uncertainty_life_history/padrino_lh_test.R")

# Life history pipeline -------------------------------------------------------
run_lh_pipeline <- function(models, overwrite = FALSE) {
  
  # Record directory ----------------------------------------------------------
  dir_records <- "uncertainty_life_history/data"
  file_records <- file.path(dir_records, "ipm_records.rds")
  file_metrics <- file.path(dir_records, "life_history_metrics.csv")
  
  dir.create(dir_records, recursive = TRUE, showWarnings = FALSE)
  
  # Load existing records -----------------------------------------------------
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
  
  update_records <- FALSE
  update_metrics <- FALSE
  
  # Process IPMs --------------------------------------------------------------
  for (path in models) {
    
    model_id <- sub("_ipmr\\.rds$", "", basename(path))
    
    # Update IPM records ------------------------------------------------------
    if (overwrite || !(model_id %in% names(ipm_records))) {
      
      ipm <- readRDS(path)
      
      species <- attr(ipm, "species")
      if (is.null(species)) {
        species <- ipm_records[[model_id]]$species
      }
      
      record <- make_ipm_record(
        ipm = ipm,
        ipm_id = model_id,
        species = species,
        P_kernels = "P",
        F_kernels = "F",
        timestep = 1,
        maturity_threshold = 0.1)
      
      ipm_records[[model_id]] <- record
      update_records <- TRUE
    }
    
    # Update life history metrics ---------------------------------------------
    if (overwrite || !(model_id %in% life_history_df$ipm_id)) {
      
      new_metrics <- calculate_life_history_metrics(
        ipm_records[[model_id]])
      
      life_history_df <- life_history_df %>%
        filter(ipm_id != model_id) %>%
        bind_rows(new_metrics)
      
      update_metrics <- TRUE
    }
  }
  
  # Save updated records ------------------------------------------------------
  if (update_records) {
    saveRDS(ipm_records, file_records)
  }
  
  if (update_metrics) {
    write.csv(life_history_df, file_metrics, row.names = FALSE)
  }
  
  return(life_history_df)
}
