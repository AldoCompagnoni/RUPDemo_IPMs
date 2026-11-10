# Script to test life history metric functions on IPMs built with ipmr

# Relies on running the life_history_metrics.R script

# Models built with ipmr can include those already stored in PADRINO, or those
  # built by hand elsewhere. The goal is to slowly expand the set of models on
  # which we are testing, to ultimately build up to the entire set we plan to
  # use in the final analysis.

# =============================================================================
#
# Workflow:
#
#   1. Create made IPM objects with make_padrino_ipm() or manually with
#      make_ipm_record().
#
#   2. Store each made IPM together with the arguments needed to calculate
#      its life-history metrics.
#
#   3. Pass a list of these records to calculate_life_history_metrics().
#
#   4. Bind the resulting one-row data frames into a table.
#
# This allows models from PADRINO and independently constructed ipmr models
# to be analyzed with the same functions.
#
# Time-based metrics (generation time, age at maturity, mature life expectancy)
# are returned in the units specified by timestep.
#
# =============================================================================

library(Rpadrino)
library(tidyverse)


# Create model records ---------------------------------------------------------

# A model record stores a made IPM together with the information needed to
  # calculate its life history metrics.

# The ipm_id can be NULL for models that were not obtained from PADRINO, or you
  # can give it some descriptive ID to help you keep track of where it came from
  # downstream (like if a species has two different models from different data
  # sources, or two different treatments, etc.).

# For simple models, you may not need to specify many arguments. Most can be
  # kept as the default values.

# For general models, you may need to make some significant decisions about
  # which subkernels belong to F vs P. In general, the subkernels that generate
  # new start of life should go into F, while the rest go into P. Ex. subkernels
  # generating seeds that go into a seedbank might go into P, if there is a
  # separate subkernel describing germination, which would go into F. These can
  # be arbitrary decisions for now but we will want to examine these in detail
  # at some point.

# We might consider adding an argument to the model record that describes growth
  # form, for later visualization.


make_ipm_record <- function(ipm,
                            ipm_id             = NULL,
                            species            = NULL,
                            P_kernels          = "P",
                            F_kernels          = "F",
                            mega_mat           = NULL,
                            timestep           = 1,
                            tolerance          = 1e-8,
                            max_age            = 100,
                            repro_var          = "poisson",
                            maturity_threshold = 0.5,
                            maturity_prob      = NULL) {
  
  list(
    ipm                = ipm,
    ipm_id             = ipm_id,
    species            = species,
    P_kernels          = P_kernels,
    F_kernels          = F_kernels,
    mega_mat           = mega_mat,
    timestep           = timestep,
    tolerance          = tolerance,
    max_age            = max_age,
    repro_var          = repro_var,
    maturity_threshold = maturity_threshold,
    maturity_prob      = maturity_prob
  )
  
}

# You can also use this function to construct a model record directly from a
  # PADRINO IPM ID. Feel free to modify this function as needed to accept
  # more complex models.

make_padrino_ipm <- function(pdb,
                             ipm_id,
                             P_kernels          = "P",
                             F_kernels          = "F",
                             mega_mat           = NULL,
                             timestep           = 1,
                             tolerance          = 1e-8,
                             max_age            = 100,
                             repro_var          = "poisson",
                             maturity_threshold = 0.5,
                             maturity_prob      = NULL) {
  
  metadata_row <- match(ipm_id, pdb$Metadata$ipm_id)
  
  if (is.na(metadata_row)) {
    stop("IPM ID not found in PADRINO: ", ipm_id)
  }
  
  proto_ipm <- pdb_make_proto_ipm(
    pdb,
    ipm_id    = ipm_id,
    det_stoch = "det"
  )
  
  made_ipm <- pdb_make_ipm(proto_ipm)
  
  make_ipm_record(
    ipm                = made_ipm[[1]],
    ipm_id             = ipm_id,
    species            = pdb$Metadata$species_accepted[metadata_row],
    P_kernels          = P_kernels,
    F_kernels          = F_kernels,
    mega_mat           = mega_mat,
    timestep           = timestep,
    tolerance          = tolerance,
    max_age            = max_age,
    repro_var          = repro_var,
    maturity_threshold = maturity_threshold,
    maturity_prob      = maturity_prob
  )
  
}

# Calculate all life history metrics -------------------------------------------

# This function calculates all of the life history metrics for a single made IPM
  # object that includes all of the necessary arguments as defined in the model 
  # record

calculate_life_history_metrics <- function(record,
                                           overrides = list()) {
  
  if (!is.list(record) || is.null(record$ipm)) {
    stop("Each model record must contain an 'ipm' element.")
  }
  
  # Model-specific settings are stored in the record.
  config <- record[c( "P_kernels", "F_kernels", "mega_mat", "timestep",
                      "tolerance", "max_age", "repro_var", "maturity_threshold",
                      "maturity_prob" )]
  
  # Optional batch-level overrides
  
  # For example, to change the timestep and maturity threshold for ALL models
    # at once:
  
  # overrides = list(
  #   timestep = 1,
  #   maturity_threshold = 0.1
  # )
  
  if (length(overrides) > 0) {
    
    unknown <- setdiff(names(overrides), names(config))
    
    if (length(unknown) > 0) {
      stop(
        "Unknown override(s): ",
        paste(unknown, collapse = ", ")
      )
    }
    
    config[names(overrides)] <- overrides
    
  }
  
  ipm <- record$ipm
  
  # Calculate metrics
  
  # Time-based metrics are multiplied by timestep inside their respective
    # functions. All other metrics remain dimensionless or are on their
    # natural demographic scale
  
  tibble(
    ipm_id = record$ipm_id %||% NA_character_,
    species = record$species %||% NA_character_,
    
    R0 = R_0( ipm,
              P_kernels = config$P_kernels,
              F_kernels = config$F_kernels,
              mega_mat  = config$mega_mat ),
    
    genT = gen_time( ipm,
                     P_kernels = config$P_kernels,
                     F_kernels = config$F_kernels,
                     mega_mat  = config$mega_mat,
                     method    = "bienvenu_legendre",
                     timestep  = config$timestep ),
    
    iter = iteroparity( ipm,
                        P_kernels = config$P_kernels,
                        F_kernels = config$F_kernels,
                        mega_mat  = config$mega_mat ),
    
    s_type = surv_type( ipm,
                        P_kernels = config$P_kernels,
                        F_kernels = config$F_kernels,
                        mega_mat  = config$mega_mat,
                        tolerance = config$tolerance,
                        max_age   = config$max_age ),
    
    pgrowth = prog_growth( ipm,
                           P_kernels = config$P_kernels,
                           F_kernels = config$F_kernels,
                           mega_mat  = config$mega_mat ),
    
    rgrowth = retro_growth( ipm,
                            P_kernels = config$P_kernels,
                            F_kernels = config$F_kernels,
                            mega_mat  = config$mega_mat ),
    
    mean_repr = mean_rep( ipm,
                         P_kernels = config$P_kernels,
                         F_kernels = config$F_kernels,
                         mega_mat  = config$mega_mat ),
    
    age_matur = age_maturity( ipm,
                              P_kernels          = config$P_kernels,
                              F_kernels          = config$F_kernels,
                              mega_mat           = config$mega_mat,
                              repro_var          = config$repro_var,
                              maturity_threshold = config$maturity_threshold,
                              timestep           = config$timestep ),
    
    mature_exp = mature_expec( ipm,
                               P_kernels          = config$P_kernels,
                               F_kernels          = config$F_kernels,
                               mega_mat           = config$mega_mat,
                               repro_var          = config$repro_var,
                               maturity_threshold = config$maturity_threshold,
                               timestep           = config$timestep )
  )
  
}

# And now to calculate for multiple IPMs

calculate_life_history_metrics_batch <- function(ipm_records,
                                                 overrides = list()) {
  
  if (!is.list(ipm_records) || length(ipm_records) == 0) {
    stop("ipm_records must be a non-empty list.")
  }
  
  results <- lapply(
    ipm_records,
    calculate_life_history_metrics,
    overrides = overrides
  )
  
  bind_rows(results)
  
}


# # Testing with models from PADRINO ---------------------------------------------
# 
# # Current set of PADRINO models for initial testing
# 
# # Feel free to expand this to include more models
# 
# id_v <- pdb$Metadata$ipm_id[c(16:19,39,49,50,58:66,68)]
# 
# 
# # Convert PADRINO IDs to analysis-ready model records
# 
# padrino_records <- lapply( id_v, 
#                            function(id) make_padrino_ipm(pdb, ipm_id = id ) )
# 
# 
# # Calculate all life-history metrics for these models
# 
# life_history_df <- calculate_life_history_metrics_batch( padrino_records,
#                                                          overrides = list(
#                                                            maturity_threshold = 0.1
#                                                          ))
# 
# life_history_df
# 
# # For the current round of testing, we are overriding the maturity threshold of
#   # 0.5 and setting it at 0.1. This means "maturity" is reached at the first 
#   # size class where the probability of reproduction is 10%. We should examine
#   # which models failed to meet the threshold of 0.5, and possibly reexamine
#   # this calculation overall to just return NA when no states meet the threshold
#   # rather than breaking the loop.
# 
# 
# 
# # To expand the testing to include models not in PADRINO, use make_ipm_record()
#   # with all of the relevant arguments defined for your IPM, then append to the
#   # padrino_records list and run the list through calculate_life_history_metrics
# 
# # my_model <- make_ipm_record(
# #   ipm       = my_made_ipm,
# #   ipm_id    = "my_id",
# #   species   = "my_species",
# #   P_kernels = "P",
# #   F_kernels = "F",
# #   timestep  = 1 )
# 
# # or for a general model:
# 
# # my_model <- make_ipm_record(
# #   ipm       = my_made_ipm,
# #   ipm_id    = "my_id",
# #   species   = "my_species",
# #   P_kernels = c("P", "survival"),
# #   F_kernels = c("repr", "seed_recruitment"),
# #   mega_mat  = "c(P + repr, seed_recruitment, SB)" )
# 
# # and then append:
#  
# # all_records <- c( padrino_records, list(my_model) )
# 
# # life_history_df <- calculate_life_history_metrics_batch( all_records )
# 
# 
# # Generation time comparison ---------------------------------------------------
# 
# # Compare the original generation-time calculation with the revised one:
#   # 1. Bienvenu-Legendre method
#   # 2. log(R0) / log(lambda)
# 
# generation_time_compare <- function( record ){
#   
#   ipm <- record$ipm
#   
#   Pm <- get_P_F_kernels( ipm,
#                          P_kernels = record$P_kernels,
#                          F_kernels = record$F_kernels,
#                          mega_mat  = record$mega_mat )$P
#   
#   Fm <- get_P_F_kernels( ipm,
#                          P_kernels = record$P_kernels,
#                          F_kernels = record$F_kernels,
#                          mega_mat  = record$mega_mat )$F
#   
#   K <- Pm + Fm
#   
#   lam <- dominant_eigen(K)
#   R0  <- R_0( ipm,
#               P_kernels = record$P_kernels,
#               F_kernels = record$F_kernels,
#               mega_mat  = record$mega_mat )
#   
#   genT_R0 <- if (
#     is.finite(lam) &&
#     is.finite(R0) &&
#     lam > 0 &&
#     R0 > 0 &&
#     !isTRUE(all.equal(lam, 1))
#   ) {
#     log(R0) / log(lam) * record$timestep
#   } else {
#     NA_real_
#   }
#   
#   tibble(
#     ipm_id                 = record$ipm_id %||% NA_character_,
#     species                = record$species %||% NA_character_,
#     genT_bienvenu_legendre = gen_time(
#       ipm,
#       P_kernels = record$P_kernels,
#       F_kernels = record$F_kernels,
#       mega_mat  = record$mega_mat,
#       method    = "bienvenu_legendre",
#       timestep  = record$timestep
#     ),
#     genT_R0 = genT_R0
#   )
#   
# }
# 
# 
# generation_time_df <- bind_rows(
#   lapply(padrino_records, generation_time_compare)
# )
# 
# # Plot the two generation-time estimates against one another
# 
# p_generation_time <- ggplot( generation_time_df,
#                              aes( x    = genT_bienvenu_legendre,
#                                   y    = genT_R0,
#                                   fill = species ) ) +
#   geom_abline( intercept = 0, slope = 1, linewidth = 0.5, color = "grey60" ) +
#   geom_point( shape = 21, size = 2.5, alpha = 0.85 ) +
#   labs( x    = "Generation time: Bienvenu & Legendre",
#         y    = "Generation time: log(R0) / log(lambda)",
#         fill = "Species" ) +
#   theme_bw()
# 
# p_generation_time
# 
# # So far, it looks like these two estimates are pretty different. We should
#   # consider which to use in the end.
# 
# 
# # Plotting iteroparity vs generation time --------------------------------------
# 
# # This preliminary plot rescales iteroparity and generation time arbitrarily
#   # based on the models currently included in the analysis. I assume that most
#   # of these models are fairly short-lived perennials, so they would in theory
#   # occupy the upper left quadrant of the final plot. As we expand the set of 
#   # models included in the exploratory analysis, the axes will adjust
#   # accordingly. But we should actually aim to select biologically meaningful
#   # thresholds for these values, rather than just rescaling around the sample
#   # medians. 
# 
# plot_data <- life_history_df %>%
#   mutate( genT_scaled = scales::rescale(genT, to = c(-1, 1), na.rm = TRUE),
#           iter_scaled = scales::rescale(iter, to = c(-1, 1), na.rm = TRUE) )
# 
# 
# p_iter_genT <- ggplot( plot_data, aes( x    = genT_scaled,
#                                        y    = iter_scaled,
#                                        fill = species ) ) +
#   geom_point( shape = 21, size = 2.5, alpha = 0.85, na.rm = TRUE ) +
#   geom_vline( xintercept = 0, linewidth = 0.4, colour = "grey60" ) +
#   geom_hline( yintercept = 0, linewidth = 0.4, colour = "grey60" ) +
#   coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
#   labs( x    = "Generation time",
#         y    = "Iteroparity",
#         fill = "Species" ) +
#   theme_bw()
# 
# p_iter_genT
# 
# 
# # This is where a growth form argument might make more sense for visualization
#   # than coloring the points by species. As we expand the diversity of models
#   # in the dataset, and the axes adjust accordingly, it would be nice to see how
#   # different growth forms occupy this space, and this might help us determine
#   # the biologically-relevant thresholds for semelparity vs iteroparity in the
#   # iteroparity metric, and short- vs. long-lived in the generation time metric.
# 
# 
# # PCA of life history metrics --------------------------------------------------
# 
# metric_cols <- c(
#   "R0",
#   "genT",
#   "iter",
#   "s_type",
#   "pgrowth",
#   "rgrowth",
#   "mean_repr",
#   "age_matur",
#   "mature_exp"
# )
# 
# metric_labels <- c(
#   R0         = "Net reproductive rate (R0)",
#   genT       = "Generation time (T)",
#   iter       = "Iteroparity (S)",
#   s_type     = "Survivorship curve type (H)",
#   pgrowth    = "Progressive growth (γ)",
#   rgrowth    = "Retrogressive growth (ρ)",
#   mean_repr  = "Mean sexual reproduction (Φ)",
#   age_matur  = "Age at maturity (Lα)",
#   mature_exp = "Mature life expectancy (Lω)"
# )
# 
# metric_process <- c(
#   R0         = "Reproduction",
#   genT       = "Turnover",
#   iter       = "Reproduction",
#   s_type     = "Longevity",
#   pgrowth    = "Growth",
#   rgrowth    = "Growth",
#   mean_repr  = "Reproduction",
#   age_matur  = "Longevity",
#   mature_exp = "Longevity"
# )
# 
# 
# # Keep only IPMs with complete, finite metric values
# 
# pca_data <- life_history_df %>%
#   select(ipm_id, species, all_of(metric_cols)) %>%
#   filter( if_all( all_of(metric_cols), ~ is.finite(.x) ) )
# 
# 
# excluded_ipms <- anti_join(
#   life_history_df %>% select(ipm_id, species),
#   pca_data %>% select(ipm_id, species),
#   by = c("ipm_id", "species")
# )
# 
# 
# if (nrow(excluded_ipms) > 0) {
#   
#   message(
#     "Excluded ",
#     nrow(excluded_ipms),
#     " IPM(s) from the PCA because one or more metrics were missing ",
#     "or non-finite."
#   )
#   
#   print(excluded_ipms)
#   
# }
# 
# 
# if (nrow(pca_data) < 3) {
#   stop("Too few complete IPMs to calculate a useful two-axis PCA.")
# }
# 
# 
# # All nine metrics should be nonnegative for this transformation.
# metric_matrix <- as.matrix(pca_data[, metric_cols])
# 
# if (any(metric_matrix < 0)) {
#   stop(
#     "At least one PCA metric is negative. ",
#     "Review the metric definitions before applying log1p()."
#   )
# }
# 
# # Transform and scale metrics
# 
# # All nine metrics are expected to be nonnegative.
# # log1p() accommodates zero values while reducing right skew.
# 
# # PCA itself is then performed on centered and unit-variance-scaled metrics.
# 
# metric_matrix <- as.matrix(pca_data[, metric_cols])
# 
# 
# if (any(metric_matrix < 0, na.rm = TRUE)) {
#   
#   negative_metrics <- pca_data %>%
#     select(ipm_id, species, all_of(metric_cols)) %>%
#     pivot_longer( cols      = all_of(metric_cols),
#                   names_to  = "metric",
#                   values_to = "value" ) %>%
#     filter( value < 0 )
#   
#   print(negative_metrics)
#   
#   stop(
#     "At least one life-history metric is negative. ",
#     "Inspect the output above before running the PCA."
#   )
#   
# }
# 
# 
# metric_matrix_transformed <- log1p(metric_matrix)
# 
# 
# # Remove metrics with zero variance after transformation
# 
# metric_sd <- apply(metric_matrix_transformed, 2, sd)
# 
# keep_metrics <- is.finite(metric_sd) & metric_sd > 0
# 
# if (any(!keep_metrics)) {
#   message(
#     "Dropping zero-variance metric(s) from the PCA: ",
#     paste(names(keep_metrics)[!keep_metrics], collapse = ", ")
#   )
# }
# 
# metric_matrix_transformed <- metric_matrix_transformed[
#   ,
#   keep_metrics,
#   drop = FALSE
# ]
# 
# if (ncol(metric_matrix_transformed) < 2) {
#   stop("Fewer than two metrics remain for PCA after variance checks.")
# }
# 
# pc <- prcomp( metric_matrix_transformed, center = TRUE, scale. = TRUE )
# 
# 
# # PCA scores
# 
# scores <- as.data.frame( pc$x[, 1:2, drop = FALSE] ) %>%
#   mutate( ipm_id  = pca_data$ipm_id,
#           species = pca_data$species )
# 
# names(scores)[1:2] <- c("PC1", "PC2")
# 
# 
# # PCA loadings
# 
# loadings <- as.data.frame( pc$rotation[, 1:2, drop = FALSE] ) %>%
#   rownames_to_column("loading") %>%
#   mutate( metric  = unname(metric_labels[loading]),
#           process = unname(metric_process[loading]) )
# 
# # Rescale loading arrows for visualiztion
# 
# score_range   <- max( abs(c(scores$PC1, scores$PC2)), na.rm = TRUE )
# 
# loading_range <- max( abs(c(loadings$PC1, loadings$PC2)), na.rm = TRUE )
# 
# arrow_scale <- 0.75 * score_range / loading_range
# 
# 
# loadings <- loadings %>% mutate( x    = 0,
#                                  y    = 0,
#                                  xend = PC1 * arrow_scale,
#                                  yend = PC2 * arrow_scale )
# 
# 
# # Variance explained by each PC
# 
# pca_variance <- summary(pc)$importance[2, ]
# 
# ex_var1 <- round(100 * pca_variance[1], 2)
# ex_var2 <- round(100 * pca_variance[2], 2)
# 
# 
# # Plot PCA scores and metric loadings
# 
# p_life_history_pca <- ggplot() +
#   geom_hline( yintercept = 0, linewidth = 0.4, color = "grey75" ) +
#   geom_vline( xintercept = 0, linewidth = 0.4, color = "grey75" ) +
#   geom_point( data = scores,
#               aes(x = PC1, y = PC2, fill = species),
#               shape = 21, size = 2.5, alpha = 0.85 ) +
#   geom_segment( data = loadings,
#                 aes( x = x, y = y,
#                      xend = xend,
#                      yend = yend,
#                      color = process ),
#                 linewidth = 0.8, 
#                 arrow = grid::arrow( length = grid::unit(0.18, "cm" ) ) ) +
#   geom_text( data = loadings,
#              aes( x = xend,
#                   y = yend,
#                   label = metric,
#                   color = process ),
#              size = 3,
#              vjust = -0.6,
#              show.legend = FALSE ) +
#   scale_color_manual( values = c( "darkblue", "gray", "darkred", "black" ) ) +
#   labs( x = paste0("PC1 (", ex_var1, "%)"),
#         y = paste0("PC2 (", ex_var2, "%)"),
#         fill = "Species",
#         color = "Life history process" ) +
#   coord_equal() +
#   theme_bw()
# 
# p_life_history_pca
# 
# 
# # Ideally these parameters should map relatively similar to Fig 2A in
#   # Salguero-Gomez et al. 2015
# 
# # I've colored the arrows following the same scheme as in Fig 2A to facilitate
#   # comparisons
# 
# # If we do add an argument describing growth form, we could consider coloring
#   # the points according to growth form instead of species as in Fig 2A, at some  
#   # point
# 
# # Note that we currently set the maturity threshold to 0.1, rather than 0.5 
#   # (which might be a more realistic threshold in the end), which results in
#   # a lower age at maturity across all models, and earlier maturity results in
#   # more time spent in the mature portion of life, so the mature life expectancy
#   # is generally overestimated. This is not necessarily evidence that a species
#   # matures earlier and should not be overinterpreted as such. We should try 
#   # using different maturity thresholds and examining why some models never
#   # reach the 0.5 threshold.

