# Populating padrino - IMP mean - pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.05


# Setting the stage ------------------------------------------------------------
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)


# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R' )
load_packages(tidyverse, ipmr, readxl, writexl, remotes, janitor)
# remotes::install_github("padrinoDB/pdbDigitUtils",force = TRUE)
library(pdbDigitUtils)


# Data -------------------------------------------------------------------------
# Species abbriviation
v_sp_abb  <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))
# Directory 
dir_pub    <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

# Prefix for all the files
v_script_prefix <- str_c(
  str_extract(v_author_year, '^[^_]+'),
  str_sub(str_extract(v_author_year, '_\\d+$'), -2, -1))
# Define prefix for two of the same author and year
if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0("^", v_author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0("^", v_author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

v_suffix <- read.csv(file.path(dir_data, 'v_suffix.csv'))

# Ipm mean and plant tracker if they not already exists
if (!file.exists(
  file.path(dir_data, paste0(
    v_script_prefix, '_', v_sp_abb, '_data_df.csv')))) {
  source(file.path(dir_data,
    paste0(v_script_prefix, '_', v_sp_abb, '_ipm_mean.R')))
}

pars          <- read.csv(
  file.path(dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_pars.csv')))
lam_mean_ipmr <- read.csv(
  file.path(dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_lambda.csv')))

v_suffix


# Populating the PADRINO database template -------------------------------------
# Store the PADRINO excel template
dir_pdb_temp <- c('.')
v_sheet_names  <- excel_sheets(file.path(dir_pdb_temp, 'pdb_template.xlsx'))
pdb          <- lapply(v_sheet_names, function(x) {
  as.data.frame(
    read_excel(file.path(dir_pdb_temp, 'pdb_template.xlsx'), sheet = x))})
names(pdb)   <- v_sheet_names


# Information ------------------------------------------------------------------
v_year_dur <- v_end_year - v_start_year + 1

pdb$Metadata[1,] <- c(
  v_ipm_id, 
  
  # Taxonomic information
  gsub(" ", "_", v_species), v_species_accepted, v_tax_genus,
  v_tax_family, v_tax_order, v_tax_class, v_tax_phylum,
  v_kingdom, v_organism_type, v_dicot_monocot, v_angio_gymno, 
  
  # Publication information
  v_authors,
  v_journal, v_pub_year, v_doi, v_corresponding_author, 
  v_email_year, v_remark, 
  v_apa_citation,
  v_demog_appendix_link,
  
  # Data collection information
  v_year_dur, v_start_year, v_start_month, v_end_year, v_end_month, v_periodicity, 
  v_population_name, v_number_populations, 
  v_lat, v_lon, v_altitude, v_country,
  v_continent, v_ecoregion,
  
  # Model information
  'A', TRUE, 'truncated_distributions', NA, NA, FALSE,
  FALSE, FALSE, FALSE, '', '', ''
)

pdb$Metadata$eviction_used <- as.logical(pdb$Metadata$eviction_used)
pdb$Metadata$duration <- as.numeric(pdb$Metadata$duration)
pdb$Metadata$periodicity <- as.numeric(pdb$Metadata$periodicity)


pdb$StateVariables[1,]      <- c(v_ipm_id, 'size', FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(
  v_ipm_id, 'size', '', pars$L, pars$U, 'P; F', '')

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(
  v_ipm_id, 'size', '', pars$mat_siz, 'midpoint', 'P; F')

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(
  v_ipm_id, 'n_size', pars$mat_siz, '')

pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(
  v_ipm_id, 'P', 'P = s * g * d_size', 'CC', 'size', 'size')

pdb$IpmKernels[2,] <- c(
  v_ipm_id, 'F', 'F= fy * d_size', 'CC', 'size', 'size')


# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(
  v_ipm_id, 'Survival', 
  if (pars$mod_su_index == 1) {
    v_s_expr <- paste('s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1)))')
    } 
  else if (pars$mod_su_index == 2) {
    v_s_expr <- paste("s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1 + surv_b2 * size_1^2)))")
    } 
  else if (pars$mod_su_index == 3) {
    v_s_expr <- paste("s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1 + surv_b2 * size_1^2 + surv_b3 * size_1^3)))")
    },
  'Evaluated', 'P')

pdb$VitalRateExpr[2,] <- c(
  v_ipm_id, 'Growth',
  if (pars$mod_gr_index  == 1) {
    v_mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1")
    } 
  else if 
  (pars$mod_gr_index  == 2) {
    v_mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1 + grow_b2 * size_1^2")
    } 
  else if (pars$mod_gr_index  == 3) {
    v_mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1 + grow_b2 * size_1^2 + grow_b3 * size_1^3")
    },
  'Evaluated', 'P')

pdb$VitalRateExpr[3,] <- c(
  v_ipm_id, 'Growth', 'g = Norm(mu_g, sd_g)', 'Substituted', 'P')

pdb$VitalRateExpr[4,] <- c(
  v_ipm_id, 'Growth', 'sd_g = sqrt(a * exp(b * size_1))', 'Evaluated', 'P')

pdb$VitalRateExpr[5,] <- c(
  v_ipm_id, 'Fecundity', 'fy = fecu_b0 * r_d', 'Evaluated', 'F')

pdb$VitalRateExpr[6,] <- c(
  v_ipm_id, 'Fecundity', 'r_d = Norm(recr_sz, recr_sd)', 'Substituted', 'F')


# Parameter Values
# Populating grow_ and surv_ dynamically
for (i in 1:nrow(pars)) {
  
  # Extract the mod_su_index and mod_gr_index  for this row
  v_mod_su_index <- pars$mod_su_index[i]
  v_mod_gr_index <- pars$mod_gr_index [i]
  
  # Get the parameters for surv and grow
  v_mod_surv_param <- paste0("surv_b", v_mod_su_index)  # E.g., "surv_b1"
  v_mod_grow_param <- paste0("grow_b", v_mod_gr_index)  # E.g., "grow_b2"
  
  # Populate for survival parameters
  if (v_mod_su_index == 0) {
    # For v_mod_su_index == 0, we add only surv_b0
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival 
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
  } else if (v_mod_su_index == 1) {
    # For v_mod_su_index == 1, we add surv_b0 and surv_b1
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
  } else if (v_mod_su_index == 2) {
    # For v_mod_su_index == 2, we add surv_b0, surv_b1, and surv_b2
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b2",         # Survival parameter surv_b2
      pars$surv_b2[i]    # Value of the survival parameter surv_b2
    )
  } else if (v_mod_su_index == 3) {
    # For v_mod_su_index == 3, we add surv_b0, surv_b1, surv_b2, and surv_b3
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b2",         # Survival parameter surv_b2
      pars$surv_b2[i]    # Value of the survival parameter surv_b2
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b3",         # Survival parameter surv_b3
      pars$surv_b3[i]    # Value of the survival parameter surv_b3
    )
  }
  
  # Now, handle the growth parameters dynamically based on v_mod_gr_index
  if (v_mod_gr_index == 0) {
    # For v_mod_gr_index == 0, we add only grow_b0
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
  } else if (v_mod_gr_index == 1) {
    # For v_mod_gr_index == 1, we add grow_b0 and grow_b1
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
  } else if (v_mod_gr_index == 2) {
    # For v_mod_gr_index == 2, we add grow_b0, grow_b1, and grow_b2
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b2",         # Growth parameter grow_b2
      pars$grow_b2[i]    # Value of the growth parameter grow_b2
    )
  } else if (v_mod_gr_index == 3) {
    # For v_mod_gr_index == 3, we add grow_b0, grow_b1, grow_b2, and grow_b3
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b2",         # Growth parameter grow_b2
      pars$grow_b2[i]    # Value of the growth parameter grow_b2
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      v_ipm_id,            # Assuming you have an `v_ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b3",         # Growth parameter grow_b3
      pars$grow_b3[i]    # Value of the growth parameter grow_b3
    )
  }
}


pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Growth', 'size', 'a', pars$a)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Growth', 'size', 'b', pars$b)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Fecundity', 'size', 'fecu_b0', pars$fecu_b0)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Fecundity', 'size', 'recr_sz', pars$recr_sz)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Fecundity', 'size', 'recr_sd', pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(
  pdb$ParameterValues$parameter_value)


# Test targets
pdb$TestTargets[1,] <- c(
  v_ipm_id, 'lambda', lam_mean_ipmr, 3)

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision <- as.numeric(pdb$TestTargets$precision)


write_xlsx(
  pdb, paste0(dir_data, '/', v_sp_abb, '_pdb.xlsx'))

pdb_test       <- read_pdb(
  paste0(dir_data, '/', v_sp_abb, '_pdb.xlsx'))

pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = 'det')
print(pdb_test_proto[[v_ipm_id]])
bg_ipm_pdb     <- make_ipm(pdb_test_proto[[v_ipm_id]])
