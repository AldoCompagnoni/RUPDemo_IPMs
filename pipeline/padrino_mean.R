# Populating padrino - IMP mean - pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.01.09


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
sp_abb  <- tolower(gsub(' ', '', paste(
  substr(unlist(strsplit(species, ' ')), 1, 2), collapse = '')))
# Directory 
pub_dir    <- file.path(paste0(author_year, '_', region_abb))
R_dir      <- file.path(pub_dir, 'R',       sp_abb)
data_dir   <- file.path(pub_dir, 'data',    sp_abb)
result_dir <- file.path(pub_dir, 'results', sp_abb)

# Prefix for all the files
script_prefix <- str_c(str_extract(author_year, '^[^_]+'), 
                       str_sub(str_extract(author_year, '_\\d+$'), -2, -1))
if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0("^", author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  script_prefix <- paste0(script_prefix, region_abb)
}

v_suffix <- read.csv(file.path(data_dir, 'v_suffix.csv'))

# Ipm mean and plant tracker if they not already exists
if (!file.exists(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_data_df.csv'))) {
  source(paste0(R_dir, '/', script_prefix, '_', sp_abb, '_ipm_mean.R'))
}

pars          <- read.csv(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_pars.csv'))
lam_mean_ipmr <- read.csv(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_lambda.csv'))

v_suffix


# Populating the PADRINO database template -------------------------------------
# Store the PADRINO excel template
YOUR_PATH   <- '.'
sheet_names <- excel_sheets(paste(YOUR_PATH, '/pdb_template.xlsx', sep = ''))
pdb         <- lapply(sheet_names, function(x) {
  as.data.frame(read_excel(paste(YOUR_PATH, '/pdb_template.xlsx', sep = ''), sheet = x))})
names(pdb)  <- sheet_names


# Information ------------------------------------------------------------------
year_dur <- start_year - end_year + 1

pdb$Metadata[1,] <- c(
  ipm_id, 
  
  # Taxonomic information
  gsub(" ", "_", species), species_accepted, tax_genus,
  tax_family, tax_order, tax_class, tax_phylum,
  kingdom, organism_type, dicot_monocot, angio_gymno, 
  
  # Publication information
  authors,
  journal, pub_year, doi, corresponding_author, 
  email_year, remark, 
  apa_citation,
  demog_appendix_link,
  
  # Data collection information
  year_dur, start_year, start_month, end_year, end_month, periodicity, 
  population_name, number_populations, 
  lat, lon, altitude, country,
  continent, ecoregion,
  
  # Model information
  'A', TRUE, 'truncated_distributions', NA, NA, FALSE,
  FALSE, FALSE, FALSE, '', '', ''
)

pdb$Metadata$eviction_used <- as.logical(pdb$Metadata$eviction_used)
pdb$Metadata$duration <- as.numeric(pdb$Metadata$duration)
pdb$Metadata$periodicity <- as.numeric(pdb$Metadata$periodicity)


pdb$StateVariables[1,]      <- c(ipm_id, 'size', FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(ipm_id, 
                               'size', 
                               '', 
                               pars$L, 
                               pars$U, 
                               'P; F', 
                               '')

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(ipm_id,
                              'size',
                              '',
                              pars$mat_siz,
                              'midpoint',
                              'P; F')

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(ipm_id,
                          'n_size',
                          pars$mat_siz,
                          '')
pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(ipm_id, 
                        'P', 
                        'P = s * g * d_size', 
                        'CC', 
                        'size', 
                        'size')

pdb$IpmKernels[2,] <- c(ipm_id, 
                        'F', 
                        'F= fy * d_size', 
                        'CC', 
                        'size', 
                        'size')


# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(ipm_id,
                           'Survival',
                           if (pars$mod_su_index == 1) {
                             s_expr <- paste('s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1)))')
                           } else if (pars$mod_su_index == 2) {
                             s_expr <- paste("s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1 + surv_b2 * size_1^2)))")
                           } else if (pars$mod_su_index == 3) {
                             s_expr <- paste("s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1 + surv_b2 * size_1^2 + surv_b3 * size_1^3)))")
                           },
                           'Evaluated',
                           'P')

pdb$VitalRateExpr[2,] <- c(ipm_id,
                           'Growth',
                           if (pars$mod_gr_index  == 1) {
                             mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1")
                           } else if (pars$mod_gr_index  == 2) {
                             mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1 + grow_b2 * size_1^2")
                           } else if (pars$mod_gr_index  == 3) {
                             mu_g_expr <- paste("mu_g = grow_b0 + grow_b1 * size_1 + grow_b2 * size_1^2 + grow_b3 * size_1^3")
                           },
                           'Evaluated',
                           'P')

pdb$VitalRateExpr[3,] <- c(ipm_id,
                           'Growth',
                           'g = Norm(mu_g, sd_g)',
                           'Substituted',
                           'P')

pdb$VitalRateExpr[4,] <- c(ipm_id,
                           'Growth',
                           'sd_g = sqrt(a * exp(b * size_1))',
                           'Evaluated',
                           'P')

pdb$VitalRateExpr[5,] <- c(ipm_id,
                           'Fecundity',
                           'fy = fecu_b0 * r_d',
                           'Evaluated',
                           'F')

pdb$VitalRateExpr[6,] <- c(ipm_id,
                           'Fecundity',
                           'r_d = Norm(recr_sz, recr_sd)',
                           'Substituted',
                           'F')


# Parameter Values
# Populating grow_ and surv_ dynamically
for (i in 1:nrow(pars)) {
  
  # Extract the mod_su_index and mod_gr_index  for this row
  su_index <- pars$mod_su_index[i]
  gr_index <- pars$mod_gr_index [i]
  
  # Get the parameters for surv and grow
  surv_param <- paste0("surv_b", su_index)  # E.g., "surv_b1"
  grow_param <- paste0("grow_b", gr_index)  # E.g., "grow_b2"
  
  # Populate for survival parameters
  if (su_index == 0) {
    # For su_index == 0, we add only surv_b0
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival 
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
  } else if (su_index == 1) {
    # For su_index == 1, we add surv_b0 and surv_b1
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
  } else if (su_index == 2) {
    # For su_index == 2, we add surv_b0, surv_b1, and surv_b2
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b2",         # Survival parameter surv_b2
      pars$surv_b2[i]    # Value of the survival parameter surv_b2
    )
  } else if (su_index == 3) {
    # For su_index == 3, we add surv_b0, surv_b1, surv_b2, and surv_b3
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b0",         # Survival parameter surv_b0
      pars$surv_b0[i]    # Value of the survival parameter surv_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b1",         # Survival parameter surv_b1
      pars$surv_b1[i]    # Value of the survival parameter surv_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b2",         # Survival parameter surv_b2
      pars$surv_b2[i]    # Value of the survival parameter surv_b2
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Survival",        # Category for survival
      "size",            # Type for size
      "surv_b3",         # Survival parameter surv_b3
      pars$surv_b3[i]    # Value of the survival parameter surv_b3
    )
  }
  
  # Now, handle the growth parameters dynamically based on gr_index
  if (gr_index == 0) {
    # For gr_index == 0, we add only grow_b0
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
  } else if (gr_index == 1) {
    # For gr_index == 1, we add grow_b0 and grow_b1
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
  } else if (gr_index == 2) {
    # For gr_index == 2, we add grow_b0, grow_b1, and grow_b2
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b2",         # Growth parameter grow_b2
      pars$grow_b2[i]    # Value of the growth parameter grow_b2
    )
  } else if (gr_index == 3) {
    # For gr_index == 3, we add grow_b0, grow_b1, grow_b2, and grow_b3
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b0",         # Growth parameter grow_b0
      pars$grow_b0[i]    # Value of the growth parameter grow_b0
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b1",         # Growth parameter grow_b1
      pars$grow_b1[i]    # Value of the growth parameter grow_b1
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b2",         # Growth parameter grow_b2
      pars$grow_b2[i]    # Value of the growth parameter grow_b2
    )
    
    pdb$ParameterValues[nrow(pdb$ParameterValues) + 1, ] <- c(
      ipm_id,            # Assuming you have an `ipm_id` variable
      "Growth",          # Category for growth
      "size",            # Type for size
      "grow_b3",         # Growth parameter grow_b3
      pars$grow_b3[i]    # Value of the growth parameter grow_b3
    )
  }
}


pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Growth',
    'size',
    'a',
    pars$a)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Growth',
    'size',
    'b',
    pars$b)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Fecundity',
    'size',
    'fecu_b0',
    pars$fecu_b0)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Fecundity',
    'size',
    'recr_sz',
    pars$recr_sz)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Fecundity',
    'size',
    'recr_sd',
    pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(
  pdb$ParameterValues$parameter_value)


# Test targets
pdb$TestTargets[1,] <- c(ipm_id,
                         'lambda',
                         lam_mean_ipmr,
                         3)

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision <- as.numeric(pdb$TestTargets$precision)


write_xlsx(pdb, 
           paste0(data_dir, '/', sp_abb, '_pdb.xlsx'))

pdb_test       <- read_pdb(
  paste0(data_dir, '/', sp_abb, '_pdb.xlsx'))

pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = 'det')
print(pdb_test_proto[[ipm_id]])
bg_ipm_pdb     <- make_ipm(pdb_test_proto[[ipm_id]])
