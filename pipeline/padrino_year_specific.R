# Populating padrino - IMP year specific - pipeline

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.11.11


# Setting the stage ------------------------------------------------------------
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)


# Packages ---------------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R' )
load_packages(tidyverse, ipmr, readxl, writexl, remotes)
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
v_suffix <- read.csv(file.path(dir_data, 'v_suffix.csv'))

# Ipm mean and plant tracker if they not already exists
if (
  !file.exists(file.path(
    dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_data_df.csv')))) {
  source(paste0(
    dir_R, '/', v_script_prefix, '_', v_sp_abb, '_ipm_year_specific.R'))
}

all_pars      <- read.csv(file.path(dir_data, 'all_pars.csv'))
pars_var_wide <- read.csv(file.path(dir_data, '2.pars_var.csv'))
lam_mean_ipmr <- read.csv(file.path(dir_data, 'lambdas_yr_vec.csv'))


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

pdb$StateVariables[1,] <- c(v_ipm_id, 'size', FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(
  v_ipm_id, 'size', '', all_pars$L, all_pars$U, 'P_yr; F_yr', '')

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(
  v_ipm_id, 'size', '', all_pars$mat_siz, 'midpoint', 'P_yr; F_yr')

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(
  v_ipm_id, 'n_size', all_pars$mat_siz, '')
pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(
  v_ipm_id, 'P_yr', 'P_yr = mu_s_yr * g_yr * d_size', 'CC', 'size', 'size')

pdb$IpmKernels[2,] <- c(
  v_ipm_id, 'F_yr', 'F_yr = fy_yr * d_size', 'CC', 'size', 'size')

# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(
  v_ipm_id, 'Survival', 
  if (all_pars$mod_su_index == 1) {
    v_mu_s_yr_expr <- paste(
      "mu_s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1)))")
    } 
  else if (all_pars$mod_su_index == 2) {
    v_mu_s_yr_expr <- paste(
      "mu_s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2)))")
    } 
  else if (all_pars$mod_su_index == 3) {
    v_mu_s_yr_expr <- paste(
      "mu_s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2 + surv_b3_yr * size_1^3)))")
    },
  'Evaluated', 'P_yr')

pdb$VitalRateExpr[2,] <- c(
  v_ipm_id, 'Growth',
  if (all_pars$mod_gr_index == 1) {
    v_mu_g_yr_expr <- paste(
      "mu_g_yr = grow_b0_yr + grow_b1_yr * size_1")
    } 
  else if (all_pars$mod_gr_index == 2) {
    v_mu_g_yr_expr <- paste(
      "mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2")
    } 
  else if (all_pars$mod_gr_index == 3) {
    v_mu_g_yr_expr <- paste(
      "mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2 + grow_b3_yr * size_1^3")
    },
  'Evaluated', 'P_yr')

pdb$VitalRateExpr[3,] <- c(
  v_ipm_id, 'Growth', 'g_yr = Norm(mu_g_yr, sd_g)', 'Substituted', 'P_yr')

pdb$VitalRateExpr[4,] <- c(
  v_ipm_id, 'Growth', 'sd_g = sqrt(a * exp(b * size_1))', 'Evaluated', 'P_yr')

pdb$VitalRateExpr[5,] <- c(
  v_ipm_id, 'Fecundity', 'fy_yr = fecu_b0_yr * r_d', 'Evaluated', 'F_yr')

pdb$VitalRateExpr[6,] <- c(
  v_ipm_id, 'Fecundity', 'r_d = Norm(recr_sz, recr_sd)', 'Substituted', 'F_yr')


# Parameter Values
for(i in 1:(length(pars_var_wide))) {
  pdb$ParameterValues[i,1] <- v_ipm_id
  pdb$ParameterValues[i,3] <- 'size'
  pdb$ParameterValues[i,4] <- names(pars_var_wide)[i]
  pdb$ParameterValues[i,5] <- as.numeric(pars_var_wide[i])
  
  if(grepl('surv', names(pars_var_wide)[i])){
    pdb$ParameterValues[i,2] <- 'Survival'
  } else {
    if(grepl('grow', names(pars_var_wide)[i])){
      pdb$ParameterValues[i,2] <- 'Growth'
    } else {pdb$ParameterValues[i,2] <- 'Fecundity'}
  }
}

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Growth', 'size', 'a', all_pars$a)
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Growth', 'size', 'b', all_pars$b)                               
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Fecundity', 'size', 'recr_sz', all_pars$recr_sz)
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(
  v_ipm_id, 'Fecundity', 'size', 'recr_sd', all_pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(
  pdb$ParameterValues$parameter_value)


# Environmental variables
pdb$ParSetIndices[1,] <- c(
  v_ipm_id, 'year', 'yr', 'lam_mean_ipmr$years', 'P_yr; F_yr', '')

# Test targets
pdb$TestTargets[1:nrow(lam_mean_ipmr),1] <- v_ipm_id
pdb$TestTargets[1:nrow(lam_mean_ipmr),2] <- 1:nrow(lam_mean_ipmr)
pdb$TestTargets[1:nrow(lam_mean_ipmr),3] <- as.numeric(lam_mean_ipmr$value)
pdb$TestTargets[1:nrow(lam_mean_ipmr),4] <- 3

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision    <- as.numeric(pdb$TestTargets$precision)


write_xlsx(
  pdb, file.path(
    dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_pdb_yr.xlsx')))

pdb_test       <- read_pdb(
  file.path(dir_data, paste0(
    v_script_prefix, '_', v_sp_abb, '_pdb_yr.xlsx')))

pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = 'det')
print(pdb_test_proto[[v_ipm_id]])
bg_ipm_pdb     <- make_ipm(pdb_test_proto[[v_ipm_id]])
