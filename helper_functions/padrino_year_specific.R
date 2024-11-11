# Populating padrino - IMP year specific - wrapper

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.11.11


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
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
sp_abb  <- 
  tolower(gsub(' ', '', paste(substr(unlist(strsplit(species, ' ')), 1, 2), 
                              collapse = '')))
# Directory 
pub_dir    <- file.path(paste0(author_year, '_', region_abb))
R_dir      <- file.path(pub_dir, 'R',      sp_abb)
data_dir   <- file.path(pub_dir, 'data',   sp_abb)
result_dir <- file.path(pub_dir, 'result', sp_abb)

# Prefix for all the files
script_prefix <- str_c(str_extract(author_year, '^[^_]+'), 
                       str_sub(str_extract(author_year, '_\\d+$'), -2, -1))

# Ipm mean and plant tracker if they not already exists
if (!file.exists(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_data_df.csv'))) {
  source(paste0(R_dir, '/', script_prefix, '_', sp_abb, '_ipm_year_specific.R'))
}

all_pars      <- read.csv(paste0(data_dir, '/all_pars.csv'))
pars_var_wide <- read.csv(paste0(data_dir, '/2.pars_var.csv'))
lam_mean_ipmr <- read.csv(paste0(data_dir, '/lambdas_yr_vec.csv'))


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

pdb$StateVariables[1,] <- c(ipm_id, 'size', FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(ipm_id, 
                               'size', 
                               '', 
                               all_pars$L, 
                               all_pars$U, 
                               'P_yr; F_yr', 
                               '')

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(ipm_id,
                              'size',
                              '',
                              all_pars$mat_siz,
                              'midpoint',
                              'P_yr; F_yr')

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(ipm_id,
                          'n_size',
                          all_pars$mat_siz,
                          '' )
pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(ipm_id, 
                        'P_yr', 
                        'P_yr = s_yr * g_yr * d_size', 
                        'CC', 
                        'size', 
                        'size')

pdb$IpmKernels[2,] <- c(ipm_id, 
                        'F_yr', 
                        'F_yr = fy_yr * d_size', 
                        'CC', 
                        'size', 
                        'size')

# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(ipm_id,
                           'Survival',
                           if (all_pars$su_mod_ys_bestfit_index == 1) {
                             s_yr_expr <- paste("s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1)))")
                           } else if (all_pars$su_mod_ys_bestfit_index == 2) {
                             s_yr_expr <- paste("s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2)))")
                           } else if (all_pars$su_mod_ys_bestfit_index == 3) {
                             s_yr_expr <- paste("s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2 + surv_b3_yr * size_1^3)))")
                           },
                           'Evaluated',
                           'P_yr')

pdb$VitalRateExpr[2,] <- c(ipm_id,
                           'Growth',
                           if (all_pars$gr_mod_ys_bestfit_index == 1) {
                             mu_g_yr_expr <- paste("mu_g_yr = grow_b0_yr + grow_b1_yr * size_1")
                           } else if (all_pars$gr_mod_ys_bestfit_index == 2) {
                             mu_g_yr_expr <- paste("mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2")
                           } else if (all_pars$gr_mod_ys_bestfit_index == 3) {
                             mu_g_yr_expr <- paste("mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2 + grow_b3_yr * size_1^3")
                           },
                           'Evaluated',
                           'P_yr')

pdb$VitalRateExpr[3,] <- c(ipm_id,
                           'Growth',
                           'g_yr = Norm(mu_g_yr, sd_g)',
                           'Substituted',
                           'P_yr')

pdb$VitalRateExpr[4,] <- c(ipm_id,
                           'Growth',
                           'sd_g = sqrt(a * exp(b * size_1))',
                           'Evaluated',
                           'P_yr')

pdb$VitalRateExpr[5,] <- c(ipm_id,
                           'Fecundity',
                           'fy_yr = fecu_b0_yr * r_d',
                           'Evaluated',
                           'F_yr')

pdb$VitalRateExpr[6,] <- c(ipm_id,
                           'Fecundity',
                           'r_d = Norm(recr_sz, recr_sd)',
                           'Substituted',
                           'F_yr')


# Parameter Values
for(i in 1:(length(pars_var_wide))) {
  pdb$ParameterValues[i,1] <- ipm_id
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

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Growth',
    'size',
    'a',
    all_pars$a)
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Growth',
    'size',
    'b',
    all_pars$b)                               
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Fecundity',
    'size',
    'recr_sz',
    all_pars$recr_sz)
pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    'Fecundity',
    'size',
    'recr_sd',
    all_pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(
  pdb$ParameterValues$parameter_value)


# Environmental variables
pdb$ParSetIndices[1,] <- c(ipm_id,
                           'year',
                           'yr',
                           'lam_mean_ipmr$years',
                           'P_yr; F_yr',
                           '')

# Test targets
pdb$TestTargets[1:nrow(lam_mean_ipmr),1] <- ipm_id
pdb$TestTargets[1:nrow(lam_mean_ipmr),2] <- 1:nrow(lam_mean_ipmr)
pdb$TestTargets[1:nrow(lam_mean_ipmr),3] <- as.numeric(lam_mean_ipmr$value)
pdb$TestTargets[1:nrow(lam_mean_ipmr),4] <- 3

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision    <- as.numeric(pdb$TestTargets$precision)


write_xlsx(pdb, 
           paste0(data_dir, '/', script_prefix, '_', sp_abb, '_pdb_yr.xlsx'))

pdb_test       <- read_pdb(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_pdb_yr.xlsx'))

pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = 'det')
print(pdb_test_proto[[ipm_id]])
bg_ipm_pdb     <- make_ipm(pdb_test_proto[[ipm_id]])
