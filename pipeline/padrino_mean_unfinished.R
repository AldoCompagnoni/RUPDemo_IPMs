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
R_dir      <- file.path(pub_dir, 'R',       sp_abb)
data_dir   <- file.path(pub_dir, 'data',    sp_abb)
result_dir <- file.path(pub_dir, 'results', sp_abb)

# Prefix for all the files
script_prefix <- str_c(str_extract(author_year, '^[^_]+'), 
                       str_sub(str_extract(author_year, '_\\d+$'), -2, -1))

# Ipm mean and plant tracker if they not already exists
if (!file.exists(
  paste0(data_dir, '/', script_prefix, '_', sp_abb, '_data_df.csv'))) {
  source(paste0(R_dir, '/', script_prefix, '_', sp_abb, '_ipm_mean.R'))
}

pars          <- read.csv(paste0(data_dir, '/pars.csv'))
lam_mean_ipmr <- read.csv(paste0(data_dir, '/lambda.csv'))


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


