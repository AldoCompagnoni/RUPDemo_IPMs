# Populating padrino - alder 2007; kansas; Sporobolus asper

# Niklas Neisse
# 2024.10.24

#


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)
options(stringsAsFactors = F)


# Packages ---------------------------------------------------------------------
# Define CRAN packages
.cran_packages <- c("tidyverse","ipmr","readxl", "writexl", "remotes") 
# Check if CRAN packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}
# Load required packages
sapply(.cran_packages, require, character.only = TRUE) 

# remotes::install_github("padrinoDB/pdbDigitUtils",force = TRUE)
library(pdbDigitUtils)

rm(list = ls())
options(stringsAsFactors = F)


# Data -------------------------------------------------------------------------
# Define the species variable
species       <- "Sporobolus asper"
sp_abb        <- tolower(
  gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2),
                      collapse = "")))

# # plantTracker 
# source(paste0('adler_2007_ks/R/', sp_abb, '/kansas_tracker_', sp_abb, '.R'))
# # overall imp 
# source(paste0('adler_2007_ks/R/', sp_abb, '/ipm_tracker_', sp_abb, '.R'))

ipm_id        <- "nnnnnn9"
pars          <- read.csv(
  paste0("adler_2007_ks/data/", sp_abb, "/pars.csv"))
lam_mean_ipmr <- read.csv(
  paste0("adler_2007_ks/data/", sp_abb, "/lambda.csv"))


# Populating the PADRINO database template -------------------------------------
# Store the PADRINO excel template
YOUR_PATH   <- "."
sheet_names <- excel_sheets(paste(YOUR_PATH, "/pdb_template.xlsx", sep = ""))
pdb         <- lapply(sheet_names, function(x) {
  as.data.frame(read_excel(paste(YOUR_PATH, "/pdb_template.xlsx", sep = ""), sheet = x))})
names(pdb)  <- sheet_names


pdb$Metadata[1,] <- c(
  ipm_id, 
  
  # Taxonomic information
  "Panicum_virgatum", "Panicum_virgatum", "Panicum",
  "Poaceae", "Poales", "Liliopsida", "Magnoliophyta",
  "Plantae", "Herbaceous", "Monocot", "angio", 
  
  # Publication information
  "Adler; Tyburczy; Laurenroth",
  "Ecology", "2007", "10.1890/0012-9658(2007)88[2673%3ALMQFKP]2.0.CO%3B2", "Adler", 
  "peter.adler@usu.edu (2024)", NA, 
  "Adler, P. B., Tyburczy, W. R., & Lauenroth, W. K. (2007). LONG‐TERM MAPPED QUADRATS FROM KANSAS PRAIRIE DEMOGRAPHIC INFORMATION FOR HERBACEOUS PLANTS: Ecological Archives E088‐161. Ecology, 88(10), 2673-2673.",
  "https://figshare.com/collections/LONG-TERM_MAPPED_QUADRATS_FROM_KANSAS_PRAIRIE_DEMOGRAPHIC_INFORMATION_FOR_HERBACEOUS_PLANTS/3299993",
  
  # Data collection information
  40, 1932, 6, 1971, 6, 1, NA, NA, 
  "38.8", "-99.3", "462", "USA",
  "n_america", "TGS",
  
  # Model information
  "A", TRUE, "truncated_distributions", NA, NA, FALSE,
  FALSE, FALSE, FALSE, "", "", ""
)

pdb$Metadata$eviction_used  <- as.logical(pdb$Metadata$eviction_used)
pdb$Metadata$duration       <- as.numeric(pdb$Metadata$duration)
pdb$Metadata$periodicity    <- as.numeric(pdb$Metadata$periodicity)

pdb$StateVariables[1,]      <- c(ipm_id, "size", FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(ipm_id, 
                               "size", 
                               "", 
                               pars$L, 
                               pars$U, 
                               "P; F", 
                               "")

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(ipm_id,
                              "size",
                              "",
                              pars$mat_siz,
                              "midpoint",
                              "P; F")

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(ipm_id,
                          "n_size",
                          pars$mat_siz,
                          "")
pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(ipm_id, 
                        "P", 
                        "P = s * g * d_size", 
                        "CC", 
                        "size", 
                        "size")

pdb$IpmKernels[2,] <- c(ipm_id, 
                        "F", 
                        "F= fy * d_size", 
                        "CC", 
                        "size", 
                        "size")

# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(ipm_id,
                           "Survival",
                           "s = 1 / (1 + exp(-(surv_b0 + surv_b1 * size_1)))",
                           "Evaluated",
                           "P")

pdb$VitalRateExpr[2,] <- c(ipm_id,
                           "Growth",
                           "mu_g = grow_b0 + grow_b1 * size_1 + grow_b2 * size_1^2",
                           "Evaluated",
                           "P")

pdb$VitalRateExpr[3,] <- c(ipm_id,
                           "Growth",
                           "g = Norm(mu_g, sd_g)",
                           "Substituted",
                           "P")

pdb$VitalRateExpr[4,] <- c(ipm_id,
                           "Growth",
                           "sd_g = sqrt(a * exp(b * size_1))",
                           "Evaluated",
                           "P")

pdb$VitalRateExpr[5,] <- c(ipm_id,
                           "Fecundity",
                           "fy = fecu_b0 * r_d",
                           "Evaluated",
                           "F")

pdb$VitalRateExpr[6,] <- c(ipm_id,
                           "Fecundity",
                           "r_d = Norm(recr_sz, recr_sd)",
                           "Substituted",
                           "F")


# CHECK -- Parameter values ####
# Parameter Values
pdb$ParameterValues[1,] <- c(ipm_id,
                             "Survival",
                             "size",
                             "surv_b0",
                             pars$surv_b0)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Survival",
                             "size",
                             "surv_b1",
                             pars$surv_b1)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Growth",
                             "size",
                             "grow_b0",
                             pars$grow_b0)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Growth",
                             "size",
                             "grow_b1",
                             pars$grow_b1)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Growth",
                             "size",
                             "grow_b2",
                             pars$grow_b2)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Growth",
                             "size",
                             "a",
                             pars$a)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Growth",
                             "size",
                             "b",
                             pars$b)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Fecundity",
                             "size",
                             "fecu_b0",
                             pars$fecu_b0)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Fecundity",
                             "size",
                             "recr_sz",
                             pars$recr_sz)

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- c(ipm_id,
                             "Fecundity",
                             "size",
                             "recr_sd",
                             pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(
  pdb$ParameterValues$parameter_value)


# Test targets
pdb$TestTargets[1,] <- c(ipm_id,
                         "lambda",
                         lam_mean_ipmr,
                         3)

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision <- as.numeric(pdb$TestTargets$precision)


write_xlsx(pdb, 
           paste0("adler_2007_ks/data/", sp_abb, "/", sp_abb, "_pdb.xlsx"))

pdb_test       <- read_pdb(
  paste0("adler_2007_ks/data/", sp_abb, "/", sp_abb, "_pdb.xlsx"))

pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = "det")
print(pdb_test_proto[[ipm_id]])
bg_ipm_pdb     <- make_ipm(pdb_test_proto[[ipm_id]])

bg_ipm_pdb
lambda(bg_ipm_pdb)
test_model(pdb_test, id = ipm_id)
