# IPM for kansas bogr year specific

# Niklas Neisse
# 2024.09.23

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


rm( list = ls() )
options( stringsAsFactors = F )




# Data -------------------------------------------------------------------------
# Define the species variable
species        <- "Bouteloua gracilis"
sp_abb         <- tolower(gsub(" ", "", paste(substr(unlist(strsplit(species, " ")), 1, 2), 
                                              collapse = "")))
ipm_id         <- "nnnnn2"
                   
all_pars       <- read.csv(paste0("adler_2007_ks/data/", sp_abb, "/all_pars.csv"))
pars_var_wide  <- read.csv(paste0("adler_2007_ks/data/", sp_abb, "/2.pars_var.csv"))
lam_mean_ipmr  <- read.csv(paste0("adler_2007_ks/data/", sp_abb, "/lambdas_yr_vec.csv"))
# lam_out        <- read.csv(paste0("adler_2007_ks/data/", sp_abb, "/lambdas_yr_df.csv"))


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
  "Bouteloua_gracilis", "Bouteloua_gracilis", "Bouteloua",
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

pdb$Metadata$eviction_used <- as.logical(pdb$Metadata$eviction_used)
pdb$Metadata$duration <- as.numeric(pdb$Metadata$duration)
pdb$Metadata$periodicity <- as.numeric(pdb$Metadata$periodicity)

pdb$StateVariables[1,] <- c(ipm_id, "size", FALSE)
pdb$StateVariables$discrete <- as.logical(pdb$StateVariables$discrete)

pdb$ContinuousDomains[1,] <- c(ipm_id, 
                               "size", 
                               "", 
                               all_pars$L, 
                               all_pars$U, 
                               "P_yr; F_yr", 
                               "")

pdb$ContinuousDomains$lower <- as.numeric(pdb$ContinuousDomains$lower)
pdb$ContinuousDomains$upper <- as.numeric(pdb$ContinuousDomains$upper)

pdb$IntegrationRules[1,] <- c(ipm_id,
                              "size",
                              "",
                              all_pars$mat_siz,
                              "midpoint",
                              "P_yr; F_yr")

pdb$IntegrationRules$n_meshpoints <- as.numeric(pdb$IntegrationRules$n_meshpoints)

pdb$StateVectors[1,] <- c(ipm_id,
                          "n_size",
                          all_pars$mat_siz,
                          "" )
pdb$StateVectors$n_bins <- as.numeric(pdb$StateVectors$n_bins)


# Ipm Kernels
pdb$IpmKernels[1,] <- c(ipm_id, 
                        "P_yr", 
                        "P_yr = s_yr * g_yr * d_size", 
                        "CC", 
                        "size", 
                        "size")

pdb$IpmKernels[2,] <- c(ipm_id, 
                        "F_yr", 
                        "F_yr = fy_yr * d_size", 
                        "CC", 
                        "size", 
                        "size")

# Vital rate expressions
pdb$VitalRateExpr[1,] <- c(ipm_id,
                           "Survival",
                           "s_yr = 1 / (1 + exp(-(surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2)))",
                           "Evaluated",
                           "P_yr")

pdb$VitalRateExpr[2,] <- c(ipm_id,
                           "Growth",
                           "mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2",
                           "Evaluated",
                           "P_yr")

pdb$VitalRateExpr[3,] <- c(ipm_id,
                            "Growth",
                            "g_yr = Norm(mu_g_yr, sd_g)",
                            "Substituted",
                            "P_yr")

pdb$VitalRateExpr[4,] <- c(ipm_id,
                           "Growth",
                           "sd_g = sqrt(a * exp(b * size_1))",
                           "Evaluated",
                           "P_yr")

pdb$VitalRateExpr[5,] <- c(ipm_id,
                           "Fecundity",
                           "fy_yr = fecu_b0_yr * r_d",
                           "Evaluated",
                           "F_yr")

pdb$VitalRateExpr[6,] <- c(ipm_id,
                           "Fecundity",
                           "r_d = Norm(recr_sz, recr_sd)",
                           "Substituted",
                           "F_yr")


# Parameter Values
for(i in 1:(length(pars_var_wide))) {
  pdb$ParameterValues[i,1] <- ipm_id
  pdb$ParameterValues[i,3] <- "size"
  pdb$ParameterValues[i,4] <- names(pars_var_wide)[i]
  pdb$ParameterValues[i,5] <- as.numeric(pars_var_wide[i])
  
  if(grepl("surv", names(pars_var_wide)[i])){
    pdb$ParameterValues[i,2] <- "Survival"
  } else {
    if(grepl("grow", names(pars_var_wide)[i])){
      pdb$ParameterValues[i,2] <- "Growth"
    } else {pdb$ParameterValues[i,2] <- "Fecundity"}
  }
}

pdb$ParameterValues[nrow(pdb$ParameterValues)+1,] <- 
  c(ipm_id,
    "Growth",
    "size",
    "a",
    all_pars$a)
pdb$ParameterValues[nrow(pdb$ParameterValues)+2,] <- 
  c(ipm_id,
     "Growth",
     "size",
     "b",
     all_pars$b)                               
pdb$ParameterValues[nrow(pdb$ParameterValues)+3,] <- 
  c(ipm_id,
     "Fecundity",
     "size",
     "recr_sz",
     all_pars$recr_sz)
pdb$ParameterValues[nrow(pdb$ParameterValues)+4,] <- 
  c(ipm_id,
     "Fecundity",
     "size",
     "recr_sd",
     all_pars$recr_sd)

pdb$ParameterValues$parameter_value <- as.numeric(pdb$ParameterValues$parameter_value)


# Environmental variables
pdb$ParSetIndices[1,] <- c(ipm_id,
                           "year",
                           "yr",
                           "34:71",
                           "P_yr; F_yr",
                           "")

# Test targets
pdb$TestTargets[1:nrow(lam_mean_ipmr),1] <- ipm_id
pdb$TestTargets[1:nrow(lam_mean_ipmr),2] <- 1:nrow(lam_mean_ipmr)
pdb$TestTargets[1:nrow(lam_mean_ipmr),3] <- as.numeric(lam_mean_ipmr$value)
pdb$TestTargets[1:nrow(lam_mean_ipmr),4] <- 3

pdb$TestTargets$target_value <- as.numeric(pdb$TestTargets$target_value)
pdb$TestTargets$precision    <- as.numeric(pdb$TestTargets$precision)


write_xlsx(pdb, paste0("adler_2007_ks/data/", sp_abb, "/bou_gra_yr_pdb.xlsx"))


pdb_test       <- read_pdb(paste0("adler_2007_ks/data/", sp_abb, "/bou_gra_yr_pdb.xlsx"))
pdb_test_proto <- pdb_make_proto_ipm(pdb_test, det_stoch = "det")
print(pdb_test_proto$nnnnn2)
bg_ipm_pdb     <- make_ipm(pdb_test_proto$nnnnn2)

bg_ipm_pdb
lambda(bg_ipm_pdb)
test_model(pdb_test, id = ipm_id)
