# IPM for kansas bocu year specific

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
.cran_packages <- c("tidyverse","ipmr","readxl", "writexl") 
# Check if CRAN packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}
# Load required packages
sapply(.cran_packages, require, character.only = TRUE) 

rm( list = ls() )
options( stringsAsFactors = F )


# Data -------------------------------------------------------------------------
all_pars       <- read.csv("adler_2007_ks/data/bocu/all_pars.csv")
lam_out        <- read.csv("adler_2007_ks/data/bocu/lambdas_yr_df.csv")
pars_var_wide  <- read.csv("adler_2007_ks/data/bocu/2.pars_var.csv")
lam_mean_ipmr  <- read.csv("adler_2007_ks/data/bocu/lambdas_yr_vec.csv")


# Populating the PADRINO database template -------------------------------------
# Store the PADRINO excel template
YOUR_PATH <- "."
sheet_names <- excel_sheets( paste( YOUR_PATH, "/pdb_template.xlsx", sep = "" ) )
pdb <- lapply( sheet_names, function( x ) {
  as.data.frame( read_excel( paste( YOUR_PATH, "/pdb_template.xlsx", sep = "" ), sheet = x ) ) } )
names( pdb ) <- sheet_names


pdb$Metadata[1,] <- c( # 1 ID
  "nnnnn1", 
  
  # Taxonomic information
  "Bouteloua_curtipendula", "Bouteloua_curtipendula", "Bouteloua",
  "Poaceae", "Poales", "Liliopsida", "Magnoliophyta",
  "Plantae", "Herbaceous", "Monocot", "angio", 
  
  # Publication information
  "Chu; Tyburczy; Laurenroth",
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

pdb$StateVariables[1,] <- c( "nnnnn1", "size", FALSE)
pdb$StateVariables$discrete <- as.logical( pdb$StateVariables$discrete )

pdb$ContinuousDomains[1,] <- c( "nnnnn1", 
                                "size", 
                                "", 
                                all_pars$L, 
                                all_pars$U, 
                                "P_yr; F_yr", 
                                "" )

pdb$ContinuousDomains$lower <- as.numeric( pdb$ContinuousDomains$lower )
pdb$ContinuousDomains$upper <- as.numeric( pdb$ContinuousDomains$upper )

pdb$IntegrationRules[1,] <- c( "nnnnn1",
                               "size",
                               "",
                               all_pars$mat_siz,
                               "midpoint",
                               "P_yr; F_yr" )

pdb$IntegrationRules$n_meshpoints <- as.numeric( pdb$IntegrationRules$n_meshpoints )

pdb$StateVectors[1,] <- c( "nnnnn1",
                           "n_size",
                           all_pars$mat_siz,
                           "" )
pdb$StateVectors$n_bins <- as.numeric( pdb$StateVectors$n_bins )


# Ipm Kernels
pdb$IpmKernels[1,] <- c( "nnnnn1", 
                         "P_yr", 
                         "P_yr = s_yr * g_yr * d_size", 
                         "CC", 
                         "size", 
                         "size" )

pdb$IpmKernels[2,] <- c( "nnnnn1", 
                         "F_yr", 
                         "F_yr = fy_yr * d_size", 
                         "CC", 
                         "size", 
                         "size" )

# Vital rate expressions
pdb$VitalRateExpr[1,] <- c( "nnnnn1",
                            "Survival",
                            "s_yr = 1 / ( 1 + exp( -( surv_b0_yr + surv_b1_yr * size_1 + surv_b2_yr * size_1^2 ) ) )",
                            "Evaluated",
                            "P_yr" )

pdb$VitalRateExpr[2,] <- c( "nnnnn1",
                            "Growth",
                            "mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2 + grow_b3_yr * size_1^3",
                            "Evaluated",
                            "P_yr" )

pdb$VitalRateExpr[3,] <- c( "nnnnn1",
                            "Growth",
                            "g_yr = Norm( mu_g_yr, sd_g )",
                            "Substituted",
                            "P_yr" )

pdb$VitalRateExpr[4,] <- c( "nnnnn1",
                            "Growth",
                            "sd_g = sqrt( a * exp( b * size_1 ) )",
                            "Evaluated",
                            "P_yr" )

pdb$VitalRateExpr[5,] <- c( "nnnnn1",
                            "Fecundity",
                            "fy_yr = fecu_b0_yr * r_d",
                            "Evaluated",
                            "F_yr" )

pdb$VitalRateExpr[6,] <- c( "nnnnn1",
                            "Fecundity",
                            "r_d = Norm( recr_sz, recr_sd )",
                            "Substituted",
                            "F_yr" )


# Parameter Values
for( i in 1:( length( pars_var_wide ) ) ) {
  pdb$ParameterValues[i,1] <- "nnnnn1"
  pdb$ParameterValues[i,3] <- "size"
  pdb$ParameterValues[i,4] <- names( pars_var_wide )[i]
  pdb$ParameterValues[i,5] <- as.numeric( pars_var_wide[i] )
  
  if( grepl( "surv", names( pars_var_wide )[i] ) ){
    pdb$ParameterValues[i,2] <- "Survival"
  } else {
    if( grepl( "grow", names( pars_var_wide )[i] ) ){
      pdb$ParameterValues[i,2] <- "Growth"
    } else { pdb$ParameterValues[i,2] <- "Fecundity" }
  }
}

pdb$ParameterValues[79,] <- c( "nnnnn1",
                               "Growth",
                               "size",
                               "a",
                               all_pars$a )
pdb$ParameterValues[80,] <- c( "nnnnn1",
                               "Growth",
                               "size",
                               "b",
                               all_pars$b )                               
pdb$ParameterValues[81,] <- c( "nnnnn1",
                               "Fecundity",
                               "size",
                               "recr_sz",
                               all_pars$recr_sz )
pdb$ParameterValues[82,] <- c( "nnnnn1",
                               "Fecundity",
                               "size",
                               "recr_sd",
                               all_pars$recr_sd )

pdb$ParameterValues$parameter_value <- as.numeric( pdb$ParameterValues$parameter_value )


# Environmental variables
pdb$ParSetIndices[1,] <- c( "nnnnn1",
                            "year",
                            "yr",
                            "1932:1971",
                            "P_yr; F_yr",
                            "" )

# Test targets
pdb$TestTargets[1:13,1] <- "nnnnn1"
pdb$TestTargets[1:13,2] <- names(lam_mean_ipmr)
pdb$TestTargets[1:13,3] <- as.numeric(lam_mean_ipmr)
pdb$TestTargets[1:13,4] <- 3

pdb$TestTargets$target_value <- as.numeric( pdb$TestTargets$target_value )
pdb$TestTargets$precision <- as.numeric( pdb$TestTargets$precision )


write_xlsx( pdb, "adler_2007_ks/data/bocu/bou_cur_yr_pdb.xlsx" )
