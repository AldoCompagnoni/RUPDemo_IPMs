# Populating padrino - Zachmann 2016 Idaho - Pseudoroegneria spicata

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.01.13

# Publication: https://doi.org/10.1890/10-0404.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker and IPM mean if the data does not exist
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the data in their respective folder


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 


# Data -------------------------------------------------------------------------
# Define publication 
author_year <- 'zachmann_2016'
# Define region abbreviation
region_abb  <- 'id'
# Define species 
species     <- 'Pseudoroegneria spicata'

# A unique identifier for each model. 
#  It is 6 alphanumeric characters with no spaces
ipm_id      <- 'nnnnn17'

# IPM-Type: 'year_specific' or 'mean'?
ipm_type    <- 'year_specific'


# Taxonomic information --------------------------------------------------------
# The accepted name of the species (here Wikipedia)
species_accepted <- gsub(" ", "_", species)
# The accepted genus
tax_genus  <- sub('_.*', '', species_accepted)
# The accepted family
tax_family <- 'Poaceae'
# The accepted order
tax_order  <- 'Poales' 
# The accepted class
tax_class  <- 'Liliopsida'
# The accepted phylum
tax_phylum <- 'Magnoliophyta' 
# The kingdom
kingdom    <- 'Plantae'
# The type of organism. For plants, this is usually something like 
#  "Herbaceous perennial", or "Shrub". For animals, this could be, for example, 
#  "mammal" or "reptile". See here for more details 
#   (but also do not hesitate to contact me if there instances that 
#   fall outside of the classification given there)
organism_type <- 'Herbaceous' 
# Whether the species is a dicotyledon or a monocotyledon 
#  (only applies to plants)
dicot_monocot <- 'Monocot'
# Whether the species is a angiosperm or a gymosperm (only applies to plants)
angio_gymno   <- 'angio'


# Publication information ------------------------------------------------------
# The last names of each author on the manuscript, separated by a semicolon
authors  <- 'Zachmann; Moffet; Adler'
# The abbreviated name of the journal that the model appears in. 
#  This follows the BIOSIS format. 
#  Exceptions are when the source is not a journal 
#  (e.g. a PhD/MSc thesis, government report). 
#  In that case, we use something like "PhD Thesis" and 
#  then include a link in the remark column
journal  <- 'Ecology'
#  The year the article was published
pub_year <- '2016'
# The DOI of the publication (NOT THE doi.org URL though!!)
doi      <- '10.1890/10-0404.1'
# The last name of the corresponding author
corresponding_author <- 'Zachmann'
# The corresponding author’s email, along with the year of publication 
#  in parentheses to denote how old (and possibly inaccessible) it is. 
#  For example, this could levisc8@gmail.com (2020). 
#  If you are able to find a more recent email address via Google, 
#  then this can also be used (this isn’t necessarily expected though).
email_year <- 'lzachmann@gmail.com (2024)'
# Any qualitative comments you may have on the model. 
#  These can range from comments to accuracy of GPS coordinates to descriptions 
#  of the different levels of a treatment that was applied
remark   <- NA
# The full APA style citation for the paper
apa_citation <- 'Zachmann, L., Moffet, C., & Adler, P. (2010). Mapped quadrats in sagebrush steppe: long‐term data for analyzing demographic rates and plant–plant interactions: Ecological Archives E091‐243. Ecology, 91(11), 3427-3427..'
# If there is one, a link to the Electronic Supplementary Material that 
#  contains further details/parameter values for the model
demog_appendix_link <- 'https://figshare.com/collections/Mapped_quadrats_in_sagebrush_steppe_long-term_data_for_analyzing_demographic_rates_and_plant_plant_interactions/3303612'


# Data collection information --------------------------------------------------
# The year that demographic data collection began. Formatted YYYY (e.g. 1990)
start_year  <- 1923
# The month of the year that demographic data collection began. 
#  This is an integer between 1 and 12, where 1 corresponds to January
start_month <- 6
#  The final year of demographic data collection. Formatted YYYY
end_year    <- 1957
# The month of the year that demographic data collection concluded
end_month   <- 6
# Indicates the time step (periodicity) for which the seasonal, annual, 
#  or multi-annual IPM was constructed. For example, 1 indicates that 
#  the IPM iteration period is 1 year; 
#  0.5 indicates that the IPM iterates once every 0.5 years or 6 months; 
#  2 indicates that the IPM iteration occurs every 2 years
periodicity <- 1
# The name of the population given by the author. 
#  For example, "Bear Creek", or "Havatselet". 
#  If the population names are missing, 
#  use sequential names in alphabetical order (e.g. "A", "B", "C", etc).
population_name <- NA
# Sometimes, a population_name may encompass multiple sub-populations that 
#  are located close by. This integer specifies the number of 
#  populations/sub-populations that are described by the model.
number_populations <- NA
# The decimal latitude of the population. 
#  Use the dms_deg function from pdbDigitUtils to generate this
lat         <- '44.2'
# The decimal longitude of the population. 
#  Use the dms_deg function from pdbDigitUtils to generate this
lon         <- '-112.1'
# The altitude above/below sea level, in meters
altitude    <- '1569'
# The ISO3 country code for the country in which the population is located. 
country     <- 'USA'
# The continent that the population is located on. 
#  Options are n_america, s_america, oceania, asia, europe and africa. 
#  Others may be added as needed
continent   <- 'n_america'
# The biome code
#  https://patrickbarks.shinyapps.io/biomes/
ecoregion   <- 'DES'


# Main code --------------------------------------------------------------------
source('pipeline/padrino_ipm_type.R')


# Parameters -------------------------------------------------------------------
# All parameters of the ipm  
all_pars

# Year specific lambda -impr- 
lam_mean_ipmr

# Padrino entry
pdb


# Run IPM with padrino --------------------------------------------------------- 
# Deterministic lambda, year specific
bg_ipm_pdb
lambda(bg_ipm_pdb)
plot(lambda(bg_ipm_pdb) ~ lam_mean_ipmr$years)

# Testing the model with padrino
test_model(pdb_test, id = ipm_id)
