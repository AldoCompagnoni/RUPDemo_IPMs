# Populating padrino - Anderson 2016 Arizona - Artemisia frigida

# Author: Diāna Spurīte
# Co    : Aspen Workman, Aldo Compagnoni, Niklas Neisse
# Email : diana.spurite@posteo.de
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.06.18

# Publication: https://doi.org/10.1890/11-0193.1


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker and IPM mean if the data does not exist
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the data in their respective folder


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 


# Data -------------------------------------------------------------------------
# Define publication 
v_author_year <- c('anderson_2016')
# Define region abbreviation
v_region_abb  <- c('mt')
# Define species 
v_species     <- c('Artemisia frigida')

# A unique identifier for each model. 
#  It is 6 alphanumeric characters with no spaces
v_ipm_id      <- c('afr162')

# IPM-Type: 'year_specific' or 'mean'?
v_ipm_type    <- c('mean')


# Taxonomic information --------------------------------------------------------
# The accepted name of the species (https://resolver.globalnames.org)
v_species_accepted <- gsub(' ', '_', v_species)
# The accepted genus
v_tax_genus  <- sub('_.*', '', v_species_accepted)
# The accepted family
v_tax_family <- c('Asteraceae')
# The accepted order
v_tax_order  <- c('Asterales')
# The accepted class
v_tax_class  <- c('Magnoliopsida')
# The accepted phylum
v_tax_phylum <- c('Tracheophyta')
# The kingdom
v_kingdom    <- c('Plantae')
# The type of organism. For plants, this is usually something like 
#  "Herbaceous perennial", or "Shrub". For animals, this could be, for example, 
#  "mammal" or "reptile". See here for more details 
#   (but also do not hesitate to contact me if there instances that 
#   fall outside of the classification given there)
v_organism_type <- c('Herbaceous')
# Whether the species is a dicotyledon or a monocotyledon 
#  (only applies to plants)
v_dicot_monocot <- c('Dicot')
# Whether the species is a angiosperm or a gymosperm (only applies to plants)
v_angio_gymno   <- c('angio')


# Publication information ------------------------------------------------------
# The last names of each author on the manuscript, separated by a semicolon
v_authors  <- c('Anderson; Vermeire; Adler')
# The abbreviated name of the journal that the model appears in. 
#  This follows the BIOSIS format. 
#  Exceptions are when the source is not a journal 
#  (e.g. a PhD/MSc thesis, government report). 
#  In that case, we use something like "PhD Thesis" and 
#  then include a link in the remark column
v_journal  <- c('Ecology')
#  The year the article was published
v_pub_year <- c('2016')
# The DOI of the publication (NOT THE doi.org URL though!!)
v_doi      <- c('10.6084/m9.figshare.c.3304113.v1')
# The last name of the corresponding author
v_corresponding_author <- c('Adler')
# The corresponding author’s email, along with the year of publication 
#  in parentheses to denote how old (and possibly inaccessible) it is. 
#  For example, this could levisc8@gmail.com (2020). 
#  If you are able to find a more recent email address via Google, 
#  then this can also be used (this isn’t necessarily expected though).
v_email_year <- c('peter.adler@usu.edu (2024)')
# Any qualitative comments you may have on the model. 
#  These can range from comments to accuracy of GPS coordinates to descriptions 
#  of the different levels of a treatment that was applied
v_remark   <- c('cattle grazing treatments with light, moderate, and heavy stocking rates of 1.24, 0.92, and 0.76 ha/animal-unit-month (two pastures in each)')
# The full APA style citation for the paper
v_apa_citation <- c('Anderson, J., Vermeire, L., & Adler, P. B. (2011). Fourteen years of mapped, permanent quadrats in a northern mixed prairie, USA. Ecology, 92(8), 1703. https://esapubs.org/archive')
# If there is one, a link to the Electronic Supplementary Material that 
#  contains further details/parameter values for the model
v_demog_appendix_link <- c('https://figshare.com/collections/Fourteen_years_of_mapped_permanent_quadrats_in_a_northern_mixed_prairie_USA/3304113')


# Data collection information --------------------------------------------------
# The year that demographic data collection began. Formatted YYYY (e.g. 1990)
v_start_year  <- c(1932)
# The month of the year that demographic data collection began. 
#  This is an integer between 1 and 12, where 1 corresponds to January
v_start_month <- c(NA)
#  The final year of demographic data collection. Formatted YYYY
v_end_year    <- c(1945)
# The month of the year that demographic data collection concluded
v_end_month   <- c(NA)
# Indicates the time step (periodicity) for which the seasonal, annual, 
#  or multi-annual IPM was constructed. For example, 1 indicates that 
#  the IPM iteration period is 1 year; 
#  0.5 indicates that the IPM iterates once every 0.5 years or 6 months; 
#  2 indicates that the IPM iteration occurs every 2 years
v_periodicity <- c(1)
# The name of the population given by the author. 
#  For example, "Bear Creek", or "Havatselet". 
#  If the population names are missing, 
#  use sequential names in alphabetical order (e.g. "A", "B", "C", etc).
v_population_name <- c(NA)
# Sometimes, a population_name may encompass multiple sub-populations that 
#  are located close by. This integer specifies the number of 
#  populations/sub-populations that are described by the model.
v_number_populations <- c(NA)
# The decimal latitude of the population. 
#  Use the dms_deg function from pdbDigitUtils to generate this
v_lat         <- c('46.4')
# The decimal longitude of the population. 
#  Use the dms_deg function from pdbDigitUtils to generate this
v_lon         <- c('-105.7')
# The altitude above/below sea level, in meters
v_altitude    <- c('1250')
# The ISO3 country code for the country in which the population is located. 
v_country     <- c('USA')
# The continent that the population is located on. 
#  Options are n_america, s_america, oceania, asia, europe and africa. 
#  Others may be added as needed
v_continent   <- c('n_america')
# The biome code
#  https://patrickbarks.shinyapps.io/biomes/
v_ecoregion   <- c('DES')


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
test_model(pdb_test, id = v_ipm_id)
