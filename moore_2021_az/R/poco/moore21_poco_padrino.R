# Populating padrino - Moore 2021 Arizona - Poa compressa

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.05.14

# Publication: https://doi.org/10.1002/ecy.3661


# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker and IPM mean if the data does not exist
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the data in their respective folder


# Setting the stage ------------------------------------------------------------
# Remove all objects in the global environment
# rm(list = ls()) 


# Data -------------------------------------------------------------------------
# Define publication 
v_author_year <- c('moore_2021')
# Define region abbreviation
v_region_abb  <- c('az')
# Define species 
v_species     <- c('Poa compressa')

# A unique identifier for each model. 
#  It is 6 alphanumeric characters with no spaces
v_ipm_id      <- c('nnn510')

# IPM-Type: 'year_specific' or 'mean'?
v_ipm_type    <- c('mean')


# Taxonomic information --------------------------------------------------------
# The accepted name of the species (https://resolver.globalnames.org)
v_species_accepted <- gsub(' ', '_', v_species)
# The accepted genus
v_tax_genus  <- sub('_.*', '', v_species_accepted)
# The accepted family
v_tax_family <- c('Poaceae')
# The accepted order
v_tax_order  <- c('Poales')
# The accepted class
v_tax_class  <- c('Liliopsida')
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
v_dicot_monocot <- c('Monocot')
# Whether the species is a angiosperm or a gymosperm (only applies to plants)
v_angio_gymno   <- c('angio')


# Publication information ------------------------------------------------------
# The last names of each author on the manuscript, separated by a semicolon
v_authors  <- c('Moore; Jenness; Laughlin; Strahan; Bakker; Dowling; Springer')
# The abbreviated name of the journal that the model appears in. 
#  This follows the BIOSIS format. 
#  Exceptions are when the source is not a journal 
#  (e.g. a PhD/MSc thesis, government report). 
#  In that case, we use something like "PhD Thesis" and 
#  then include a link in the remark column
v_journal  <- c('Ecology')
#  The year the article was published
v_pub_year <- c('2022')
# The DOI of the publication (NOT THE doi.org URL though!!)
v_doi      <- c('10.1002/ecy.3661')
# The last name of the corresponding author
v_corresponding_author <- c('Moore')
# The corresponding author’s email, along with the year of publication 
#  in parentheses to denote how old (and possibly inaccessible) it is. 
#  For example, this could levisc8@gmail.com (2020). 
#  If you are able to find a more recent email address via Google, 
#  then this can also be used (this isn’t necessarily expected though).
v_email_year <- c('margaret.moore@nau.edu (2025)')
# Any qualitative comments you may have on the model. 
#  These can range from comments to accuracy of GPS coordinates to descriptions 
#  of the different levels of a treatment that was applied
v_remark   <- c(NA)
# The full APA style citation for the paper
v_apa_citation <- c('Moore, M. M., Jenness, J. S., Laughlin, D. C., Strahan, R. T., Bakker, J. D., Dowling, H. E., & Springer, J. D. (2022). Cover and density of southwestern ponderosa pine understory plants in permanent chart quadrats (2002–2020). Ecology, 103(5), e3661.')
# If there is one, a link to the Electronic Supplementary Material that 
#  contains further details/parameter values for the model
v_demog_appendix_link <- c('https://www.fs.usda.gov/rds/archive/catalog/RDS-2021-0092')


# Data collection information --------------------------------------------------
# The year that demographic data collection began. Formatted YYYY (e.g. 1990)
v_start_year  <- c(2002)
# The month of the year that demographic data collection began. 
#  This is an integer between 1 and 12, where 1 corresponds to January
v_start_month <- c(NA)
#  The final year of demographic data collection. Formatted YYYY
v_end_year    <- c(2022)
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
v_lat         <- c('35.184')
# The decimal longitude of the population. 
#  Use the dms_deg function from pdbDigitUtils to generate this
v_lon         <- c('-111.767')
# The altitude above/below sea level, in meters
v_altitude    <- c('2207.4')
# The ISO3 country code for the country in which the population is located. 
v_country     <- c('USA')
# The continent that the population is located on. 
#  Options are n_america, s_america, oceania, asia, europe and africa. 
#  Others may be added as needed
v_continent   <- c('n_america')
# The biome code
#  https://patrickbarks.shinyapps.io/biomes/
v_ecoregion   <- c('TCF')


# Main code --------------------------------------------------------------------
source('pipeline/padrino_ipm_type.R')


# Parameters -------------------------------------------------------------------
# All parameters of the ipm  
pars

# Year specific lambda -impr- 
lam_mean_ipmr

# Padrino entry
pdb


# Run IPM with padrino --------------------------------------------------------- 
# Deterministic lambda
bg_ipm_pdb
lambda(bg_ipm_pdb)
lam_mean_ipmr

# Testing the model with padrino
test_model(pdb_test, id = v_ipm_id)
