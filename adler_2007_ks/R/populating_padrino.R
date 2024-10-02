# Populating the PADRINO database template -------------------------------------
# Store the PADRINO excel template
YOUR_PATH <- "."
sheet_names <- excel_sheets( paste( YOUR_PATH, "/pdb_template.xlsx", sep = "" ) )
pdb <- lapply( sheet_names, function( x ) {
  as.data.frame( read_excel( paste( YOUR_PATH, "/pdb_template.xlsx", sep = "" ), sheet = x ) ) } )
names( pdb ) <- sheet_names


pdb$Metadata[1,] <- c( # 1 ID
  "nnnnn1", 
  
  # Taxonomic information: 2-12
  "Bouteloua_curtipendula", "Bouteloua_curtipendula", "Bouteloua",
  "Poaceae", "Poales", "Liliopsida", "Magnoliophyta",
  "Plantae", "Herbaceous", "Monocot", "angio", 
  
  # Publication information: 13-21
  "Chu; Tyburczy; Laurenroth",
  "Ecology", "2007", "10.1890/0012-9658(2007)88[2673%3ALMQFKP]2.0.CO%3B2", "Adler", 
  "peter.adler@usu.edu (2024)", NA, 
  "Adler, P. B., Tyburczy, W. R., & Lauenroth, W. K. (2007). LONG‐TERM MAPPED QUADRATS FROM KANSAS PRAIRIE DEMOGRAPHIC INFORMATION FOR HERBACEOUS PLANTS: Ecological Archives E088‐161. Ecology, 88(10), 2673-2673.",
  "https://figshare.com/collections/LONG-TERM_MAPPED_QUADRATS_FROM_KANSAS_PRAIRIE_DEMOGRAPHIC_INFORMATION_FOR_HERBACEOUS_PLANTS/3299993",
  
  # Data collection information : 
  40, 1932, 6, 2071, 6, 1, NA, NA, 
  "38.8", "-99.3", "462", "USA",
  "n_america", "TGS",
  
  # Model information
  "A", TRUE, "truncated_distributions", "P_yr; F_yr", NA, FALSE,
  FALSE, FALSE, FALSE, "", "", ""
)

pdb$Metadata$eviction_used <- as.logical(pdb$Metadata$eviction_used)
pdb$Metadata$duration <- as.numeric(pdb$Metadata$duration)
pdb$Metadata$periodicity <- as.numeric(pdb$Metadata$periodicity)