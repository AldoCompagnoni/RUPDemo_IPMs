# Clean data - Anderson 2016 Montana

# Author: Niklas Neisse
# Co    : Diāna Spurīte, Aspen Workman, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.04.15

# Publication: https://doi.org/10.1890/11-0193.1


# Packages ---------------------------------------------------------------------
# library(sf) #ver 1.0-1.2
# library(plantTracker) #ver 1.1.0
library(tidyverse)
library(janitor)


# Directories ------------------------------------------------------------------
dir_pub <- file.path('anderson_2016_mt')
dir_dat <- file.path(dir_pub, 'data')
dir_qud <- file.path(dir_dat, 'quadrat_data')
dir_og  <- file.path(dir_qud, 'original_data')
dir_shp <- file.path(dir_qud, 'shapefiles')


# Data -------------------------------------------------------------------------
# Species list
df_sp_list_og <- {
  path_qud <- file.path(dir_qud, 'species_list.csv')
  path_og  <- file.path(dir_og,  'species_list.csv')
  
  if (file.exists(path_og)) {
    # Load original data
    read_delim(path_og, delim = ';') %>% clean_names()
  } else {
    # Read and preserve raw content without modification
    raw_lines <- readLines(path_qud)
    dir.create(dir_og, showWarnings = FALSE, recursive = TRUE)
    writeLines(raw_lines, path_og)
    
    # Then load it as data frame
    read_delim(path_qud, delim = ';') %>% clean_names()
  }
}

# The published data seperates the genus from the species name with a ';'
#  and has thus a missmatch of column names and entries
df_sp_list <- df_sp_list_og %>% 
  {
    v_spec_col_names <- names(.)[!names(.) %in% c('x6')]
    bind_cols(
      tibble(!!v_spec_col_names[1] := str_c(.[[1]], .[[2]], sep = ' ')),
      .[, -c(1,2)]
    ) %>%
      set_names(v_spec_col_names)} 
write.csv(df_sp_list, file.path(dir_qud, 'species_list.csv'),
          row.names = F)




