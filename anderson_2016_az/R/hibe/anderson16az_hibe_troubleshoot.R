library(plotly)
library(shiny)
library(gridExtra)

# Specifications ---------------------------------------------------------------
# Define publication 
v_author_year      <- c('anderson_2016')
# Define region abbreviation
v_region_abb       <- c('az')
# Define growth form (grass, forb, shrub, c4)
v_gr_form          <- c('grass')
# Customized delimiter for `read_delim` function, comma is predefined
v_custom_delimiter <- c(',')


# Main pipelines ---------------------------------------------------------------
source('pipeline/plant_tracker_01.R')

# Select the x_th species (target species)
head(sp_list, 10)
target_spec <- sp_list %>% .[c(4),]  

# Second part of the plant tracker pipeline ------------------------------------

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2024.10.29

# Code adapted from: https://github.com/aestears/plantTracker
# Adapted from plantTracker How to (Stears et al. 2022)


# Specifications ----------------------------------------------------------------
# Define the species variable and abbreviation
v_species <- target_spec[1,1] 
v_sp_abb  <- v_species %>% 
  strsplit(' ') %>% 
  pluck(1) %>% 
  substr( 1, 2 ) %>% 
  paste(collapse = '') %>% 
  str_replace_all(' ', '') %>% 
  tolower

# Directory --------------------------------------------------------------------
# Normalize v_gr_form
v_gr_form_lower <- tolower(v_gr_form)
# Determine the folder suffix
if (exists('v_model_spec')) {
  if (v_model_spec == 'mpm') {
    v_folder_suffix <- paste0(v_sp_abb, '_mpm')
  } else if (v_model_spec == 'ipm') {
    v_folder_suffix <- v_sp_abb
  } else {
    stop("Unknown v_model_spec value: must be 'ipm' or 'mpm'")
  }
} else if (v_gr_form_lower %in% c('forb', 'shrub', 'density')) {
  v_folder_suffix <- paste0(v_sp_abb, '_mpm')
} else {
  v_folder_suffix <- v_sp_abb
}

# Construct paths
dir_R      <- file.path(dir_pub, 'R', v_folder_suffix)
dir_data   <- file.path(dir_pub, 'data', v_folder_suffix)
dir_result <- file.path(dir_pub, 'results', v_folder_suffix)


# Data -------------------------------------------------------------------------
# Read anf format quadrat inventory file 
# Identify quadrat inventory file 
quad_file <- list.files(
  path = dir_qud, full.names = TRUE,
  pattern = 'quad_inventory|quadrat_inventory'
)[1]

if (grepl("\\.RData$", quad_file)) {
  quad_inv <- readRDS(quad_file)
} else {
  quad_inv <- read_delim(
    quad_file, delim = v_delimiter, escape_double = FALSE, trim_ws = TRUE) %>%
    as.data.frame %>% 
    as.list()
}

# Remove NAs
inv        <- lapply(
  X = quad_inv, FUN = function(x) x[is.na(x) == FALSE])

# Replace dots in names with dashes
names(inv) <- gsub('\\.','-',names(inv))  

# Read spatial data (polygon for each species per quadrat)
data <- {
  # Check for files ending with '_quadrats_filtered.rds' in the specified directory
  rds_files <- list.files( 
    dir_qud, pattern = '_quadrats_filtered\\.rds$', full.names = TRUE)
  # If files matching the pattern are found, read the first one
  if (length(rds_files) > 0) {
    readRDS(file = rds_files[1])
  } else {
    # If not fallback to the default file
    readRDS(file = list.files(
      dir_qud, pattern = 'all_filtered\\.rds$', full.names = TRUE)[1])
  }
} %>%
  # Remove the 'type' column
  select(-any_of(c('type'))) %>% 
  clean_names() #%>%   
# This is to accommodate older versions of 'quadrat' scripts
#set_names(c('species', 'site', 'quad', 'year', 'geometry'))


# Subset data for the target species
dat_target_spec <- data %>%
  subset(species %in% target_spec$species) %>%
  filter(
    if (exists("mod_plot") && length(mod_plot) > 0) !(quad %in% mod_plot)
    else TRUE
  ) %>%
  select(
    species, site, 
    matches("quad"),  # selects quad or quadrat
    year, geometry,
    everything()      # all remaining columns
  ) %>% 
  rename(Species = species,
         Site    = site,
         Quad    = matches("quad"),
         Year    = year)

# set the buffer to 0.05 or 5 depending on the unit of measure used in GIS file
buff <- v_buff <- if_else(
  st_bbox(dat_target_spec)[3] < 1.1 | st_bbox(dat_target_spec)[3] > 100, 
  0.05, 5)

# Forbs are not "clonal"; we assume non-forbs are clonal.
v_clonal  <- if_else(v_gr_form == 'forb', F, T)

# # Prepare data for the trackSpp function
# datTrackSpp <- trackSpp(
#   dat = dat_target_spec,
#   inv = inv,
#   # Dormancy flag
#   dorm         = 1,
#   # Buffer size
#   buff         = v_buff,
#   # Allow for clonal tracking
#   clonal       = v_clonal,
#   # Buffer for genet
#   buffGenet    = v_buff,
#   # Aggregate by genet
#   aggByGenet   = TRUE,
#   # Flag potential issues
#   flagSuspects = TRUE)   
# 
# # store (takes too long!)
# saveRDS( datTrackSpp, 
#          paste0(dir_data, '/', gsub('20','',v_author_year), 
#                                     v_region_abb,
#                 '_',v_sp_abb,'_raw.rds') )

# track species
datTrackSpp <- readRDS( paste0(dir_data, '/', 
                               gsub('20','',v_author_year),
                               v_region_abb,
                               '_',v_sp_abb,'_raw.rds') )

# CHECK -- Adaptions -----------------------------------------------------------

# Removal of certain years if unspecified nothing is removed
v_years_re       <- c()
# Define size threshold
v_size_threshold <- c(-10.7)
# Set a complexity to the growth and survival model 
# (NULL = highest AIC / 0 = intercept / 1 = linear / 2 = quadratic / 3 = cubic)
v_mod_set_gr     <- c()
v_mod_set_su     <- c()


# explore the data -------------------------------------------------------------

# load packages
source('helper_functions/load_packages.R')
load_packages( tidyverse, patchwork, skimr, ipmr, binom, bbmle, janitor)


# Specification ----------------------------------------------------------------


v_script_prefix <- str_c(
  str_extract(v_author_year, '^[^_]+'),
  str_sub(str_extract(v_author_year, '_\\d+$'), -2, -1))
# Define prefix for two of the same author and year
if (
  length(
    list.dirs(
      full.names = TRUE, recursive = FALSE)[grepl(
        paste0('^', v_author_year), basename(
          list.dirs(full.names = TRUE, recursive = FALSE)))]
  ) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

# Define suffix for plot outputs
v_suffix     <- ''
if (length(v_years_re)       > 0) {v_suffix <- paste0(v_suffix, '_yr1')}
if (length(v_size_threshold) > 0) {v_suffix <- paste0(v_suffix, '_st1')}
if (length(v_mod_set_gr)     > 0) {v_suffix <- paste0(v_suffix, '_gr1')}
if (length(v_mod_set_su)     > 0) {v_suffix <- paste0(v_suffix, '_su1')}

# Define graph subtitle
v_ggp_suffix <- paste(
  paste0(toupper(substr(v_script_prefix, 1, 1)), 
         substr(v_script_prefix, 2, nchar(v_script_prefix))), '/', 
  v_species, 
  '\n Size threshold:', 
  ifelse(is.null(v_size_threshold), !is.null(v_size_threshold), v_size_threshold),
  '\n Model complexity altered in growth / survival:', 
  ifelse(is.null(v_mod_set_gr), !is.null(v_mod_set_gr), v_mod_set_gr), '/', 
  ifelse(is.null(v_mod_set_su), !is.null(v_mod_set_su), v_mod_set_su),
  '\n Years removed:',
  ifelse(is.null(v_years_re), !is.null(v_years_re), paste(v_years_re, collapse = ', ')))


# Directory --------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

# create R, data, and results directory if they do not exist 
if (!dir.exists(paste0(dir_pub, '/R'))) {
  dir.create(paste0(dir_pub, '/R'))}
if (!dir.exists(paste0(dir_pub, '/data'))) {
  dir.create(paste0(dir_pub, '/data'))}
if (!dir.exists(paste0(dir_pub, '/results'))) {
  dir.create(paste0(dir_pub, '/results'))}

# create species-specific directories if they do not exist
if (!dir.exists(dir_R     )) {dir.create(dir_R     )}
if (!dir.exists(dir_data  )) {dir.create(dir_data  )}
if (!dir.exists(dir_result)) {dir.create(dir_result)}


# Save the suffix --------------------------------------------------------------
write.csv(v_suffix, file.path(dir_data, 'v_suffix.csv'), row.names = F)


# Plant tracker if its not already exists --------------------------------------
if(!file.exists(paste0(dir_data, '/',v_script_prefix, '_',v_sp_abb, '.csv'))) {
  source(paste0(dir_R, '/',v_script_prefix, '_',v_sp_abb, '_tracker.R'))
}


# Data -------------------------------------------------------------------------
# Read and clean the species data
df <- read.csv(paste0(dir_data, '/',v_script_prefix, '_',v_sp_abb, '.csv')) %>% 
  clean_names() %>% 
  filter(species == v_species) %>%
  select(-c(suspect, near_edge, site)) %>%
  mutate(across(c(quad), as.factor)) %>%
  rename(size_t0  = basal_area_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1,
         track_id = track_id) %>%
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),
         logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3)

# Implement size threshold
if (length(v_size_threshold) > 0) {
  df <- df %>% 
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0)) %>% 
    filter(logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

# Growth data frame
grow_df <- df %>% 
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>% 
  select(quad, track_id, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3)

# Growth analysis
g_gr_overall <- grow_df %>% 
  mutate( id_quad_year = paste(track_id,quad,year,
                               sep="; ") ) %>% 
  ggplot( aes( x = logsize_t0, 
               y = logsize_t1,
               group = id_quad_year) ) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') 

# check the plot
g_gr_overall

# Run
ggplotly(g_gr_overall)


# SHINYAPP ---------------------------------------------------------------------

# First, plot using plotly (shiny below doesn't seem to work otherwise)
ggplotly(g_gr_overall)

# Create the "page" for visualization
ui <- fluidPage(
  h3("Hover a point and look below:"),
  plotlyOutput("plt"),
  verbatimTextOutput("row")
)

# set up the way to visualize the plot
server <- function(input, output, session){

  # build the plot and tell Plotly we'll listen for hover
  output$plt <- renderPlotly({
    (ggplotly(g_gr_overall,tooltip = "text")) |>
      event_register("plotly_hover")
  })

  # whenever you hover, grab the row and copy it
  row_r <- eventReactive(event_data("plotly_hover"), {
    h <- event_data("plotly_hover")
    req(h)                          # ignore NULL on app start
    out <- grow_df[h$pointNumber + 1,
                   c('quad','year','track_id') ]  # row in the original data
    clipr::write_clip(out)                 # copies to system clipboard
    out
  })

  output$row <- renderPrint(row_r())
}

# finally launch shiny (select a dot, and copy-paste it typing Ctrl+C)
shinyApp(ui, server)


# Final tests ------------------------------------------------------------------ 

# set bounding box ()
bbox <- st_bbox(datTrackSpp)

# Plot the polygons
plot_polygons <- function( df ){
  ggplot( df ) +
    geom_sf( fill = "lightblue") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"])) +
    theme_bw() 
    
}

# plot individual 
#   before (single indiv in year t0)
#   after (all individuals in a single plot. PlantTracker doesn't allow
#   us to track the "geometry" as well
plot_before_after <- function( track_id, quad, yr ){
  
  # TARGET INDIVIDUAL IN YEAR T0!
  year_t0 <- datTrackSpp %>% 
    subset( trackID == track_id ) %>% 
    subset( Quad == quad ) %>% 
    subset( Year == yr ) %>% 
    select( trackID, geometry ) %>% 
    plot_polygons() +
    labs( title = 'Year t0' )
    
  # This is all of the individuals in the plot the next year. 
  # Unfortunately plant tracker doesn't give us the plot geometries in year t1
  feature_year_t1 <- dat_target_spec %>% 
    subset( Quad == quad ) %>% 
    subset( Year == yr+1 ) %>% 
    select( geometry ) %>% 
    plot_polygons +
    labs( title = 'Year t1' )
    
  # compare before and after    
  grid.arrange( year_t0, 
                feature_year_t1,
                nrow = 1)

}


# identify "suspect" individuals 
ggplotly(g_gr_overall)

# IMPOTANT 	Copy paste here:
  # ID NUMBER
  # QUAD NUMBER
  # YEAR (without the "19" in front of it)
  
# Use these in the function to check before and after
#   IMPORTANT: useyear need be without "19" in front of it!!!
plot_before_after('HILBEL_29_1', 'C1P', 29)



# TO DO: show the transition directly from clicking the figure -----------------

# Exapmle code that works:

# Create the "page" for visualization
ui <- fluidPage(
  h3("Click a point to see the full profile"),
  fluidRow(
    column(6, plotlyOutput("scatter")),
    column(6, uiOutput("detail_ui"))  # placeholder that appears only after a click
  )
)

# set up the way to visualize the plot
server <- function(input, output, session) {
  
  #---- 1. base scatter plot ---------------------------------------
  output$scatter <- renderPlotly({
    # Build a ggplot object first (so you can facet, add aesthetics, etc.)
    base_plot <- ggplot(df, aes(wt, mpg, text = paste("model:", model))) +
      geom_point(size = 3)
    
    # Convert to plotly; assign a source ID so we know **which** plot
    # the click came from (helpful if you have multiple plots).
    (base_plot) |>
      ggplotly(tooltip = "text", source = "scatterSrc") |>
      layout(dragmode = "select")  # enabling box select also works
  })
  
  #---- 2. listen for a click --------------------------------------
  observeEvent(event_data("plotly_click", source = "scatterSrc"), {
    click <- event_data("plotly_click", source = "scatterSrc")
    req(click)                       # guard against NULL
    
    # Look up the row corresponding to pointNumber
    row <- df[click$pointNumber + 1, ]
    
    # Reshape that row to long format for a bar chart
    bar_data <- row %>%
      select(-model) %>%             # drop non‑numeric if needed
      pivot_longer(everything(),
                   names_to  = "variable",
                   values_to = "value")
    
    #---- 3. render the bar chart ----------------------------------
    output$detail <- renderPlotly({
      plot_ly(bar_data,
              x = ~variable, y = ~value,
              type = "bar") |>
        layout(title = paste0("All variables for ", row$model),
               xaxis = list(title = "Variable"),
               yaxis = list(title = "Value"))
    })
    
    # Insert (or replace) the UI slot with the new plotlyOutput
    output$detail_ui <- renderUI({
      plotlyOutput("detail")
    })
    
  }, ignoreNULL = TRUE)
  
}
