# Shiny app data - Growth Analysis Version 1
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.18

# Version 1: Simple click-to-show polygons (before/after)

# ---- PACKAGES ----
source('helper_functions/load_packages.R')
load_packages(tidyverse, sf, shiny, plotly, clipr, gridExtra)

# ---- SPECIFICATION ----
v_author_year    <- c('anderson_2016')
v_region_abb     <- c('az')
v_species        <- c('Hilaria belangeri')
v_size_threshold <- c(-10.7)

# Species abbreviation
v_sp_abb <- tolower(
  gsub(' ', '', paste(
    substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Script prefix
v_script_prefix <- str_c(
  str_extract(v_author_year, '^[^_]+'),
  str_sub(str_extract(v_author_year, '_\\d+$'), -2, -1))
if (
  length(
    list.dirs(full.names = TRUE, recursive = FALSE)[
      grepl(paste0('^', v_author_year),
            basename(list.dirs(full.names = TRUE, recursive = FALSE)))
    ]) > 1
) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

# ---- DIRECTORIES ----
dir_pub    <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R      <- file.path(dir_pub, 'R',       v_sp_abb)
dir_data   <- file.path(dir_pub, 'data',    v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

# ---- DATA ----
df <- readRDS(file.path(
  dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_raw.rds'))) %>%
  select(-c(suspect, near_edge, site)) %>%
  mutate(across(c(quad, track_id), as.factor)) %>%
  rename(size_t0  = basal_area_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1) %>%
  mutate(logsize_t0   = log(size_t0),
         logsize_t1   = log(size_t1),
         logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3)

# Size threshold
if (length(v_size_threshold) > 0) {
  df <- df %>%
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0)) %>%
    filter(logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

# Growth dataframe
grow_df <- df %>%
  subset(size_t0 != 0) %>%
  subset(size_t1 != 0) %>%
  select(quad, track_id, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry) %>%
  mutate(id_quad_year = paste0(track_id, "_", quad, "_", year))

# ---- Growth scatter ----
g_gr_overall <- grow_df %>%
  ggplot(aes(x = logsize_t0, y = logsize_t1, group = id_quad_year)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw()

# ---- Bounding box ----
bbox <- st_bbox(df)

# Function to plot polygons with fixed bbox
plot_polygons <- function(df) {
  ggplot(df) +
    geom_sf(fill = "lightblue") +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"])
    ) +
    theme_bw()
}

# Function to generate before/after comparison
plot_before_after <- function(t_id, quad, yr) {
  year_t0 <- df %>%
    filter(track_id == t_id, quad == quad, year == yr) %>%
    select(track_id, geometry) %>%
    plot_polygons() +
    labs(title = paste("Year", yr))
  
  feature_year_t1 <- df %>%
    filter(track_id == t_id, quad == quad, year == yr + 1) %>%
    select(track_id, geometry) %>%
    plot_polygons() +
    labs(title = paste("Year", yr + 1))
  
  grid.arrange(year_t0, feature_year_t1, nrow = 1)
}

# ---- SHINY APP ----
ui <- fluidPage(
  h3("Click a point to see polygons for Year t0 and Year t1:"),
  plotlyOutput("plt"),
  verbatimTextOutput("row"),
  plotOutput("before_after_plot", height = "400px")
)

server <- function(input, output, session) {
  # Main scatter plot
  output$plt <- renderPlotly({
    ggplotly(g_gr_overall, tooltip = "text") |>
      event_register("plotly_click")
  })
  
  # Capture clicked row
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df[click$pointNumber + 1, c("track_id", "quad", "year")]
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint(row_r())
  
  # Polygons before/after
  output$before_after_plot <- renderPlot({
    r <- row_r()
    req(r)
    yr_num <- as.numeric(r$year)
    plot_before_after(t_id = r$track_id, quad = r$quad, yr = yr_num)
  })
}

shinyApp(ui, server)
