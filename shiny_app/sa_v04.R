# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.13

# Version 04: introduce linked zoom

# ---- PACKAGES ----
source('helper_functions/load_packages.R')
load_packages(tidyverse, sf, shiny, plotly, gridExtra, clipr)

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

dir.create(dir_R, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_result, recursive = TRUE, showWarnings = FALSE)

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

# Bounding box for polygons
bbox <- st_bbox(df)

# Function to make an interactive polygon map with fixed bbox
plot_polygons_plotly <- function(df, title, src) {
  # Handle empty data
  if (nrow(df) == 0 || !"geometry" %in% names(df)) {
    p <- ggplot() +
      annotate("text", x = mean(bbox[c("xmin", "xmax")]),
               y = mean(bbox[c("ymin", "ymax")]),
               label = "No data available") +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"])
      ) +
      labs(title = title) +
      theme_bw()
    return(ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout"))
  }
  
  p <- ggplot(df) +
    geom_sf(fill = "lightblue") +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"])
    ) +
    labs(title = title) +
    theme_bw()
  ggplotly(p, dynamicTicks = TRUE, source = src) |>
    event_register("plotly_relayout")
}

# Function to generate t0 and t1 polygon maps
plot_before_after_plotly <- function(track_id_val, quad_val, yr_val) {
  
  # t0 polygons
  year_t0 <- df %>%
    filter(
      .data[["track_id"]] == track_id_val,
      .data[["quad"]] == quad_val,
      .data[["year"]] == yr_val
    ) %>%
    select(track_id, geometry)
  
  # t1 polygons
  feature_year_t1 <- df %>%
    filter(
      .data[["track_id"]] == track_id_val,
      .data[["quad"]] == quad_val,
      .data[["year"]] == yr_val + 1
    ) %>%
    select(track_id, geometry)
  
  list(
    t0 = plot_polygons_plotly(year_t0, paste("year", yr_val), src = "t0"),
    t1 = plot_polygons_plotly(feature_year_t1, paste("year", yr_val + 1), src = "t1")
  )
}


# ------------------ UI ------------------
ui <- fluidPage(
  h3("Click a point to see polygons for year t0 and year t1:"),
  
  fluidRow(
    column(8, plotlyOutput("plt")),
    column(4, verbatimTextOutput("row"))
  ),
  
  fluidRow(
    column(6, plotlyOutput("map_t0", height = "400px")),
    column(6, plotlyOutput("map_t1", height = "400px"))
  )
)

# ---------------- SERVER ----------------
server <- function(input, output, session) {
  
  yr_val <- reactiveVal(NULL)  # store yr_num globally for syncing
  
  # Main scatter plot
  output$plt <- renderPlotly({
    ggplotly(g_gr_overall, tooltip = "text") |>
      event_register("plotly_click")
  })
  
  # Get clicked row
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df[click$pointNumber + 1, c("track_id", "quad", "year")]
    clipr::write_clip(out)  # optional clipboard copy
    out
  })
  
  # Show metadata
  output$row <- renderPrint(row_r())
  
  # Render t0 and t1 maps when clicked
  observeEvent(row_r(), {
    r <- row_r()
    req(r)
    yr_num <- as.numeric(r$year)
    yr_val(yr_num)  # store for sync
    
    plots <- plot_before_after_plotly(r$track_id, r$quad, yr_num)
    
    output$map_t0 <- renderPlotly(plots$t0)
    output$map_t1 <- renderPlotly(plots$t1)
  })
  
  # SYNC ZOOM between t0 and t1
  observeEvent(event_data("plotly_relayout", source = "t0"), {
    rel <- event_data("plotly_relayout", source = "t0")
    if (!is.null(rel)) {
      plotlyProxy("map_t1", session) %>%
        plotlyProxyInvoke("relayout", rel)
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(event_data("plotly_relayout", source = "t1"), {
    rel <- event_data("plotly_relayout", source = "t1")
    if (!is.null(rel)) {
      plotlyProxy("map_t0", session) %>%
        plotlyProxyInvoke("relayout", rel)
    }
  }, ignoreNULL = TRUE)
}

shinyApp(ui, server)
