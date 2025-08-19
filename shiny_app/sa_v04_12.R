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

# ---- POLYGON + BUFFER + OTHERS FUNCTION ----
plot_polygons_with_others <- function(df_selected, quad_val, year_val, show_buffer = TRUE, show_others = TRUE, src) {
  
  if (nrow(df_selected) == 0 || !"geometry" %in% names(df_selected)) {
    p <- ggplot() +
      annotate("text", x = mean(bbox[c("xmin", "xmax")]),
               y = mean(bbox[c("ymin", "ymax")]),
               label = "No data available") +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
               ylim = c(bbox["ymin"], bbox["ymax"])) +
      theme_bw()
    return(ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout"))
  }
  
  df_selected <- st_cast(df_selected, "MULTIPOLYGON")
  
  # Buffer distance for selected individual
  if (show_buffer) {
    buff <- if_else(st_bbox(df_selected)[3] < 1.1 | st_bbox(df_selected)[3] > 100, 0.05, 5)
    df_buffer <- st_buffer(df_selected, dist = buff)
  }
  
  # Other polygons for the same plot & year
  df_others <- df %>%
    filter(quad == quad_val, year == year_val, !(track_id %in% df_selected$track_id)) %>%
    st_cast("MULTIPOLYGON")
  
  p <- ggplot()
  
  if (show_others && nrow(df_others) > 0) {
    p <- p + geom_sf(data = df_others, fill = "lightgrey", alpha = 0.3)
  }
  
  if (show_buffer) p <- p + geom_sf(data = df_buffer, fill = "lightblue")
  p <- p + geom_sf(data = df_selected, fill = "dodgerblue") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"])) +
    theme_bw()
  
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

# Function to generate t0 and t1 polygon maps
plot_before_after_plotly <- function(track_id_val, quad_val, yr_val, show_buffer, show_others) {
  
  year_t0 <- df %>%
    filter(track_id == track_id_val, quad == quad_val, year == yr_val) %>%
    select(track_id, geometry)
  
  year_t1 <- df %>%
    filter(track_id == track_id_val, quad == quad_val, year == yr_val + 1) %>%
    select(track_id, geometry)
  
  list(
    t0 = plot_polygons_with_others(year_t0, quad_val, yr_val, show_buffer, show_others, src = "t0"),
    t1 = plot_polygons_with_others(year_t1, quad_val, yr_val + 1, show_buffer, show_others, src = "t1")
  )
}

# ------------------ UI ------------------
ui <- fluidPage(
  h3("Interactive Growth & Polygon Viewer"),
  
  fluidRow(
    # Left column: clicked point summary
    column(3,
           verbatimTextOutput("row")
    ),
    
    # Right column: growth plot, polygons, toggles
    column(9,
           plotlyOutput("plt", height = "300px"),
           fluidRow(
             column(6, plotlyOutput("map_t0", height = "300px")),
             column(6, plotlyOutput("map_t1", height = "300px"))
           ),
           fluidRow(
             column(6, checkboxInput("show_buffer", "Show buffer", value = TRUE)),
             column(6, checkboxInput("show_others", "Show other individuals", value = TRUE))
           )
    )
  )
)

# ---------------- SERVER ----------------
server <- function(input, output, session) {
  
  yr_val <- reactiveVal(NULL)
  
  # Scatter plot
  output$plt <- renderPlotly({
    g_gr_overall <- grow_df %>%
      mutate(tooltip = paste0(
        "id_quad_year: ", id_quad_year,
        "<br>track_id: ", as.character(track_id),
        "<br>logsize_t0: ", round(as.numeric(logsize_t0), 3),
        "<br>logsize_t1: ", round(as.numeric(logsize_t1), 3)
      )) %>%
      ggplot(aes(x = logsize_t0, y = logsize_t1, text = tooltip, group = id_quad_year)) +
      geom_point(alpha = 0.7, pch = 16, size = 0.7, color = "#D55E00") +
      theme_bw(base_size = 11) +
      labs(x = "Log(Size t0)", y = "Log(Size t1)") +
      theme(panel.grid.major = element_line(color = "grey85"),
            panel.grid.minor = element_line(color = "grey92"),
            legend.position = "none")
    
    ggplotly(g_gr_overall, tooltip = "text") |> event_register("plotly_click")
  })
  
  # Clicked row
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df[click$pointNumber + 1, c("track_id", "quad", "year")]
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint(row_r())
  
  # Render t0/t1 maps
  observeEvent(list(row_r(), input$show_buffer, input$show_others), {
    r <- row_r()
    req(r)
    yr_val(as.numeric(r$year))
    plots <- plot_before_after_plotly(r$track_id, r$quad, as.numeric(r$year), input$show_buffer, input$show_others)
    
    output$map_t0 <- renderPlotly(plots$t0)
    output$map_t1 <- renderPlotly(plots$t1)
  })
  
  # Sync zoom
  observeEvent(event_data("plotly_relayout", source = "t0"), {
    rel <- event_data("plotly_relayout", source = "t0")
    if (!is.null(rel)) {
      plotlyProxy("map_t1", session) %>% plotlyProxyInvoke("relayout", rel)
    }
  })
  
  observeEvent(event_data("plotly_relayout", source = "t1"), {
    rel <- event_data("plotly_relayout", source = "t1")
    if (!is.null(rel)) {
      plotlyProxy("map_t0", session) %>% plotlyProxyInvoke("relayout", rel)
    }
  })
}

shinyApp(ui, server)
