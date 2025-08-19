library(shiny)
library(plotly)
library(sf)
library(ggplot2)
library(clipr)

# ----------------------------------------------------------
# YOUR EXISTING DATA OBJECTS (assumed to exist in environment):
# g_gr_overall   = ggplot object for main scatter
# grow_df        = dataframe with columns: track_id, quad, year
# datTrackSpp    = sf object with polygons for individuals
# dat_target_spec= sf object with polygons for next year
# ----------------------------------------------------------

# Bounding box for polygons
bbox <- st_bbox(datTrackSpp)

# Function to make an interactive polygon map with fixed bbox
plot_polygons_plotly <- function(df, title, src) {
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
plot_before_after_plotly <- function(track_id, quad, yr) {
  
  # t0 polygons — filtered by track_id, quad, and year
  year_t0 <- datTrackSpp %>%
    filter(trackID == track_id, Quad == quad, Year == yr) %>%
    select(trackID, geometry)
  
  # t1 polygons — filtered by same track_id, quad, and year+1
  feature_year_t1 <- datTrackSpp %>%
    filter(trackID == track_id, Quad == quad, Year == yr + 1) %>%
    select(trackID, geometry)
  
  list(
    t0 = plot_polygons_plotly(year_t0, paste("Year", yr), src = "t0"),
    t1 = plot_polygons_plotly(feature_year_t1, paste("Year", yr + 1), src = "t1")
  )
}

# ------------------ UI ------------------
ui <- fluidPage(
  h3("Click a point to see polygons for Year t0 and Year t1:"),
  
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
    if (yr_num > 1900) yr_num <- yr_num - 1900
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
