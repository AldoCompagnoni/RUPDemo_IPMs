# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.18

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
  select(-c(suspect, near_edge)) %>%
  mutate(
    across(c(quad, track_id), as.factor),
    site = sub("^.{3}_([A-Za-z]+)_.*$", "\\1", track_id) # extract site from track_id
  ) %>%
  rename(size_t0  = basal_area_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1) %>%
  group_by(track_id, quad) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    # check the next row's year
    next_year = lead(year),
    # keep size_t1 only if it's consecutive year
    size_t1   = if_else(!is.na(next_year) & next_year == year + 1, size_t1, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    logsize_t0   = log(size_t0),
    logsize_t1   = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    track_id_quad = paste0(track_id, "_", quad)
  )

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
  select(site, quad, track_id, track_id_quad, year, size_t0, survives, size_t1, 
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry) %>%
  mutate(id_quad_year = paste0(track_id, "_", quad, "_", year))

# ---- PLOTS ----
bbox <- st_bbox(df)

plot_polygons <- function(df) {
  ggplot(df) +
    geom_sf(fill = "lightblue") +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"])
    ) +
    theme_bw() +
    theme(axis.title = element_blank())
}

plot_before_after <- function(track_id_val, quad_val, yr) {
  t0_poly <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr)
  t1_poly <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr + 1)
  
  year_t0 <- plot_polygons(t0_poly) + labs(title = paste("Year", yr))
  
  year_t1 <- if (nrow(t1_poly) > 0) {
    plot_polygons(t1_poly) + labs(title = paste("Year", yr + 1))
  } else {
    ggplot() +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
               ylim = c(bbox["ymin"], bbox["ymax"])) +
      theme_bw() +
      theme(axis.title = element_blank(), panel.grid = element_blank()) +
      labs(title = paste("Year", yr + 1)) +
      annotate("text",
               x = mean(c(bbox["xmin"], bbox["xmax"])),
               y = mean(c(bbox["ymin"], bbox["ymax"])),
               label = "Polygon missing",
               size = 4, color = "darkred", hjust = 0.5, vjust = 0.5)
  }
  
  grid.arrange(year_t0, year_t1, nrow = 1)
}

# ---- SHINY APP ----
ui <- fluidPage(
  fluidRow(
    column(
      width = 4,
      # ---- Filters panel (foldable) ----
      tags$div(class="panel panel-default",
               tags$div(class="panel-heading",
                        tags$h5(class="panel-title",
                                tags$a(href="#collapseFilters", "data-toggle"="collapse",
                                       tags$strong("Filters"))
                        )
               ),
               tags$div(id="collapseFilters", class="panel-collapse collapse",  # starts closed
                        tags$div(class="panel-body",
                                 selectizeInput("track_query", "Select track_id_quad:",
                                                choices = c("", sort(unique(grow_df$track_id_quad))),
                                                selected = "", multiple = FALSE),
                                 selectizeInput("site_query", "Filter by site:",
                                                choices = c("", sort(unique(grow_df$site))),
                                                selected = "", multiple = FALSE),
                                 selectizeInput("quad_query", "Filter by quad:",
                                                choices = c("", sort(unique(as.character(grow_df$quad)))),
                                                selected = "", multiple = FALSE),
                                 sliderInput("year_query", "Select year(s):",
                                             min = min(grow_df$year, na.rm=TRUE),
                                             max = max(grow_df$year, na.rm=TRUE),
                                             value = c(min(grow_df$year, na.rm=TRUE),
                                                       max(grow_df$year, na.rm=TRUE)),
                                             step = 1, sep = "")
                        )
               )
      ),
      # You can add other sidebar elements below the foldable panel
      verbatimTextOutput("row")
    ),
    column(
      width = 8,
      plotlyOutput("plt", height = "400px"),
      plotOutput("before_after_plot", height = "400px")
    )
  )
)



server <- function(input, output, session) {
  
  grow_df_filtered <- reactive({
    dat <- grow_df
    if (!is.null(input$track_query) && input$track_query != "") {
      dat <- dat %>% filter(track_id_quad == input$track_query)
    }
    if (!is.null(input$site_query) && input$site_query != "") {
      dat <- dat %>% filter(site == input$site_query)
    }
    if (!is.null(input$quad_query) && input$quad_query != "") {
      dat <- dat %>% filter(as.character(quad) == input$quad_query)
    }
    if (!is.null(input$year_query) && length(input$year_query) == 2) {
      dat <- dat %>% filter(year >= input$year_query[1], year <= input$year_query[2])
    }
    dat
  })
  
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat) > 0, "No rows match the current filter."))
    
    g <- dat %>%
      mutate(tooltip = paste0(
        "track_id: ", as.character(track_id),
        "<br>quad: ", as.character(quad),
        "<br>year: ", year,
        "<br>logsize_t0: ", round(logsize_t0, 3),
        "<br>logsize_t1: ", round(logsize_t1, 3)
      )) %>%
      ggplot(aes(x = logsize_t0, y = logsize_t1, text = tooltip, group = id_quad_year)) +
      geom_point(alpha = 0.7, pch = 16, size = 0.7, color = "#D55E00") +
      theme_bw(base_size = 11) +
      labs(x = "Log(Size t0)", y = "Log(Size t1)") +
      theme(panel.grid.major = element_line(color = "grey85"),
            panel.grid.minor = element_line(color = "grey92"),
            legend.position = "none")
    
    ggplotly(g, tooltip = "text") |> event_register("plotly_click")
  })
  
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df_filtered()[click$pointNumber + 1,
                              c("track_id", "quad", "year", "logsize_t0", "logsize_t1")]
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint({
    r <- row_r()
    req(r)
    cat("track_id:", r$track_id, "\n",
        "quad:", r$quad, "\n",
        "year:", r$year, "\n",
        "log(size_t0):", round(r$logsize_t0, 3), "\n",
        "log(size_t1):", round(r$logsize_t1, 3), "\n")
  })
  
  output$before_after_plot <- renderPlot({
    r <- row_r()
    req(r)
    plot_before_after(track_id_val = as.character(r$track_id),
                      quad_val = as.character(r$quad),
                      yr = as.numeric(r$year))
  })
}

shinyApp(ui, server)
