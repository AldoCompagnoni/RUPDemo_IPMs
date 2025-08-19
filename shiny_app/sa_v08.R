# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.14

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
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>%
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
    theme(
      axis.title = element_blank()
    )
}

plot_before_after <- function(track_id_val, quad_val, yr) {
  
  t0_poly <- df %>% 
    filter(track_id == track_id_val, quad == quad_val, year == yr)
  
  t1_poly <- df %>%
    filter(track_id == track_id_val, quad == quad_val, year == yr + 1)
  
  # Year t0
  year_t0 <- plot_polygons(t0_poly) + 
    labs(title = paste("Year", yr))
  
  # Year t1
  if (nrow(t1_poly) > 0) {
    year_t1 <- plot_polygons(t1_poly) +
      labs(title = paste("Year", yr + 1))
  } else {
    year_t1 <- ggplot() +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"])
      ) +
      theme_bw() +
      theme_void() +
      labs(title = paste("Year", yr + 1)) +
      annotate(
        "text",
        x = mean(c(bbox["xmin"], bbox["xmax"])),
        y = mean(c(bbox["ymin"], bbox["ymax"])),
        label = "Polygon missing",
        size = 4,
        color = "darkred",  
        hjust = 0.5, vjust = 0.5
      )
  }
  
  grid.arrange(year_t0, year_t1, nrow = 1)
}

# ---- SHINY APP ----
ui <- fluidPage(
  fluidRow(
    column(
      width = 4,
      selectizeInput(
        "track_query",
        "Select or type track_id:",
        choices = c("", sort(unique(as.character(grow_df$track_id)))),
        selected = "",
        multiple = FALSE,
        options = list(
          placeholder = 'Type to search…',
          create = FALSE
        )
      ),
      verbatimTextOutput("row"),
      hr(),
      selectizeInput(
        "exclude_input",
        "Exclude id_quad_year:",
        choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
        selected = "",
        multiple = FALSE,
        options = list(
          placeholder = 'Type to exclude…',
          create = FALSE
        )
      ),
      actionButton("add_exclude", "Add to exclude list", class = "btn btn-sm btn-danger"),
      br(), br(),
      
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          tags$h5(
            class = "panel-title",
            tags$a(
              href = "#collapseExcluded",
              "data-toggle" = "collapse",
              tags$small("Excluded id_quad_years")
            )
          )
        ),
        tags$div(
          id = "collapseExcluded",
          class = "panel-collapse collapse",
          tags$div(
            class = "panel-body",
            uiOutput("excluded_list")
          )
        )
      ),
      br(),
      downloadButton("download_excluded", "Download excluded CSV", class = "btn btn-sm btn-primary")
    ),
    
    column(
      width = 8,
      plotlyOutput("plt", height = "400px"),
      plotOutput("before_after_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  excluded_ids <- reactiveVal(character())
  
  observeEvent(input$add_exclude, {
    req(input$exclude_input)
    excluded_ids(unique(c(excluded_ids(), input$exclude_input)))
  })
  
  output$excluded_list <- renderUI({
    ids <- excluded_ids()
    if (length(ids) == 0) return("None")
    tagList(
      lapply(ids, function(id) {
        fluidRow(
          column(8, id),
          column(4, actionButton(
            paste0("remove_", id),
            label = "x",
            style = "padding: 2px 6px; font-size: 10px; height: 22px; line-height: 1;"
          ))
        )
      })
    )
  })
  
  observe({
    ids <- excluded_ids()
    lapply(ids, function(id) {
      observeEvent(input[[paste0("remove_", id)]], {
        excluded_ids(setdiff(excluded_ids(), id))
      }, ignoreInit = TRUE)
    })
  })
  
  grow_df_filtered <- reactive({
    q <- input$track_query
    dat <- grow_df
    if (!is.null(q) && trimws(q) != "") {
      dat <- dat %>%
        mutate(track_id_chr = as.character(track_id)) %>%
        filter(str_detect(track_id_chr, fixed(q, ignore_case = TRUE)))
    }
    dat <- dat %>% filter(!id_quad_year %in% excluded_ids())
    dat
  })
  
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat) > 0, "No rows match the current filter."))
    
    g <- dat %>%
      mutate(
        tooltip = paste0(
          "id_quad_year: ", id_quad_year,
          "<br>track_id: ", as.character(track_id),
          "<br>logsize_t0: ", round(as.numeric(logsize_t0), 3),
          "<br>logsize_t1: ", round(as.numeric(logsize_t1), 3)
        )
      ) %>%
      ggplot(aes(
        x = logsize_t0,
        y = logsize_t1,
        text = tooltip,
        group = id_quad_year
      )) +
      geom_point(
        alpha = 0.7,
        pch = 16,
        size = 0.7,
        color = "#D55E00"
      ) +
      theme_bw(base_size = 11) +
      labs(
        x = "Log(Size t0)",
        y = "Log(Size t1)"
      ) +
      theme(
        panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_line(color = "grey92"),
        legend.position = "none"
      )
    
    ggplotly(g, tooltip = "text") |> event_register("plotly_click")
  })
  
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df_filtered()[click$pointNumber + 1,
                              c("track_id", "quad", "year", "logsize_t0", "logsize_t1", "id_quad_year")]
    out$track_id   <- as.character(out$track_id)
    out$quad       <- as.character(out$quad)
    out$logsize_t0 <- as.numeric(out$logsize_t0)
    out$logsize_t1 <- as.numeric(out$logsize_t1)
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint({
    r <- row_r()
    req(r)
    cat(
      "track_id:", r$track_id, "\n",
      "quad:", r$quad, "\n",
      "year:", r$year, "\n",
      "log(size_t0):", round(r$logsize_t0, 3), "\n",
      "log(size_t1):", round(r$logsize_t1, 3), "\n",
      "id_quad_year:", r$id_quad_year, "\n"
    )
  })
  
  output$before_after_plot <- renderPlot({
    r <- row_r()
    req(r)
    plot_before_after(
      track_id_val = as.character(r$track_id),
      quad_val     = as.character(r$quad),
      yr           = as.numeric(r$year)
    )
  })
  
  output$download_excluded <- downloadHandler(
    filename = function() {
      paste0("excluded_id_quad_years_", Sys.Date(), ".csv")
    },
    content = function(file) {
      ids <- excluded_ids()
      write_csv(tibble(id_quad_year = ids), file)
    }
  )
}

shinyApp(ui, server)
