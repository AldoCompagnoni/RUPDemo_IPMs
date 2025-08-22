# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.22

# ---- PACKAGES ----
source('helper_functions/load_packages.R')
load_packages(tidyverse, sf, shiny, plotly, gridExtra, clipr)

# ---- SPECIFICATION ----
v_author_year <- 'anderson_2016'
v_region_abb <- 'az'
v_species <- 'Hilaria belangeri'
v_size_threshold <- -10.7

# Species abbreviation
v_sp_abb <- tolower(gsub(' ', '', paste(substr(unlist(strsplit(v_species, ' ')), 1, 2), collapse = '')))

# Script prefix
v_script_prefix <- str_c(
  str_extract(v_author_year, '^[^_]+'),
  str_sub(str_extract(v_author_year, '_\\d+$'), -2, -1)
)
if (length(list.dirs(full.names = TRUE, recursive = FALSE)[
  grepl(paste0('^', v_author_year), basename(list.dirs(full.names = TRUE, recursive = FALSE)))
]) > 1) {
  v_script_prefix <- paste0(v_script_prefix, v_region_abb)
}

# ---- DIRECTORIES ----
dir_pub <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R <- file.path(dir_pub, 'R', v_sp_abb)
dir_data <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

dir.create(dir_R, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_result, recursive = TRUE, showWarnings = FALSE)

# ---- DATA ----
df <- readRDS(file.path(dir_data, paste0(v_script_prefix, '_', v_sp_abb, '_raw.rds'))) %>%
  select(-c(suspect, near_edge)) %>%
  mutate(across(c(quad, track_id), as.factor)) %>%
  rename(size_t0 = basal_area_genet,
         size_t1 = size_tplus1,
         survives = survives_tplus1) %>%
  group_by(track_id, quad) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(next_year = lead(year),
         size_t1 = if_else(!is.na(next_year) & next_year == year + 1, size_t1, NA_real_)) %>%
  ungroup() %>%
  mutate(logsize_t0 = log(size_t0),
         logsize_t1 = log(size_t1),
         logsize_t0_2 = logsize_t0^2,
         logsize_t0_3 = logsize_t0^3,
         site = sub("^.{4}.*?_(.*?)_.*$", "\\1", as.character(track_id)),
         track_id_quad = paste0(track_id, "_", quad),
         id_quad_year = paste0(track_id, "_", quad, "_", year))

if (!is.null(v_size_threshold)) {
  df <- df %>%
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0),
           logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

grow_df <- df %>%
  subset(size_t0 != 0 & size_t1 != 0) %>%
  select(site, quad, track_id, track_id_quad, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry, id_quad_year)

bbox <- st_bbox(df)

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
      verbatimTextOutput("row")
    ),
    column(
      width = 8,
      plotlyOutput("plt", height = "500px")
    )
  )
)

server <- function(input, output, session) {
  
  output$plt <- renderPlotly({
    g <- grow_df %>%
      mutate(
        track_id_chr = as.character(track_id),
        highlight = case_when(
          is.null(input$track_query) || trimws(input$track_query) == "" ~ "all",
          track_id_chr == input$track_query ~ "selected",
          TRUE ~ "other"
        ),
        tooltip = paste0(
          "track_id: ", track_id_chr,
          "<br>logsize_t0: ", round(logsize_t0, 3),
          "<br>logsize_t1: ", round(logsize_t1, 3),
          "<br>year: ", year
        )
      ) %>%
      ggplot(aes(x = logsize_t0, y = logsize_t1, text = tooltip)) +
      geom_point(aes(alpha = highlight, color = highlight), size = 1.2) +
      scale_alpha_manual(values = c("selected" = 1, "other" = 0.15, "all" = 0.7)) +
      scale_color_manual(values = c("selected" = "#D55E00", "other" = "grey60", "all" = "#D55E00")) +
      theme_bw(base_size = 12) +
      labs(x = "Log(Size t0)", y = "Log(Size t1)") +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_line(color = "grey92")
      )
    
    ggplotly(g, tooltip = "text") |> event_register("plotly_click")
  })
  
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df[click$pointNumber + 1,
                   c("track_id", "quad", "year", "logsize_t0", "logsize_t1")]
    out$track_id  <- as.character(out$track_id)
    out$quad      <- as.character(out$quad)
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
      "log(size_t1):", round(r$logsize_t1, 3), "\n"
    )
  })
}

shinyApp(ui, server)
