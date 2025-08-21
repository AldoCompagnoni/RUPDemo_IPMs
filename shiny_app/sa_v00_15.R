# Shiny app data - Growth Analysis Version 2 (Only Scatter + Info)
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.18

# Version 00: Only growth scatterplot + clicked point details

# ---- PACKAGES ----
source('helper_functions/load_packages.R')
load_packages(tidyverse, sf, shiny, plotly, clipr)

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
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3) %>%
  mutate(id_quad_year = paste0(track_id, "_", quad, "_", year))

# ---- Growth scatter ----
g_gr_overall <- grow_df %>%
  ggplot(aes(x = logsize_t0, y = logsize_t1, group = id_quad_year)) +
  geom_point(alpha = 0.5, pch = 16, size = 0.7, color = 'red') +
  theme_bw()

# ---- SHINY APP ----
ui <- fluidPage(
  h3("Click a point to see details:"),
  sliderInput("jit_w", "Jitter width (x-axis):",
              min = 0, max = 0.5, value = 0, step = 0.01),
  sliderInput("jit_h", "Jitter height (y-axis):",
              min = 0, max = 0.5, value = 0, step = 0.01),
  plotlyOutput("plt"),
  verbatimTextOutput("row")
)

server <- function(input, output, session) {
  # Main scatter plot with jitter
  output$plt <- renderPlotly({
    set.seed(123)  # ensures jitter is always the same
    g <- grow_df %>%
      ggplot(aes(x = logsize_t0, y = logsize_t1, group = id_quad_year)) +
      geom_jitter(alpha = 0.5, pch = 16, size = 0.7, color = 'red',
                  width = input$jit_w, height = input$jit_h) +
      theme_bw()
    
    ggplotly(g, tooltip = "text") |>
      event_register("plotly_click")
  })
  
  # Capture clicked row
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    out <- grow_df[click$pointNumber + 1, c("track_id", "quad", "year", "size_t0", "size_t1")]
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint(row_r())
}

shinyApp(ui, server)
