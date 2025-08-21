# ---- PACKAGES ----
source('helper_functions/load_packages.R')
load_packages(tidyverse, shiny, plotly, clipr)

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
        options = list(placeholder = 'Type to search…', create = FALSE)
      ),
      verbatimTextOutput("row"),
      hr(),
      selectizeInput(
        "misfit_input",
        "Add id_quad_year to misfit list:",
        choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
        selected = "",
        multiple = FALSE,
        options = list(placeholder = 'Type to add…', create = FALSE)
      ),
      
      # Only button now
      actionButton("add_misfit", "Add to misfit list", class = "btn btn-sm btn-danger"),
      br(), br(),
      
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          tags$h5(class = "panel-title",
                  tags$a(href = "#collapseMisfits", "data-toggle" = "collapse",
                         tags$small("Misfit id_quad_years")))
        ),
        tags$div(
          id = "collapseMisfits", class = "panel-collapse collapse",
          tags$div(class = "panel-body", uiOutput("misfit_list"))
        )
      ),
      br(),
      downloadButton("download_misfits", "Download misfit CSV", class = "btn btn-sm btn-primary"),
      br(), br(),
      fileInput(
        "upload_misfits",
        "Upload misfit CSV:",
        accept = c(".csv"),
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      )
    ),
    
    column(
      width = 8,
      plotlyOutput("plt", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  # Store misfits as tibble with id + status
  misfit_tbl <- reactiveVal(tibble(id_quad_year = character(), status = character()))
  
  # Add misfit with default status = "include"
  observeEvent(input$add_misfit, {
    req(input$misfit_input)
    new_row <- tibble(
      id_quad_year = input$misfit_input,
      status       = "include"
    )
    old_tbl <- misfit_tbl()
    updated <- bind_rows(new_row, old_tbl %>% filter(id_quad_year != input$misfit_input))
    misfit_tbl(updated)
  })
  
  # Misfit list UI with editable status (THE ONLY PLACE TO CHANGE STATUS)
  output$misfit_list <- renderUI({
    df <- misfit_tbl()
    if (nrow(df) == 0) return("None")
    tagList(
      lapply(seq_len(nrow(df)), function(i) {
        id <- df$id_quad_year[i]
        status <- df$status[i]
        fluidRow(
          column(5, id),
          column(3,
                 selectInput(
                   paste0("status_", id),
                   NULL,
                   choices = c("include", "exclude"),
                   selected = status,
                   width = "100%"
                 )
          ),
          column(4, actionButton(
            paste0("remove_", id),
            label = "x",
            style = "padding: 2px 6px; font-size: 10px; height: 22px; line-height:1;"
          ))
        )
      })
    )
  })
  
  # Update status dynamically
  observe({
    df <- misfit_tbl()
    lapply(df$id_quad_year, function(id) {
      observeEvent(input[[paste0("status_", id)]], {
        df2 <- misfit_tbl()
        if (id %in% df2$id_quad_year) {
          df2$status[df2$id_quad_year == id] <- input[[paste0("status_", id)]]
          misfit_tbl(df2)
        }
      }, ignoreInit = TRUE)
    })
  })
  
  # Remove button
  observe({
    df <- misfit_tbl()
    lapply(df$id_quad_year, function(id) {
      observeEvent(input[[paste0("remove_", id)]], {
        df2 <- misfit_tbl()
        misfit_tbl(df2 %>% filter(id_quad_year != id))
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
    dat
  })
  
  # Plot
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat) > 0, "No rows match the current filter."))
    
    df <- misfit_tbl()
    dat <- dat %>%
      mutate(
        status = ifelse(id_quad_year %in% df$id_quad_year,
                        df$status[match(id_quad_year, df$id_quad_year)],
                        NA),
        is_misfit = !is.na(status),
        color = case_when(
          status == "exclude" ~ "grey",
          status == "include" ~ "#1b9e77",
          TRUE ~ "#D55E00"
        ),
        tooltip = paste0(
          "id_quad_year: ", id_quad_year,
          "<br>track_id: ", as.character(track_id),
          "<br>logsize_t0: ", round(as.numeric(logsize_t0), 3),
          "<br>logsize_t1: ", round(as.numeric(logsize_t1), 3),
          ifelse(is_misfit, paste0("<br><b>[", status, "]</b>"), "")
        )
      )
    
    g <- ggplot(dat, aes(x = logsize_t0, y = logsize_t1, text = tooltip)) +
      geom_point(aes(color = color), alpha = 0.7, pch = 16, size = 0.7) +
      scale_color_identity() +
      theme_bw(base_size = 11) +
      labs(x = "Log(Size t0)", y = "Log(Size t1)") +
      theme(panel.grid.major = element_line(color = "grey85"),
            panel.grid.minor = element_line(color = "grey92"),
            legend.position = "none")
    
    ggplotly(g, tooltip = "text") |> event_register("plotly_click")
  })
  
  # Output row info
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
  
  # Download CSV
  output$download_misfits <- downloadHandler(
    filename = function() {
      paste0("misfit_id_quad_years_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(misfit_tbl(), file)
    }
  )
  
  # Upload CSV
  observeEvent(input$upload_misfits, {
    req(input$upload_misfits)
    tryCatch({
      new_df <- readr::read_csv(input$upload_misfits$datapath, show_col_types = FALSE)
      stopifnot(all(c("id_quad_year", "status") %in% colnames(new_df)))
      old_df <- misfit_tbl()
      merged <- bind_rows(new_df, old_df %>% filter(!id_quad_year %in% new_df$id_quad_year))
      misfit_tbl(merged)
    }, error = function(e) {
      showNotification("Error: CSV must have columns 'id_quad_year' and 'status'.", type = "error")
    })
  })
}

shinyApp(ui, server)
