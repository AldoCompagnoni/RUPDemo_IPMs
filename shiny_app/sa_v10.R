# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.15

# Version 10: Extends csv to dataframe with a comment section

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
  
  if (nrow(t1_poly) > 0) {
    year_t1 <- plot_polygons(t1_poly) + labs(title = paste("Year", yr + 1))
  } else {
    year_t1 <- ggplot() +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"])
      ) +
      theme_bw() + theme_void() +
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
        "exclude_input",
        "Exclude id_quad_year:",
        choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
        selected = "",
        multiple = FALSE,
        options = list(placeholder = 'Type to exclude…', create = FALSE)
      ),
      actionButton("add_exclude", "Add to exclude list", class = "btn btn-sm btn-danger"),
      br(), br(),
      
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          tags$h5(class = "panel-title",
                  tags$a(href = "#collapseExcluded", "data-toggle" = "collapse",
                         tags$small("Excluded id_quad_years & comments")))
        ),
        tags$div(
          id = "collapseExcluded", class = "panel-collapse collapse",
          tags$div(class = "panel-body", uiOutput("excluded_list"))
        )
      ),
      br(),
      downloadButton("download_excluded", "Download excluded CSV", class = "btn btn-sm btn-primary"),
      br(), br(),
      fileInput(
        "upload_excluded",
        "Upload excluded CSV:",
        accept = c(".csv"),
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      )
    ),
    
    column(
      width = 8,
      plotlyOutput("plt", height = "400px"),
      plotOutput("before_after_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  # --- Excluded data: pure data frame (no HTML), ensures clean CSV ---
  excluded_df <- reactiveVal(tibble(id_quad_year = character(), comment = character()))
  
  # Add selected ID to exclude list (no comment yet)
  observeEvent(input$add_exclude, {
    req(input$exclude_input)
    cur <- excluded_df()
    if (!(input$exclude_input %in% cur$id_quad_year)) {
      excluded_df(bind_rows(cur, tibble(id_quad_year = input$exclude_input, comment = "")))
    }
  })
  
  # ---- Excluded list UI (ID + comment text + Edit + Remove) ----
  output$excluded_list <- renderUI({
    df_ex <- excluded_df()
    if (nrow(df_ex) == 0) return("None")
    
    # helper to keep IDs safe for input names
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    tagList(
      lapply(seq_len(nrow(df_ex)), function(i) {
        this_id  <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        
        tags$div(
          style = "padding:6px 0; border-bottom:1px solid #eee;",
          fluidRow(
            column(7,
                   tags$b(this_id),
                   tags$div(style = "font-size: 90%; color:#555; margin-top:2px;",
                            ifelse(nchar(this_com) > 0, this_com, tags$em("(no comment)")))
            ),
            column(5, style = "text-align:right;",
                   actionButton(paste0("edit_", sid), "Edit comment",
                                class = "btn btn-xs btn-primary",
                                style = "margin-right:6px;"),
                   actionButton(paste0("remove_", sid), "Remove",
                                class = "btn btn-xs btn-default")
            )
          )
        )
      })
    )
  })
  
  # ---- Edit comment (modal) ----
  observe({
    df_ex <- excluded_df()
    if (nrow(df_ex) == 0) return(NULL)
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    lapply(seq_len(nrow(df_ex)), function(i) {
      local({  # important: create a new environment for each row
        this_id  <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        
        # Edit button
        observeEvent(input[[paste0("edit_", sid)]], {
          showModal(modalDialog(
            title = paste("Edit comment for", this_id),
            textAreaInput("modal_comment_input", "Comment", value = this_com, width = "100%", rows = 6),
            footer = tagList(
              modalButton("Cancel"),
              actionButton("modal_save_comment", "Save", class = "btn btn-primary")
            ),
            easyClose = FALSE
          ))
          
          # Scoped Save handler
          observeEvent(input$modal_save_comment, {
            df2 <- excluded_df()
            df2$comment[df2$id_quad_year == this_id] <- input$modal_comment_input %||% ""
            excluded_df(df2)
            removeModal()
          }, once = TRUE, ignoreInit = TRUE)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # ---- Remove excluded ID ----
  observe({
    df_ex <- excluded_df()
    if (nrow(df_ex) == 0) return(NULL)
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    lapply(seq_len(nrow(df_ex)), function(i) {
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        
        observeEvent(input[[paste0("remove_", sid)]], {
          df2 <- excluded_df()
          df2 <- df2[df2$id_quad_year != this_id, ]
          excluded_df(df2)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # ---- Filtered data for plot (exclude rows present in excluded_df) ----
  grow_df_filtered <- reactive({
    q <- input$track_query
    dat <- grow_df
    if (!is.null(q) && trimws(q) != "") {
      dat <- dat %>%
        mutate(track_id_chr = as.character(track_id)) %>%
        filter(str_detect(track_id_chr, fixed(q, ignore_case = TRUE)))
    }
    # apply exclusions
    ex_ids <- excluded_df()$id_quad_year
    dat %>% filter(!id_quad_year %in% ex_ids)
  })
  
  # ---- Scatter plot ----
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat) > 0, "No rows match the current filter."))
    
    g <- dat %>%
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
    
    ggplotly(g, tooltip = "text") |> event_register("plotly_click")
  })
  
  # ---- Row detail from click + clipboard copy ----
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
  
  # ---- Before/After polygons ----
  output$before_after_plot <- renderPlot({
    r <- row_r()
    req(r)
    plot_before_after(as.character(r$track_id), as.character(r$quad), as.numeric(r$year))
  })
  
  # ---- Download excluded (clean CSV) ----
  output$download_excluded <- downloadHandler(
    filename = function() {
      paste0("excluded_id_quad_years_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # ensure empty strings instead of NA in CSV
      df_out <- excluded_df() %>%
        mutate(across(everything(), ~replace_na(., "")))
      readr::write_csv(df_out, file, na = "")
    }
  )
  
  # ---- Upload excluded (merge by id, keep latest comment; ensure columns exist) ----
  observeEvent(input$upload_excluded, {
    req(input$upload_excluded)
    tryCatch({
      new_df <- readr::read_csv(input$upload_excluded$datapath, show_col_types = FALSE)
      if (!all(c("id_quad_year", "comment") %in% names(new_df))) {
        stop("CSV must have columns 'id_quad_year' and 'comment'.")
      }
      new_df <- new_df %>%
        mutate(
          id_quad_year = as.character(id_quad_year),
          comment = replace_na(as.character(comment), "")
        ) %>%
        distinct(id_quad_year, .keep_all = TRUE)
      
      # merge with existing, prefer uploaded comment when there is a conflict
      cur <- excluded_df()
      merged <- cur %>%
        anti_join(new_df, by = "id_quad_year") %>%
        bind_rows(new_df) %>%
        arrange(id_quad_year)
      
      excluded_df(merged)
      showNotification("Excluded list updated from CSV.", type = "message")
    }, error = function(e) {
      showNotification(
        paste("Error reading CSV:", e$message),
        type = "error"
      )
    })
  })
}

shinyApp(ui, server)
