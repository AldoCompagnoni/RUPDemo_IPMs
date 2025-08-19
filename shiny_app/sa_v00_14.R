# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.19

# Version 00_14: Misfits list with two-row layout

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
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry) %>%
  mutate(id_quad_year = paste0(track_id, "_", quad, "_", year))

# ---- SHINY APP ----
ui <- fluidPage(
  fluidRow(
    column(
      width = 5,
      selectizeInput(
        "exclude_input",
        "Add to misfits list (id_quad_year):",
        choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
        selected = "",
        multiple = FALSE,
        options = list(placeholder = 'Type to add…', create = FALSE)
      ),
      actionButton("add_exclude", "Add", class = "btn btn-sm btn-danger"),
      br(), br(),
      
      tags$div(
        class = "panel panel-default",
        tags$div(
          class = "panel-heading",
          tags$h5(class = "panel-title",
                  tags$a(href = "#collapseExcluded", "data-toggle" = "collapse",
                         tags$small("Misfits list")))
        ),
        tags$div(
          id = "collapseExcluded", class = "panel-collapse collapse in",
          tags$div(class = "panel-body", uiOutput("excluded_list"))
        )
      ),
      br(),
      downloadButton("download_excluded", "Download misfits CSV", class = "btn btn-sm btn-primary"),
      br(), br(),
      fileInput(
        "upload_excluded",
        "Upload misfits CSV:",
        accept = c(".csv"),
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      )
    ),
    
    column(
      width = 7,
      plotlyOutput("plt", height = "600px"),
      br(),
      verbatimTextOutput("clicked_id")
    )
  )
)

server <- function(input, output, session) {
  
  # --- Misfits data frame ---
  misfits_df <- reactiveVal(
    tibble(id_quad_year = character(),
           comment = character(),
           status  = character(),
           group   = character())
  )
  
  # Add ID to misfits list
  observeEvent(input$add_exclude, {
    req(input$exclude_input)
    cur <- misfits_df()
    if (!(input$exclude_input %in% cur$id_quad_year)) {
      misfits_df(bind_rows(cur, tibble(
        id_quad_year = input$exclude_input,
        comment = "",
        status = "include",
        group = "A"
      )))
    }
  })
  
  # ---- Misfits list UI ----
  output$excluded_list <- renderUI({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return("None")
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    tagList(
      lapply(seq_len(nrow(df_ex)), function(i) {
        this_id  <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        this_stat <- df_ex$status[i]
        this_group <- df_ex$group[i]
        sid <- safe_id(this_id)
        
        tags$div(
          style = "padding:8px 0; border-bottom:1px solid #ddd;",
          # Row 1
          fluidRow(
            column(5,
                   tags$b(this_id),
                   tags$div(style = "font-size: 90%; color:#555; margin-top:2px;",
                            ifelse(nchar(this_com) > 0, this_com, tags$em("(no comment)")))
            ),
            column(7, style = "text-align:right;",
                   actionButton(paste0("edit_", sid), "Edit comment",
                                class = "btn btn-xs btn-primary",
                                style = "margin-right:6px;"),
                   actionButton(paste0("remove_", sid), "Remove",
                                class = "btn btn-xs btn-default")
            )
          ),
          # Row 2
          fluidRow(
            column(6,
                   selectInput(paste0("status_", sid), "Status:",
                               choices = c("include", "exclude"),
                               selected = this_stat,
                               width = "100%")
            ),
            column(6,
                   selectInput(paste0("group_", sid), "Group:",
                               choices = c("A", "B", "C"),
                               selected = this_group,
                               width = "100%")
            )
          )
        )
      })
    )
  })
  
  # ---- Edit comment (modal) ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    lapply(seq_len(nrow(df_ex)), function(i) {
      local({
        this_id  <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        
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
          
          observeEvent(input$modal_save_comment, {
            df2 <- misfits_df()
            df2$comment[df2$id_quad_year == this_id] <- input$modal_comment_input %||% ""
            misfits_df(df2)
            removeModal()
          }, once = TRUE, ignoreInit = TRUE)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # ---- Remove from misfits ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    lapply(seq_len(nrow(df_ex)), function(i) {
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        
        observeEvent(input[[paste0("remove_", sid)]], {
          df2 <- misfits_df()
          df2 <- df2[df2$id_quad_year != this_id, ]
          misfits_df(df2)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # ---- Sync status + group inputs ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    
    lapply(seq_len(nrow(df_ex)), function(i) {
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        
        observeEvent(input[[paste0("status_", sid)]], {
          df2 <- misfits_df()
          df2$status[df2$id_quad_year == this_id] <- input[[paste0("status_", sid)]]
          misfits_df(df2)
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("group_", sid)]], {
          df2 <- misfits_df()
          df2$group[df2$id_quad_year == this_id] <- input[[paste0("group_", sid)]]
          misfits_df(df2)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  grow_df_filtered <- reactive({
    ex_ids <- misfits_df() %>% filter(status == "exclude") %>% pull(id_quad_year)
    grow_df %>% filter(!id_quad_year %in% ex_ids)
  })
  
  # ---- Capture clicked point ----
  selected_point <- reactiveVal(NULL)
  
  observeEvent(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    sel <- grow_df_filtered()[click$pointNumber + 1, ]
    selected_point(sel)
    clipr::write_clip(sel[, c("track_id","quad","year")])  # copy to clipboard
  })
  
  # ---- Scatter plot ----
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat) > 0, "No rows match the current filter."))
    
    # Flag misfits
    misfits_ids <- misfits_df()$id_quad_year
    dat <- dat %>%
      mutate(is_misfit = id_quad_year %in% misfits_ids)
    
    # Tooltip
    dat <- dat %>%
      mutate(tooltip = paste0(
        "id_quad_year: ", id_quad_year,
        "<br>track_id: ", as.character(track_id),
        "<br>logsize_t0: ", round(as.numeric(logsize_t0),3),
        "<br>logsize_t1: ", round(as.numeric(logsize_t1),3)
      ))
    
    # Base plot
    g <- ggplot(dat, aes(x = logsize_t0, y = logsize_t1, group = id_quad_year)) +
      geom_point(aes(color = is_misfit, alpha = is_misfit), size = 0.7, pch = 16) +
      scale_color_manual(values = c("FALSE"="#D55E00","TRUE"="grey50")) +
      scale_alpha_manual(values = c("FALSE"=0.7,"TRUE"=0.3)) +
      theme_bw(base_size = 11) +
      labs(x="Log(Size t0)", y="Log(Size t1)") +
      theme(panel.grid.major = element_line(color = "grey85"),
            panel.grid.minor = element_line(color = "grey92"),
            legend.position = "none")
    
    # Highlight selected point
    sel <- selected_point()
    if(!is.null(sel)){
      g <- g + geom_point(data = sel, aes(x = logsize_t0, y = logsize_t1),
                          color = "black", size = 2, shape = 17)
    }
    
    # Convert to plotly with tooltip mapping
    ggplotly(g, tooltip = c("x","y")) %>%
      event_register("plotly_click")
  })
  
  # ---- Show clicked row ----
  output$clicked_id <- renderPrint({
    selected_point()
  })
  
  # ---- Download misfits ----
  output$download_excluded <- downloadHandler(
    filename = function() {
      paste0("misfits_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df_out <- misfits_df() %>% mutate(across(everything(), ~replace_na(., "")))
      readr::write_csv(df_out, file, na = "")
    }
  )
  
  # ---- Upload misfits ----
  observeEvent(input$upload_excluded, {
    req(input$upload_excluded)
    tryCatch({
      new_df <- readr::read_csv(input$upload_excluded$datapath, show_col_types = FALSE)
      if (!all(c("id_quad_year", "comment", "status", "group") %in% names(new_df))) {
        stop("CSV must have columns 'id_quad_year', 'comment', 'status', 'group'.")
      }
      new_df <- new_df %>%
        mutate(id_quad_year = as.character(id_quad_year),
               comment = replace_na(as.character(comment), ""),
               status = ifelse(status %in% c("include","exclude"), status, "include"),
               group = ifelse(group %in% c("A","B","C"), group, "A")) %>%
        distinct(id_quad_year, .keep_all = TRUE)
      
      cur <- misfits_df()
      merged <- cur %>%
        anti_join(new_df, by = "id_quad_year") %>%
        bind_rows(new_df) %>%
        arrange(id_quad_year)
      misfits_df(merged)
      showNotification("Misfits list updated from CSV.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
    })
  })
}

shinyApp(ui, server)
