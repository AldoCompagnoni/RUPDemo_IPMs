# Shiny app data

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.09.05

# Version 21: Colored "other individual" polygons by track_id, relinked zoom on
  # maps, defined "groups" for misfits list

# ----------------------------------------------------------------------
# PACKAGES
# ----------------------------------------------------------------------
source('helper_functions/load_packages.R')
load_packages(tidyverse, sf, shiny, plotly, gridExtra, clipr)

# ----------------------------------------------------------------------
# SPECIFICATION
# ----------------------------------------------------------------------
v_author_year <- 'anderson_2016'
v_region_abb  <- 'mt'
v_species     <- 'Bouteloua gracilis'

v_size_threshold <- NULL

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

# ----------------------------------------------------------------------
# DIRECTORIES
# ----------------------------------------------------------------------
dir_pub    <- file.path(paste0(v_author_year, '_', v_region_abb))
dir_R      <- file.path(dir_pub, 'R', v_sp_abb)
dir_data   <- file.path(dir_pub, 'data', v_sp_abb)
dir_result <- file.path(dir_pub, 'results', v_sp_abb)

dir.create(dir_R, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_result, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# DATA
# ----------------------------------------------------------------------
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

# Apply size threshold (if provided)
if (!is.null(v_size_threshold)) {
  df <- df %>%
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0),
           logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

# Growth data (non-zero)
grow_df <- df %>%
  subset(size_t0 != 0 & size_t1 != 0) %>%
  select(site, quad, track_id, track_id_quad, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry, id_quad_year)

bbox <- st_bbox(df)

# ----------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------
# Generate a placeholder plot when no spatial data is available.
empty_plotly <- function(src) {
  p <- ggplot() + 
    # no geometries (geom_blank) and no axes (theme_void)
    geom_blank() + theme_void()
  # Convert to an interactive plotly object
  # Posts an event hook for plotly_relayout so interactions are still handled
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

# Draw polygons (plant outlines) spatially for a given quadrat and year
plot_polygons_with_others <- function(df_selected, quad_val, year_val, show_buffer = TRUE, show_others = TRUE, src) {
  if (nrow(df_selected) == 0 || !"geometry" %in% names(df_selected)) {
    # If no data or missing geometry, show message “No data available” with bounding box extents
    p <- ggplot() +
      annotate("text", x = mean(bbox[c("xmin","xmax")]), y = mean(bbox[c("ymin","ymax")]),
               label = "No data available") +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
      theme_bw()
    return(ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout"))
  }
  
  # Cast and buffer (use a small, fixed buffer to avoid bbox-dependent errors)
  suppressWarnings({ df_selected <- 
    # Geometry casting
    st_cast(df_selected, "MULTIPOLYGON") })
  if (show_buffer) {
    # Buffer of 0.05 since the plot is at x and ymax == 1
    df_buffer <- tryCatch(st_buffer(df_selected, dist = 0.05), error = function(e) NULL)
  } else df_buffer <- NULL
  
  # Extract other plants in the quadrat for the same year
  df_others <- df %>% filter(quad == quad_val, year == year_val, !(track_id %in% df_selected$track_id))
  # Cast them to polygons
  if (nrow(df_others) > 0) suppressWarnings({ df_others <- st_cast(df_others, "MULTIPOLYGON") })
  
  # Final layer stack
  p <- ggplot()
  if (show_buffer && !is.null(df_buffer)) p <- p + geom_sf(data = df_buffer, fill = "lightblue", alpha = 0.3)
  if (show_others && nrow(df_others) > 0) p <- p + geom_sf(data = df_others, aes( fill = track_id ), alpha = 0.5 ) + guides(fill = "none" )
  p <- p + geom_sf(data = df_selected, fill = "dodgerblue") +
    # Ensure plots always align with bounding box of original dataset
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
    theme_bw()
  
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}


# Define how individuals are shown in the scatter
.state_labels <- c(
  `0` = "Filtered out",
  `1` = "Default",
  `2` = "Selected",
  `3` = "Misfit include",
  `4` = "Misfit exclude"
)

.state_colors <- c(
  `0` = "grey70",      # grayed out
  `1` = "#D55E00",     # default orange/red
  `2` = "#8B5CF6",     # purple (violet-500)
  `3` = "#FDAE6B",     # lighter orange
  `4` = "#FFFFFF"       # white
)


# ----------------------------------------------------------------------
# SHINY UI
# ----------------------------------------------------------------------
# The UI is organized in fluidPage → fluidRow → 2 columns
# We define the UI at the colum level
ui <- fluidPage(
  fluidRow(
    # Column 1 (Filters + Misfits + Controls)
    column(width = 5,
           # Filter panel:
           tags$div(style = "font-size:90%;", class = "panel panel-default",
                    tags$div(class = "panel-heading",
                             tags$h5(class = "panel-title",
                                     tags$a(href = "#collapseFilters", "data-toggle" = "collapse",
                                            tags$strong("Filters (set state column)")
                                     )
                             )
                    ),
                    tags$div(id = "collapseFilters", class = "panel-collapse collapse",
                             tags$div(class = "panel-body",
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
                                                  min = min(grow_df$year), max = max(grow_df$year),
                                                  value = c(min(grow_df$year), max(grow_df$year)), step = 1, sep = "")
                             )
                    )
           ),
           hr(),
           
           # Selected Point Info box:
           tags$div(class = "panel panel-default",
                    tags$div(class = "panel-heading",
                             tags$h5(class = "panel-title", tags$strong("Selected point (purple = state 2)"))
                    ),
                    tags$div(class = "panel-body", verbatimTextOutput("row"))
           ),
           
           # Misfit inputs:
           selectizeInput("misfit_input", "Add to misfits list (id_quad_year):",
                          choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
                          selected = "", multiple = FALSE),
           actionButton("add_misfit", "Add", class = "btn btn-sm btn-danger"),
           br(), br(),
           
           # Misfits panel:
           tags$div(class = "panel panel-default",
                    tags$div(class = "panel-heading",
                             tags$h5(class = "panel-title",
                                     tags$a(href = "#collapseMisfits", "data-toggle" = "collapse",
                                            tags$small("Misfits list (state 3 = include, 4 = exclude)")
                                     )
                             )
                    ),
                    tags$div(id = "collapseMisfits", class = "panel-collapse collapse",
                             tags$div(class = "panel-body", uiOutput("misfits_list"))
                    )
           ),
           br(),
           
           # Download misfits CSV:
           downloadButton("download_misfits", "Download misfits CSV", class = "btn btn-sm btn-primary"),
           br(), br(),
           
           # Upload misfits CSV:
           fileInput("upload_misfits", "Upload misfits CSV:", accept = c(".csv"),
                     buttonLabel = "Browse...", placeholder = "No file selected")
    ),
    
    # Column 2 (Plots + Maps)
    column(width = 7,
           
           # Jitter control panel:
           tags$div(class = "panel panel-default", style = "font-size:85%; margin-bottom:10px;",
                    tags$div(class = "panel-heading",
                             tags$h5(class = "panel-title",
                                     tags$a(href = "#collapseJitter", "data-toggle" = "collapse",
                                            tags$strong("Jitter controls")))
                    ),
                    tags$div(id = "collapseJitter", class = "panel-collapse collapse",
                             tags$div(class = "panel-body",
                                      fluidRow(
                                        column(6, sliderInput("jit_w", "Width (x-axis):", min = 0, max = 0.5, value = 0, step = 0.01)),
                                        column(6, sliderInput("jit_h", "Height (y-axis):", min = 0, max = 0.5, value = 0, step = 0.01))
                                      )
                             )
                    )
           ),
           
           # Scatter plot:
           plotlyOutput("plt", height = "400px"),
           br(),
           
           # Geometry plots
           tags$div(class = "panel panel-default", style = "padding:10px; margin-top:12px;",
                    
                    # Another level of fluidRow → 2 columns
                    fluidRow(
                      column(6,
                             # Another level of fluidRow → 2 columns
                             fluidRow(column(6, uiOutput("map_t0_title")),
                                      # Buffer checkboxes:
                                      column(6, checkboxInput("show_buffer", "Show buffer", value = TRUE))),
                             plotlyOutput("map_t0", height = "400px")
                      ),
                      column(6,
                             fluidRow(column(6, uiOutput("map_t1_title")),
                                      column(6, checkboxInput("show_others", "Show other individuals", value = TRUE))),
                             plotlyOutput("map_t1", height = "400px")
                      )
                    )
           )
    )
  )
)


# ----------------------------------------------------------------------
# SHINY SERVER
# ----------------------------------------------------------------------
server <- function(input, output, session) {
  
  # -------------------------
  # Misfits reactive storage
  # -------------------------

  # A reactive state holder for the misfit list
  misfits_df <- reactiveVal(tibble(id_quad_year = character(), comment = character(), status = character(), group = character()))
  
  # When user clicks actionButton("add_misfit" [...])
  observeEvent(input$add_misfit, {
    # Check for non-null
    req(input$misfit_input)
    cur <- misfits_df()
    # If not already in misfits_df, a new row is added
    if (!(input$misfit_input %in% cur$id_quad_year)) {
      misfits_df(bind_rows(tibble(id_quad_year = input$misfit_input, comment = "", status = "include", group = "A"), cur))
    }
  })
  
  # -------------------------
  # Misfits list UI
  # -------------------------
  output$misfits_list <- renderUI({
    df_ex <- misfits_df()
    # If the list is empty: Returns “None”
    if (nrow(df_ex) == 0) return("None")
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    tagList(
      lapply(seq_len(nrow(df_ex)), function(i){
        this_id <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        this_stat <- df_ex$status[i]
        this_group <- df_ex$group[i]
        sid <- safe_id(this_id)
        
        # A visually distinct div for each item (with a border and padding)
        tags$div(
          style = "padding:8px 0; border-bottom:1px solid #ddd;",
          # First row
          fluidRow(
            column(5, tags$b(this_id),
                   tags$div(style = "font-size:90%; color:#555; margin-top:2px;",
                            ifelse(nchar(this_com) > 0, this_com, tags$em("(no comment)")))) ,
            column(7, style = "text-align:right;",
                   actionButton(paste0("edit_", sid), "Edit comment", class = "btn btn-xs btn-primary", style = "margin-right:6px;"),
                   actionButton(paste0("remove_", sid), "Remove", class = "btn btn-xs btn-default"))
          ),
          fluidRow(
            column(6, selectInput(paste0("status_", sid), "Status:", choices = c("include", "exclude"), selected = this_stat, width = "100%")),
            column(6, selectInput(paste0("group_", sid), "Group:", choices = c("A", "AM", "AP", "B", "BM", "C", "CM", "S", "M", "Other..."), selected = this_group, width = "100%"))
          )
        )
      })
    )
  })
  
  # Edit comment modal
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("edit_", sid)]], {
          showModal(modalDialog(title = paste("Edit comment for", this_id),
                                textAreaInput("modal_comment_input", "Comment", value = this_com, width = "100%", rows = 6),
                                footer = tagList(modalButton("Cancel"), actionButton("modal_save_comment", "Save", class = "btn btn-primary")),
                                easyClose = FALSE))
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
  
  # Remove misfit
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("remove_", sid)]], {
          df2 <- misfits_df()
          misfits_df(df2[df2$id_quad_year != this_id, ])
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # Sync status & group
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    lapply(seq_len(nrow(df_ex)), function(i){
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
  
  # -------------------------
  # Selected point (id only)
  # -------------------------
  selected_id <- reactiveVal(NULL)
  
  observeEvent(event_data("plotly_click", source = "A"), {
    click <- event_data("plotly_click", source = "A")
    req(click$key)
    selected_id(as.character(click$key))
  })
  
  # -------------------------
  # Build state-coded dataset
  # -------------------------
  # State rules:
  # 1. Everyone starts as 1 (Default, orange)
  # 2. Filters set non-matching rows to 0 (grayed out)
  # 3. Misfits override filters: include -> 3 (lighter orange), exclude -> 4 (white)
  # 4. Clicked point overrides all -> 2 (purple)
  state_data <- reactive({
    dat <- grow_df %>% mutate(state_code = 1L) # start all at 1
    
    # Active filters
    matches <- rep(TRUE, nrow(dat))
    if (!is.null(input$track_query) && input$track_query != "") {
      matches <- matches & (dat$track_id_quad == input$track_query)
    }
    if (!is.null(input$site_query) && input$site_query != "") {
      matches <- matches & (dat$site == input$site_query)
    }
    if (!is.null(input$quad_query) && input$quad_query != "") {
      matches <- matches & (as.character(dat$quad) == input$quad_query)
    }
    if (!is.null(input$year_query) && length(input$year_query) == 2) {
      matches <- matches & (dat$year >= input$year_query[1] & dat$year <= input$year_query[2])
    }
    
    # Apply filter effect: set to 0 when not matching
    dat$state_code[!matches] <- 0L
    
    # Apply misfits override
    mf <- misfits_df()
    if (nrow(mf) > 0) {
      idx <- match(dat$id_quad_year, mf$id_quad_year)
      has_mf <- !is.na(idx)
      if (any(has_mf)) {
        dat$state_code[has_mf & mf$status[idx] == "include"] <- 3L
        dat$state_code[has_mf & mf$status[idx] == "exclude"] <- 4L
      }
    }
    
    # Apply selection override (purple)
    sid <- selected_id()
    if (!is.null(sid)) {
      dat$state_code[dat$id_quad_year == sid] <- 2L
    }
    
    # Build aesthetics
    dat <- dat %>%
      mutate(
        state_code_f = factor(as.character(state_code), levels = c("0","1","2","3","4")),
        alpha_val = ifelse(state_code == 0L, 0.2, 0.9),
        state_label = .state_labels[as.character(state_code)],
        tooltip = paste0(
          "id_quad_year: ", id_quad_year,
          "<br>logsize_t0: ", round(logsize_t0, 3),
          "<br>logsize_t1: ", round(logsize_t1, 3),
          "<br>size_t0: ", round(size_t0, 3),
          "<br>size_t1: ", round(size_t1, 3),
          "<br>State: ", state_label
        )
      )
    
    dat
  })
  
  # -------------------------
  # Scatter plot
  # -------------------------
  output$plt <- renderPlotly({
    dat <- state_data()
    
    g <- ggplot(dat, aes(x = logsize_t0, y = logsize_t1)) +
      geom_jitter(aes(color = state_code_f, alpha = alpha_val, text = tooltip, key = id_quad_year),
                  width = input$jit_w, height = input$jit_h, size = 0.7) +
      scale_color_manual(values = .state_colors, breaks = names(.state_labels), labels = .state_labels) +
      scale_alpha_identity() +
      theme_bw(base_size = 11) +
      theme(panel.grid.major = element_line(color = "grey85"),
            panel.grid.minor = element_line(color = "grey92"),
            legend.position = "none") +
      labs(x = "Log(Size t0)", y = "Log(Size t1)")
    
    ggplotly(g, tooltip = "text", source = "A") |> event_register("plotly_click")
  })
  
  # -------------------------
  # Selected row printout
  # -------------------------
  output$row <- renderPrint({
    sid <- selected_id()
    if (is.null(sid)) {
      cat("No point selected")
    } else {
      sel <- grow_df %>% filter(id_quad_year == sid) %>% st_drop_geometry() %>%
        select(id_quad_year, track_id, quad, year, logsize_t0, logsize_t1, size_t0, size_t1, survives)
      print(sel)
    }
  })
  
  # -------------------------
  # Maps (t0 and t1)
  # -------------------------
  output$map_t0_title <- renderUI({ tags$strong("Year t0") })
  output$map_t1_title <- renderUI({ tags$strong("Year t1") })
  
  output$map_t0 <- renderPlotly({
    sid <- selected_id()
    if (is.null(sid)) return(empty_plotly("t0"))
    sel <- grow_df %>% filter(id_quad_year == sid) %>% st_drop_geometry() %>% slice(1)
    plot_polygons_with_others(df %>% filter(track_id == sel$track_id, quad == sel$quad, year == sel$year),
                              sel$quad, sel$year, show_buffer = input$show_buffer, show_others = input$show_others, src = "t0")
  })
  
  output$map_t1 <- renderPlotly({
    sid <- selected_id()
    if (is.null(sid)) return(empty_plotly("t1"))
    sel <- grow_df %>% filter(id_quad_year == sid) %>% st_drop_geometry() %>% slice(1)
    plot_polygons_with_others(df %>% filter(track_id == sel$track_id, quad == sel$quad, year == sel$year + 1),
                              sel$quad, sel$year + 1, show_buffer = input$show_buffer, show_others = input$show_others, src = "t1")
  })
  
  # Link zoom between maps
  observeEvent(event_data("plotly_relayout", source="t0"),{
    rel <- event_data("plotly_relayout", source="t0")
    if(!is.null(rel)) plotlyProxy("map_t1", session) %>% plotlyProxyInvoke("relayout", rel)
    }, ignoreNULL=TRUE)
  
  observeEvent(event_data("plotly_relayout", source="t1"),{
    rel <- event_data("plotly_relayout", source="t1")
    if(!is.null(rel)) plotlyProxy("map_t0", session) %>% plotlyProxyInvoke("relayout", rel)
    }, ignoreNULL=TRUE)
  
  # -------------------------
  # Download / Upload misfits
  # -------------------------
  output$download_misfits <- downloadHandler(
    filename = function() { paste0("misfits_", Sys.Date(), ".csv") },
    content = function(file) { write_csv(misfits_df(), file) }
  )
  
  observeEvent(input$upload_misfits, {
    req(input$upload_misfits)
    tryCatch({
      new_df <- readr::read_csv(input$upload_misfits$datapath, show_col_types = FALSE)
      stopifnot(all(c("id_quad_year", "status") %in% names(new_df)))
      old_df <- misfits_df()
      merged <- bind_rows(new_df, old_df %>% filter(!id_quad_year %in% new_df$id_quad_year))
      misfits_df(merged)
    }, error = function(e) {
      showNotification("Error: CSV must have columns 'id_quad_year' and 'status'.", type = "error")
    })
  })
}

# ----------------------------------------------------------------------
# RUN APP
# ----------------------------------------------------------------------
shinyApp(ui, server)
