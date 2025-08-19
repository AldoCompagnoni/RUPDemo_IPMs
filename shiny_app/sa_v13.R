# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri
# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.18

# Version 13: Added site, track_id_quad, quad, and year filters/slider

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
  mutate(across(c(quad, track_id), as.factor)) %>%
  rename(size_t0  = basal_area_genet,
         size_t1  = size_tplus1,
         survives = survives_tplus1) %>%
  group_by(track_id, quad) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    next_year = lead(year),
    size_t1   = if_else(!is.na(next_year) & next_year == year + 1, size_t1, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    logsize_t0   = log(size_t0),
    logsize_t1   = log(size_t1),
    logsize_t0_2 = logsize_t0^2,
    logsize_t0_3 = logsize_t0^3,
    site = sub("^.{3}_([A-Za-z]+)_.*$", "\\1", track_id),        # extract site
    track_id_quad = paste0(track_id, "_", quad)                  # combined track+quad
  )

# Size threshold
if (length(v_size_threshold) > 0) {
  df <- df %>%
    filter(logsize_t0 > v_size_threshold | is.na(logsize_t0)) %>%
    filter(logsize_t1 > v_size_threshold | is.na(logsize_t1))
}

# Growth dataframe
grow_df <- df %>%
  subset(size_t0 != 0 & size_t1 != 0) %>%
  select(site, quad, track_id, track_id_quad, year, size_t0, survives, size_t1,
         logsize_t0, logsize_t1, logsize_t0_2, logsize_t0_3, geometry) %>%
  mutate(id_quad_year = paste0(track_id, "_", quad, "_", year))

# Bounding box
bbox <- st_bbox(df)

# ---- HELPERS ----
empty_plotly <- function(src) {
  p <- ggplot() + geom_blank() + theme_void()
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

plot_polygons_with_others <- function(df_selected, quad_val, year_val,
                                      show_buffer = TRUE, show_others = TRUE, src) {
  if (nrow(df_selected) == 0 || !"geometry" %in% names(df_selected)) {
    p <- ggplot() +
      annotate("text",
               x = mean(bbox[c("xmin","xmax")]),
               y = mean(bbox[c("ymin","ymax")]),
               label = "No data available") +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
               ylim = c(bbox["ymin"], bbox["ymax"])) +
      theme_bw()
    return(ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout"))
  }
  
  df_selected <- st_cast(df_selected, "MULTIPOLYGON")
  
  if (show_buffer) {
    buff <- if_else(st_bbox(df_selected)[3] < 1.1 | st_bbox(df_selected)[3] > 100, 0.05, 5)
    df_buffer <- st_buffer(df_selected, dist = buff)
  }
  
  df_others <- df %>%
    filter(quad == quad_val, year == year_val, !(track_id %in% df_selected$track_id))
  if (nrow(df_others) > 0) df_others <- st_cast(df_others, "MULTIPOLYGON")
  
  p <- ggplot()
  if (show_buffer) p <- p + geom_sf(data = df_buffer, fill = "lightblue", alpha = 0.3)
  if (show_others && nrow(df_others) > 0) p <- p + geom_sf(data = df_others, fill = "lightgrey", alpha = 0.5)
  p <- p + geom_sf(data = df_selected, fill = "dodgerblue") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"])) +
    theme_bw()
  
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

plot_before_after_plotly <- function(track_id_val, quad_val, yr_val,
                                     show_buffer, show_others) {
  
  year_t0 <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr_val) %>%
    select(track_id, geometry)
  year_t1 <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr_val+1) %>%
    select(track_id, geometry)
  
  list(
    t0 = plot_polygons_with_others(year_t0, quad_val, yr_val, show_buffer, show_others, src = "t0"),
    t1 = plot_polygons_with_others(year_t1, quad_val, yr_val+1, show_buffer, show_others, src = "t1")
  )
}

# ---- SHINY UI ----
ui <- fluidPage(
  fluidRow(
    column(
      width = 4,
      
      # ---- FILTERS PANEL ----
      tags$div(class="panel panel-default",
               tags$div(class="panel-heading",
                        tags$h5(class="panel-title",
                                tags$a(href="#collapseFilters","data-toggle"="collapse",
                                       tags$strong("Filters")))),
               tags$div(id="collapseFilters", class="panel-collapse collapse",  # starts collapsed
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
                                             step = 1, sep = ""),
                                 verbatimTextOutput("row")
                        )
               )
      ),
      
      hr(),
      
      # ---- EXCLUDE PANEL (already foldable) ----
      selectizeInput("exclude_input", "Exclude id_quad_year:",
                     choices = c("", sort(unique(as.character(grow_df$id_quad_year)))),
                     selected = "", multiple = FALSE),
      actionButton("add_exclude", "Add to exclude list", class="btn btn-sm btn-danger"),
      br(), br(),
      tags$div(class="panel panel-default",
               tags$div(class="panel-heading",
                        tags$h5(class="panel-title",
                                tags$a(href="#collapseExcluded","data-toggle"="collapse",
                                       tags$small("Excluded id_quad_years & comments")))),
               tags$div(id="collapseExcluded", class="panel-collapse collapse",  # also collapsed by default
                        tags$div(class="panel-body", uiOutput("excluded_list")))),
      br(),
      downloadButton("download_excluded", "Download excluded CSV", class="btn btn-sm btn-primary"),
      br(), br(),
      fileInput("upload_excluded", "Upload excluded CSV:",
                accept=c(".csv"),
                buttonLabel="Browse...", placeholder="No file selected")
    ),
    column(
      width = 8,
      plotlyOutput("plt", height = "400px"),
      fluidRow(
        column(6, plotlyOutput("map_t0", height = "400px")),
        column(6, plotlyOutput("map_t1", height = "400px"))
      ),
      fluidRow(
        column(6, checkboxInput("show_buffer", "Show buffer", value=TRUE)),
        column(6, checkboxInput("show_others", "Show other individuals", value=TRUE))
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  output$map_t0 <- renderPlotly(empty_plotly("t0"))
  output$map_t1 <- renderPlotly(empty_plotly("t1"))
  
  excluded_df <- reactiveVal(tibble(id_quad_year=character(), comment=character()))
  
  observeEvent(input$add_exclude, {
    req(input$exclude_input)
    cur <- excluded_df()
    if (!(input$exclude_input %in% cur$id_quad_year)) {
      excluded_df(bind_rows(cur, tibble(id_quad_year=input$exclude_input, comment="")))
    }
  })
  
  output$excluded_list <- renderUI({
    df_ex <- excluded_df()
    if (nrow(df_ex) == 0) return("None")
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    tagList(
      lapply(seq_len(nrow(df_ex)), function(i){
        this_id <- df_ex$id_quad_year[i]; this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        tags$div(
          style="padding:6px 0; border-bottom:1px solid #eee;",
          fluidRow(
            column(7,
                   tags$b(this_id),
                   tags$div(style="font-size: 90%; color:#555; margin-top:2px;",
                            ifelse(nchar(this_com)>0,this_com,tags$em("(no comment)")))),
            column(5, style="text-align:right;",
                   actionButton(paste0("edit_",sid),"Edit comment", class="btn btn-xs btn-primary", style="margin-right:6px;"),
                   actionButton(paste0("remove_",sid),"Remove", class="btn btn-xs btn-default")
            )
          )
        )
      })
    )
  })
  
  observe({
    df_ex <- excluded_df()
    if (nrow(df_ex)==0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("edit_",sid)]],{
          showModal(modalDialog(
            title=paste("Edit comment for",this_id),
            textAreaInput("modal_comment_input","Comment", value=this_com, width="100%", rows=6),
            footer=tagList(modalButton("Cancel"), actionButton("modal_save_comment","Save",class="btn btn-primary")),
            easyClose=FALSE
          ))
          observeEvent(input$modal_save_comment,{
            df2 <- excluded_df()
            df2$comment[df2$id_quad_year==this_id] <- input$modal_comment_input %||% ""
            excluded_df(df2); removeModal()
          }, once=TRUE, ignoreInit=TRUE)
        }, ignoreInit=TRUE)
        observeEvent(input[[paste0("remove_",sid)]],{
          df2 <- excluded_df()
          df2 <- df2[df2$id_quad_year != this_id, ]
          excluded_df(df2)
        }, ignoreInit=TRUE)
      })
    })
  })
  
  grow_df_filtered <- reactive({
    dat <- grow_df
    if (!is.null(input$track_query) && input$track_query != "") {
      dat <- dat %>% filter(track_id_quad == input$track_query)
    }
    if (!is.null(input$site_query) && input$site_query != "") {
      dat <- dat %>% filter(site == input$site_query)
    }
    if (!is.null(input$quad_query) && input$quad_query != "") {
      dat <- dat %>% filter(as.character(quad)==input$quad_query)
    }
    if (!is.null(input$year_query) && length(input$year_query)==2) {
      dat <- dat %>% filter(year >= input$year_query[1], year <= input$year_query[2])
    }
    ex_ids <- excluded_df()$id_quad_year
    dat %>% filter(!id_quad_year %in% ex_ids)
  })
  
  output$plt <- renderPlotly({
    dat <- grow_df_filtered()
    validate(need(nrow(dat)>0, "No rows match the current filter."))
    g <- dat %>%
      mutate(tooltip=paste0(
        "id_quad_year: ", id_quad_year,
        "<br>track_id: ", as.character(track_id),
        "<br>quad: ", quad,
        "<br>year: ", year,
        "<br>logsize_t0: ", round(as.numeric(logsize_t0),3),
        "<br>logsize_t1: ", round(as.numeric(logsize_t1),3)
      )) %>%
      ggplot(aes(x=logsize_t0, y=logsize_t1, text=tooltip, group=id_quad_year)) +
      geom_point(alpha=0.7, pch=16, size=0.7, color="#D55E00") +
      theme_bw(base_size=11) +
      labs(x="Log(Size t0)", y="Log(Size t1)") +
      theme(panel.grid.major=element_line(color="grey85"),
            panel.grid.minor=element_line(color="grey92"),
            legend.position="none")
    ggplotly(g, tooltip="text") |> event_register("plotly_click")
  })
  
  row_r <- eventReactive(event_data("plotly_click"), {
    click <- event_data("plotly_click"); req(click)
    out <- grow_df_filtered()[click$pointNumber+1,
                              c("track_id","quad","year","logsize_t0","logsize_t1","id_quad_year")]
    out$track_id <- as.character(out$track_id)
    out$quad <- as.character(out$quad)
    out$logsize_t0 <- as.numeric(out$logsize_t0)
    out$logsize_t1 <- as.numeric(out$logsize_t1)
    clipr::write_clip(out)
    out
  })
  
  output$row <- renderPrint({
    r <- row_r(); req(r)
    cat("track_id:", r$track_id,"\n",
        "quad:", r$quad,"\n",
        "year:", r$year,"\n",
        "log(size_t0):", round(r$logsize_t0,3),"\n",
        "log(size_t1):", round(r$logsize_t1,3),"\n",
        "id_quad_year:", r$id_quad_year,"\n")
  })
  
  observeEvent(list(row_r(), input$show_buffer, input$show_others),{
    r <- row_r(); req(r)
    plots <- plot_before_after_plotly(r$track_id, r$quad, as.numeric(r$year),
                                      input$show_buffer, input$show_others)
    output$map_t0 <- renderPlotly(plots$t0)
    output$map_t1 <- renderPlotly(plots$t1)
  })
  
  observeEvent(event_data("plotly_relayout", source="t0"),{
    rel <- event_data("plotly_relayout", source="t0")
    if(!is.null(rel)) plotlyProxy("map_t1", session) %>% plotlyProxyInvoke("relayout", rel)
  }, ignoreNULL=TRUE)
  observeEvent(event_data("plotly_relayout", source="t1"),{
    rel <- event_data("plotly_relayout", source="t1")
    if(!is.null(rel)) plotlyProxy("map_t0", session) %>% plotlyProxyInvoke("relayout", rel)
  }, ignoreNULL=TRUE)
  
  output$download_excluded <- downloadHandler(
    filename=function(){paste0("excluded_id_quad_years_",Sys.Date(),".csv")},
    content=function(file){
      df_out <- excluded_df() %>% mutate(across(everything(), ~replace_na(.,"")))
      readr::write_csv(df_out, file, na="")
    }
  )
  
  observeEvent(input$upload_excluded,{
    req(input$upload_excluded)
    tryCatch({
      new_df <- readr::read_csv(input$upload_excluded$datapath, show_col_types=FALSE)
      if(!all(c("id_quad_year","comment") %in% names(new_df)))
        stop("CSV must have columns 'id_quad_year' and 'comment'.")
      new_df <- new_df %>% mutate(id_quad_year=as.character(id_quad_year),
                                  comment=replace_na(as.character(comment),"")) %>%
        distinct(id_quad_year, .keep_all=TRUE)
      cur <- excluded_df()
      merged <- cur %>% anti_join(new_df, by="id_quad_year") %>% bind_rows(new_df) %>% arrange(id_quad_year)
      excluded_df(merged)
      showNotification("Excluded list updated from CSV.", type="message")
    }, error=function(e){showNotification(paste("Error reading CSV:", e$message), type="error")})
  })
}

shinyApp(ui, server)
