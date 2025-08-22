# Shiny app data - Anderson 2016 Arizona - Hilaria belangeri

# Author: Niklas Neisse (neisse.n@protonmail.com)
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.08.21

# Version 16: Integrated misfit include/exclude

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

# ---- HELPERS ----
empty_plotly <- function(src) {
  p <- ggplot() + geom_blank() + theme_void()
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

plot_polygons_with_others <- function(df_selected, quad_val, year_val, show_buffer=TRUE, show_others=TRUE, src) {
  if (nrow(df_selected) == 0 || !"geometry" %in% names(df_selected)) {
    p <- ggplot() +
      annotate("text", x = mean(bbox[c("xmin","xmax")]), y = mean(bbox[c("ymin","ymax")]),
               label = "No data available") +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
      theme_bw()
    return(ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout"))
  }
  df_selected <- st_cast(df_selected, "MULTIPOLYGON")
  if (show_buffer) {
    buff <- if_else(st_bbox(df_selected)[3] < 1.1 | st_bbox(df_selected) > 100, 0.05, 5)
    df_buffer <- st_buffer(df_selected, dist = buff)
  }
  df_others <- df %>% filter(quad == quad_val, year == year_val, !(track_id %in% df_selected$track_id))
  if (nrow(df_others) > 0) df_others <- st_cast(df_others, "MULTIPOLYGON")
  p <- ggplot()
  if (show_buffer) p <- p + geom_sf(data = df_buffer, fill = "lightblue", alpha = 0.3)
  if (show_others && nrow(df_others) > 0) p <- p + geom_sf(data = df_others, fill = "lightgrey", alpha = 0.5)
  p <- p + geom_sf(data = df_selected, fill = "dodgerblue") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
    theme_bw()
  ggplotly(p, dynamicTicks = TRUE, source = src) |> event_register("plotly_relayout")
}

plot_before_after_plotly <- function(track_id_val, quad_val, yr_val, show_buffer, show_others) {
  year_t0 <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr_val) %>% select(track_id, geometry)
  year_t1 <- df %>% filter(track_id == track_id_val, quad == quad_val, year == yr_val+1) %>% select(track_id, geometry)
  list(
    t0 = plot_polygons_with_others(year_t0, quad_val, yr_val, show_buffer, show_others, src = "t0"),
    t1 = plot_polygons_with_others(year_t1, quad_val, yr_val+1, show_buffer, show_others, src = "t1")
  )
}

# ---- SHINY UI ----
ui <- fluidPage(
  fluidRow(
    column(width=5,
           tags$div(style="font-size:90%;", class="panel panel-default",
                    tags$div(class="panel-heading",
                             tags$h5(class="panel-title", tags$a(href="#collapseFilters","data-toggle"="collapse", tags$strong("Filters")))),
                    tags$div(id="collapseFilters", class="panel-collapse collapse",
                             tags$div(class="panel-body",
                                      selectizeInput("track_query", "Select track_id_quad:", choices = c("", sort(unique(grow_df$track_id_quad))), selected="", multiple=FALSE),
                                      selectizeInput("site_query", "Filter by site:", choices = c("", sort(unique(grow_df$site))), selected="", multiple=FALSE),
                                      selectizeInput("quad_query", "Filter by quad:", choices = c("", sort(unique(as.character(grow_df$quad)))), selected="", multiple=FALSE),
                                      sliderInput("year_query", "Select year(s):", min=min(grow_df$year), max=max(grow_df$year),
                                                  value=c(min(grow_df$year), max(grow_df$year)), step=1, sep="")
                             )
                    )
           ),
           hr(),
           tags$div(class="panel panel-default",
                    tags$div(class="panel-heading", tags$h5(class="panel-title", tags$strong("Selected point"))),
                    tags$div(class="panel-body", verbatimTextOutput("row"))
           ),
           selectizeInput("misfit_input", "Add to misfits list (id_quad_year):", choices = c("", sort(unique(as.character(grow_df$id_quad_year)))), selected = "", multiple = FALSE),
           actionButton("add_misfit", "Add", class = "btn btn-sm btn-danger"),
           br(), br(),
           tags$div(class="panel panel-default",
                    tags$div(class="panel-heading", tags$h5(class="panel-title", tags$a(href="#collapseMisfits","data-toggle"="collapse", tags$small("Misfits list")))),
                    tags$div(id="collapseMisfits", class="panel-collapse collapse",
                             tags$div(class="panel-body", uiOutput("misfits_list"))
                    )
           ),
           br(),
           downloadButton("download_misfits", "Download misfits CSV", class = "btn btn-sm btn-primary"),
           br(), br(),
           fileInput("upload_misfits", "Upload misfits CSV:", accept = c(".csv"), buttonLabel = "Browse...", placeholder="No file selected")
    ),
    column(width=7,
           tags$div(class = "panel panel-default", style="font-size:85%; margin-bottom:10px;",
                    tags$div(class="panel-heading", tags$h5(class="panel-title", tags$a(href="#collapseJitter","data-toggle"="collapse", tags$strong("Jitter controls")))),
                    tags$div(id="collapseJitter", class="panel-collapse collapse",
                             tags$div(class="panel-body", fluidRow(
                               column(6, sliderInput("jit_w","Width (x-axis):", min=0,max=0.5,value=0,step=0.01)),
                               column(6, sliderInput("jit_h","Height (y-axis):", min=0,max=0.5,value=0,step=0.01))
                             ))
                    )
           ),
           checkboxInput("highlight_point", "Highlight selected point", value=TRUE),
           plotlyOutput("plt", height="400px"),
           br(),
           tags$div(class="panel panel-default", style="padding:10px; margin-top:12px;",
                    fluidRow(
                      column(6, fluidRow(column(6, uiOutput("map_t0_title")), column(6, checkboxInput("show_buffer","Show buffer",value=TRUE))), plotlyOutput("map_t0", height="400px")),
                      column(6, fluidRow(column(6, uiOutput("map_t1_title")), column(6, checkboxInput("show_others","Show other individuals",value=TRUE))), plotlyOutput("map_t1", height="400px"))
                    )
           )
    )
  )
)

# ---- SHINY SERVER (updated scatter plot with fixed tooltip & misfit handling) ----
server <- function(input, output, session) {
  
  # ---- Misfits reactive ----
  misfits_df <- reactiveVal(tibble(id_quad_year=character(), comment=character(), status=character(), group=character()))
  
  observeEvent(input$add_misfit, {
    req(input$misfit_input)
    cur <- misfits_df()
    if (!(input$misfit_input %in% cur$id_quad_year)) {
      misfits_df(bind_rows(tibble(id_quad_year=input$misfit_input, comment="", status="include", group="A"), cur))
    }
  })
  
  # ---- Misfits list UI ----
  output$misfits_list <- renderUI({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return("None")
    safe_id <- function(x) gsub("[^A-Za-z0-9_]","_",x)
    tagList(
      lapply(seq_len(nrow(df_ex)), function(i){
        this_id <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        this_stat <- df_ex$status[i]
        this_group <- df_ex$group[i]
        sid <- safe_id(this_id)
        tags$div(
          style="padding:8px 0; border-bottom:1px solid #ddd;",
          fluidRow(
            column(5, tags$b(this_id),
                   tags$div(style="font-size:90%; color:#555; margin-top:2px;", ifelse(nchar(this_com)>0,this_com,tags$em("(no comment)")))),
            column(7, style="text-align:right;",
                   actionButton(paste0("edit_",sid),"Edit comment", class="btn btn-xs btn-primary", style="margin-right:6px;"),
                   actionButton(paste0("remove_",sid),"Remove", class="btn btn-xs btn-default"))
          ),
          fluidRow(
            column(6, selectInput(paste0("status_",sid),"Status:", choices=c("include","exclude"), selected=this_stat, width="100%")),
            column(6, selectInput(paste0("group_",sid),"Group:", choices=c("A","B","C"), selected=this_group, width="100%"))
          )
        )
      })
    )
  })
  
  # ---- Edit comment modal ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex) == 0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]","_",x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        this_com <- df_ex$comment[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("edit_",sid)]],{
          showModal(modalDialog(title=paste("Edit comment for",this_id),
                                textAreaInput("modal_comment_input","Comment",value=this_com,width="100%",rows=6),
                                footer=tagList(modalButton("Cancel"),actionButton("modal_save_comment","Save",class="btn btn-primary")),
                                easyClose=FALSE))
          observeEvent(input$modal_save_comment,{
            df2 <- misfits_df()
            df2$comment[df2$id_quad_year==this_id] <- input$modal_comment_input %||% ""
            misfits_df(df2)
            removeModal()
          }, once=TRUE, ignoreInit=TRUE)
        }, ignoreInit=TRUE)
      })
    })
  })
  
  # ---- Remove misfit ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex)==0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]","_",x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("remove_",sid)]],{
          df2 <- misfits_df()
          misfits_df(df2[df2$id_quad_year!=this_id,])
        }, ignoreInit=TRUE)
      })
    })
  })
  
  # ---- Sync status & group ----
  observe({
    df_ex <- misfits_df()
    if (nrow(df_ex)==0) return(NULL)
    safe_id <- function(x) gsub("[^A-Za-z0-9_]","_",x)
    lapply(seq_len(nrow(df_ex)), function(i){
      local({
        this_id <- df_ex$id_quad_year[i]
        sid <- safe_id(this_id)
        observeEvent(input[[paste0("status_",sid)]],{
          df2 <- misfits_df()
          df2$status[df2$id_quad_year==this_id] <- input[[paste0("status_",sid)]]
          misfits_df(df2)
        }, ignoreInit=TRUE)
        observeEvent(input[[paste0("group_",sid)]],{
          df2 <- misfits_df()
          df2$group[df2$id_quad_year==this_id] <- input[[paste0("group_",sid)]]
          misfits_df(df2)
        }, ignoreInit=TRUE)
      })
    })
  })
  
  # ---- Filtered data (with misfit exclusion) ----
  grow_df_filtered <- reactive({
    dat <- grow_df
    if (!is.null(input$track_query) && input$track_query!="") dat <- dat %>% filter(track_id_quad==input$track_query)
    if (!is.null(input$site_query) && input$site_query!="") dat <- dat %>% filter(site==input$site_query)
    if (!is.null(input$quad_query) && input$quad_query!="") dat <- dat %>% filter(as.character(quad)==input$quad_query)
    if (!is.null(input$year_query) && length(input$year_query)==2) dat <- dat %>% filter(year>=input$year_query[1], year<=input$year_query[2])
    ex_ids <- misfits_df() %>% filter(status=="exclude") %>% pull(id_quad_year)
    dat %>% filter(!id_quad_year %in% ex_ids)
  })
  
  # ---- Selected point ----
  selected_point <- reactiveVal(NULL)
  observeEvent(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    req(click)
    sel <- grow_df_filtered()[click$pointNumber+1,]
    selected_point(sel)
    clipr::write_clip(sel[,c("track_id","quad","year")])
  })
  
  # ---- Scatter plot with correct tooltip & misfit highlight ----
  output$plt <- renderPlotly({
    suppressWarnings({
      
      dat <- grow_df_filtered()
      validate(need(nrow(dat)>0,"No rows match the current filter."))
      
      misfit_ids <- misfits_df()$id_quad_year
      dat <- dat %>% mutate(is_misfit = id_quad_year %in% misfit_ids)
      
      # Build tooltip
      dat <- dat %>% mutate(tooltip = paste0(
        "id_quad_year: ", id_quad_year,
        "<br>logsize_t0: ", round(logsize_t0,3),
        "<br>logsize_t1: ", round(logsize_t1,3),
        "<br>size_t0: ", round(size_t0,3),
        "<br>size_t1: ", round(size_t1,3),
        ifelse(is_misfit, paste0("<br>Status: ", misfits_df()$status[match(id_quad_year, misfits_df()$id_quad_year)]), "")
      ))
      
      # Plot
      g <- ggplot(dat, aes(x=logsize_t0, y=logsize_t1)) +
        geom_jitter(aes(
          color = ifelse(is_misfit,
                         misfits_df()$status[match(id_quad_year, misfits_df()$id_quad_year)],
                         "#D55E00"),
          text = tooltip
        ),
        width = input$jit_w, height = input$jit_h,
        size = 0.7, alpha = 0.7) +
        scale_color_manual(values=c("include"="#1b9e77","exclude"="grey","#D55E00"="#D55E00")) +
        theme_bw(base_size=11) +
        theme(panel.grid.major=element_line(color="grey85"), panel.grid.minor=element_line(color="grey92"), legend.position="none") +
        labs(x="Log(Size t0)",y="Log(Size t1)")
      
      # Highlight selected point
      sel <- selected_point()
      if (!is.null(sel) && nrow(sel)>0 && isTRUE(input$highlight_point)) {
        g <- g + geom_point(data=sel, aes(x=logsize_t0, y=logsize_t1), color="black", size=2, shape=17)
      }
      
      ggplotly(g, tooltip="text") |> event_register("plotly_click")
    })
  })
  
  # ---- Selected row output ----
  output$row <- renderPrint({
    sel <- selected_point()
    if (is.null(sel)) {
      cat("No point selected")
    } else {
      sel_simple <- sel %>%
        st_set_geometry(NULL) %>%   # remove geometry
        select(id_quad_year, logsize_t0, logsize_t1, size_t0, size_t1, survives)  
      print(sel_simple)
    }
  })
  
  # ---- Maps ----
  output$map_t0_title <- renderUI({ tags$strong("Year t0") })
  output$map_t1_title <- renderUI({ tags$strong("Year t1") })
  
  output$map_t0 <- renderPlotly({
    sel <- selected_point()
    if (is.null(sel)) return(empty_plotly("t0"))
    plot_polygons_with_others(df %>% filter(track_id==sel$track_id, quad==sel$quad, year==sel$year),
                              sel$quad, sel$year, show_buffer=input$show_buffer, show_others=input$show_others, src="t0")
  })
  
  output$map_t1 <- renderPlotly({
    sel <- selected_point()
    if (is.null(sel)) return(empty_plotly("t1"))
    plot_polygons_with_others(df %>% filter(track_id==sel$track_id, quad==sel$quad, year==sel$year+1),
                              sel$quad, sel$year+1, show_buffer=input$show_buffer, show_others=input$show_others, src="t1")
  })
  
  # ---- Download CSV ----
  output$download_misfits <- downloadHandler(
    filename = function() { paste0("misfits_",Sys.Date(),".csv") },
    content = function(file) { write_csv(misfits_df(), file) }
  )
  
  # ---- Upload CSV ----
  observeEvent(input$upload_misfits,{
    req(input$upload_misfits)
    tryCatch({
      new_df <- readr::read_csv(input$upload_misfits$datapath, show_col_types = FALSE)
      stopifnot(all(c("id_quad_year","status") %in% names(new_df)))
      old_df <- misfits_df()
      merged <- bind_rows(new_df, old_df %>% filter(!id_quad_year %in% new_df$id_quad_year))
      misfits_df(merged)
    }, error = function(e){showNotification("Error: CSV must have columns 'id_quad_year' and 'status'.", type="error")})
  })
  
}

# ---- RUN APP ----
shinyApp(ui, server)
