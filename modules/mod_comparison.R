# modules/mod_comparison.R

mod_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Protein Comparison"),
    helpText("Enter multiple UniProt IDs for comparison."),
    helpText("Domain is always shown. Other features are optional."),
    # rhandsontable::rHandsontableOutput(ns("uid_table")),
    div(style = "width: 400px; overflow-x: hidden;",
        rhandsontable::rHandsontableOutput(ns("uid_table"))
    ),
    br(),
    actionButton(ns("add_uid"), "Add Row", icon = icon("plus")),
    hr(),
    checkboxGroupInput(ns("features"), "Select features to show:",
                       choices = c("Missense", "PTM", "AlphaMissense"),
                       selected = NULL,
                       inline = TRUE),
    checkboxInput(ns("fixed_scale"), "Use fixed x-axis scale (shared across proteins)", value = FALSE),
    actionButton(ns("run"), "Run Comparison", class = "btn-primary"),
    hr(),
    downloadButton(ns("download_plot"), "Download Comparison Plot (PNG)", icon = icon("download")),
    br(),
    uiOutput(ns("comparison_legend"))
  )
}

mod_comparison_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    uid_table_data <- reactiveVal(
      data.frame(UniProtID = "", stringsAsFactors = FALSE)
    )
    
    observeEvent(input$add_uid, {
      df <- uid_table_data()
      df <- rbind(df, data.frame(UniProtID = "", stringsAsFactors = FALSE))
      uid_table_data(df)
    })
    
    output$uid_table <- rhandsontable::renderRHandsontable({
      rhandsontable::rhandsontable(uid_table_data(), rowHeaders = NULL) %>%
        rhandsontable::hot_cols(colWidths = 200)
    })
    
    observeEvent(input$uid_table, {
      tbl <- tryCatch(rhandsontable::hot_to_r(input$uid_table), error = function(e) NULL)
      if (!is.null(tbl)) uid_table_data(tbl)
    })
    
    comparison_data <- eventReactive(input$run, {
      df <- uid_table_data() # users enter Uniprot ids
      uids <- na.omit(df$UniProtID)
      if (length(uids) == 0) return(NULL)
      
      lapply(uids, function(uid) {
        protein_data <- get_uniprot_protein_data(uid)
        if (is.null(protein_data)) return(NULL)
        
        list(
          uid = uid,
          protein_data = protein_data,
          sequence = get_protein_sequence(protein_data),
          domain_df = extract_uniprot_domains(protein_data),
          ptm_df = if ("PTM" %in% input$features) get_ptm_data(protein_data) else NULL,
          missense_df = if ("Missense" %in% input$features) {
            tx_df <- tryCatch(get_transcripts_by_uniprot(uid), error = function(e) NULL)
            tx_id <- tx_df$ensembl_transcript_id[which(tx_df$transcript_is_canonical == 1)][1]
            if (!is.null(tx_id) && !is.na(tx_id)) {
              get_gnomad_variants_for_transcript(tx_id)
              # cat("Comparison transcript:", tx_id)
            } else {
              NULL
            }
          } else NULL,
          am_df = if ("AlphaMissense" %in% input$features) get_alphamissense_data(uid) else NULL
        )
      }) |> purrr::compact()
    })
    
    domain_color_map <- reactive({
      cds <- comparison_data()
      all_domains <- do.call(rbind, lapply(cds, function(d) d$domain_df))
      unique_names <- unique(all_domains$description)
      setNames(scales::hue_pal()(length(unique_names)), unique_names)
    })
    
    text_plot <- function(label) {
      ggplot() +
        annotate("text", x = 1, y = 1, label = label, hjust = 1, size = 4, fontface = "bold") +
        theme_void() +
        xlim(0, 1.1)
    }
    
    output$comparison_plot_ui <- renderUI({
      req(comparison_data())
      plotOutput(ns("comparison_plot"), height = paste0(200 * length(comparison_data()), "px"))
    })
    
    output$comparison_plot <- renderPlot({
      req(comparison_data())
      data <- comparison_data()
      max_len <- max(sapply(data, function(d) nchar(d$sequence)))
      
      row_plots <- lapply(data, function(d) {
        label <- text_plot(d$uid)
        main_plot <- render_1d_plot_comparison(
          seq_val = d$sequence,
          domain_df = d$domain_df,
          missense_df = d$missense_df,
          ptm_df = d$ptm_df,
          am_df = d$am_df,
          title = NULL,
          fixed_x_max = if (input$fixed_scale) max_len else NULL,
          domain_color_map = domain_color_map()
        )
        label + main_plot + patchwork::plot_layout(ncol = 2, widths = c(0.06, 0.94))
      })
      patchwork::wrap_plots(row_plots, ncol = 1)
    })
    
    output$comparison_legend <- renderUI({
      tagList(
        tags$div(style = "padding: 10px; font-size: 13px;",
                 tags$h5("Legend"),
                 tags$ul(style = "list-style-type: none; padding-left: 0; margin-bottom: 0;",
                         
                         # Protein bar
                         tags$li(
                           tags$span(style = "display: inline-block; width: 12px; height: 12px;
                             background-color: #d3d3d3; border: 1px solid #ccc;
                             margin-right: 6px;"),
                           "Grey bar – Protein length"
                         ),
                         
                         # Domain block
                         tags$li(
                           tags$span(style = "display: inline-block; width: 12px; height: 12px;
                             background-color: #fca6a6; border: 1px solid #000;
                             margin-right: 6px;"),
                           "Pink block – Domain"
                         ),
                         
                         # Missense line
                         tags$li(
                           tags$span(style = "display: inline-block; width: 12px; height: 2px;
                             background-color: slateblue; margin-right: 6px; vertical-align: middle;"),
                           "Blue line – Missense variant"
                         ),
                         
                         # PTM
                         tags$li("Colored dots – PTM sites"),
                         
                         tags$li(
                           tags$span(style = "display: inline-block; width: 10px; height: 10px; background-color: #1f77b4; margin-right: 6px; border-radius: 50%; display: inline-block;"), 
                           "Phosphorylation"
                         ),
                         tags$li(
                           tags$span(style = "display: inline-block; width: 10px; height: 10px; background-color: #ff7f0e; margin-right: 6px; border-radius: 50%; display: inline-block;"), 
                           "Acetylation"
                         ),
                         tags$li(
                           tags$span(style = "display: inline-block; width: 10px; height: 10px; background-color: #2ca02c; margin-right: 6px; border-radius: 50%; display: inline-block;"), 
                           "Succinylation"
                         ),
                         tags$li(
                           tags$span(style = "display: inline-block; width: 10px; height: 10px; background-color: #d62728; margin-right: 6px; border-radius: 50%; display: inline-block;"), 
                           "Methylation"
                         ),
                         tags$li(
                           tags$span(style = "display: inline-block; width: 10px; height: 10px; background-color: #9467bd; margin-right: 6px; border-radius: 50%; display: inline-block;"), 
                           "Other"
                         ),
                         
                         # AlphaMissense
                         tags$li(
                           tags$span(style = "display: inline-block; width: 12px; height: 12px;
                             background-color: darkgreen; margin-right: 6px;"),
                           "Green bar – AlphaMissense average score"
                         )
                 )
        )
      )
    })
    
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("comparison_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(comparison_data())
        data <- comparison_data()
        max_len <- max(sapply(data, function(d) nchar(d$sequence)))
        
        row_plots <- lapply(data, function(d) {
          label <- text_plot(d$uid)
          main_plot <- render_1d_plot_comparison(
            seq_val = d$sequence,
            domain_df = d$domain_df,
            missense_df = d$missense_df,
            ptm_df = d$ptm_df,
            am_df = d$am_df,
            title = NULL,
            fixed_x_max = if (input$fixed_scale) max_len else NULL,
            domain_color_map = domain_color_map()
          )
          label + main_plot + patchwork::plot_layout(ncol = 2, widths = c(0.08, 0.92))
        })
        
        final_plot <- patchwork::wrap_plots(row_plots, ncol = 1)
        
        ggsave(file, final_plot, width = 12, height = 3 * length(row_plots), dpi = 300, bg = "white")
      }
    )
  })
}