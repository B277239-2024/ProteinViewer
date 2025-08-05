# modules/mod_comparison.R

mod_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Protein Comparison"),
    helpText("Enter multiple UniProt IDs for comparison. Domain is always shown. Other features are optional."),
    rhandsontable::rHandsontableOutput(ns("uid_table")),
    actionButton(ns("add_uid"), "Add Row", icon = icon("plus")),
    checkboxGroupInput(ns("features"), "Select features to show:",
                       choices = c("Missense", "PTM", "AlphaMissense"),
                       selected = NULL,
                       inline = TRUE),
    actionButton(ns("run"), "Run Comparison", class = "btn-primary")
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
      rhandsontable::rhandsontable(uid_table_data(), rowHeaders = NULL)
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
            tx_id <- tx_df$ensembl_transcript_id[tx_df$transcript_is_canonical == 1][1]
            if (!is.null(tx_id)) get_gnomad_variants_for_transcript(tx_id) else NULL
          } else NULL,
          am_df = if ("AlphaMissense" %in% input$features) get_alphamissense_data(uid) else NULL
        )
      }) |> purrr::compact()
    })
    
    output$comparison_plot_ui <- renderUI({
      req(comparison_data())
      plotOutput(ns("comparison_plot"), height = paste0(300 * length(comparison_data()), "px"))
    })
    
    output$comparison_plot <- renderPlot({
      req(comparison_data())
      plots <- lapply(comparison_data(), function(d) {
        render_1d_plot_comparison(
          seq_val = d$sequence,
          domain_df = d$domain_df,
          missense_df = d$missense_df,
          ptm_df = d$ptm_df,
          am_df = d$am_df,
          title = d$uid
        )
      })
      patchwork::wrap_plots(plots, ncol = 1)
    })
  })
}