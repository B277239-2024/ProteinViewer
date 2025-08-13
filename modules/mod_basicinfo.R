# modules/mod_basicinfo.R

mod_basicinfo_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("seq_source"), "Sequence Source",
                 choices = c("API" = "api", "Upload File" = "upload"),
                 selected = "api", inline = TRUE),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'api'", ns("seq_source")),
      textInput(ns("pid"), "Enter UniProt ID", value = "Q9Y6K1")
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("seq_source")),
      tagList(
        fileInput(ns("fasta_upload"), "Upload FASTA File", accept = c(".fasta", ".fa"),
                  placeholder = "Upload a FASTA file"),
        helpText("FASTA File provides the sequence information. "),
        hr(),
        fileInput(ns("gnomad_upload"), "Upload gnomAD-like CSV", accept = c(".csv"),
                  placeholder = "Optional variant CSV with HGVSp column"),
        helpText("CSV must contain HGVSp column like 'p.Arg12Cys'.", 
                 br(),
                 "CSV File provides the variants information. "),
        hr()
      )
    ),
    actionButton(ns("fetch"), "Fetch info"),
    hr(),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'api'", ns("seq_source")),
      selectInput(ns("selected_transcript"), "Select Transcript", choices = NULL),
      helpText("The default transcript is the first canonical one.", 
               br(),
               "You can choose a different transcript id from the dropdown.")
    ),
    tags$hr(),
    checkboxInput(ns("use_custom_domain"), "Use your own domain info", value = FALSE),
    conditionalPanel(
      condition = sprintf("input['%s'] == true", ns("use_custom_domain")),
      tagList(
        rhandsontable::rHandsontableOutput(ns("custom_domain_table")),
        fluidRow(
          column(6, actionButton(ns("add_row"), "Add Row", icon = icon("plus"))),
          column(6, actionButton(ns("clear_table"), "Clear Table", icon = icon("trash")))
        ),
        helpText("Enter columns: Description (text), Start (integer), End (integer)",
                 br(),
                 "Rows with missing values will be ignored. ")
      )
    ),
    helpText("Domain info from API is retrieved from UniProt.",br(),
             "For more accurate domain definitions, visit:",br(),
             "https://alphafold.ebi.ac.uk/entry/your_uniprot_id", br(),
             "and refer to TED Domains and Predicted Aligned Error (PAE).", br(),
             "You may enter domains manually if needed.")
  )
}

mod_basicinfo_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    canonical_tx <- reactiveVal(NULL)
    
    domain_table_data <- reactiveVal(
      data.frame(
        description = as.character(c("", "")),
        start = as.integer(c(NA, NA)),
        end = as.integer(c(NA, NA)),
        stringsAsFactors = FALSE
      )
    )
    
    observeEvent(input$add_row, {
      df <- domain_table_data()
      df <- rbind(df, data.frame(
        description = "", start = NA, end = NA, stringsAsFactors = FALSE
      ))
      domain_table_data(df)
    })
    
    observeEvent(input$clear_table, {
      domain_table_data(
        data.frame(
          description = as.character(c("", "")),
          start = as.integer(c(NA, NA)),
          end = as.integer(c(NA, NA)),
          stringsAsFactors = FALSE
        )
      )
    })
    
    output$custom_domain_table <- rhandsontable::renderRHandsontable({
      req(input$use_custom_domain)
      rhandsontable::rhandsontable(domain_table_data(), rowHeaders = NULL) %>%
        rhandsontable::hot_cols(colWidths = 110)
    })
    
    # Fetch UniProt JSON
    protein_data <- eventReactive(input$fetch, {
      if (input$seq_source == "upload") {
        # No need to fetch from API
        return(NULL)
      }
      req(input$pid)
      url <- paste0("https://rest.uniprot.org/uniprotkb/", input$pid, ".json")
      res <- httr::GET(url)
      if (httr::status_code(res) == 200) {
        httr::content(res, as = "parsed", type = "application/json")
      } else {
        NULL
      }
    })
    
    # Read FASTA
    uploaded_fasta_seq <- reactive({
      req(input$fasta_upload)
      tryCatch({
        fasta_lines <- readLines(input$fasta_upload$datapath)
        seq_lines <- fasta_lines[!grepl("^>", fasta_lines)]
        paste(seq_lines, collapse = "")
      }, error = function(e) {
        showNotification("Failed to read FASTA file", type = "error")
        return(NULL)
      })
    })
    
    # Unified sequence output
    current_sequence <- reactive({
      if (input$seq_source == "upload") {
        return(uploaded_fasta_seq())
      } else {
        req(protein_data())
        protein_data()$sequence$value
      }
    })
    
    # Transcripts for API mode
    get_transcripts_by_uniprot <- function(uniprot_id, canonical_tx) {
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      result <- biomaRt::getBM(
        attributes = c("uniprotswissprot", "ensembl_transcript_id", "external_gene_name", "transcript_biotype", "transcript_is_canonical"),
        filters = "uniprotswissprot",
        values = uniprot_id,
        mart = ensembl
      )
      protein_tx <- subset(result, transcript_biotype == "protein_coding")
      protein_tx <- unique(protein_tx)
      canonical_id <- protein_tx$ensembl_transcript_id[which(protein_tx$transcript_is_canonical == 1)][1]
      if (!is.null(canonical_id)) canonical_tx(canonical_id)
      return(protein_tx)
    }
    
    transcript_table <- reactive({
      req(input$pid)
      get_transcripts_by_uniprot(input$pid, canonical_tx)
    })
    
    observeEvent(transcript_table(), {
      tx_df <- transcript_table()
      if (nrow(tx_df) == 0) return()
      choices <- setNames(tx_df$ensembl_transcript_id, 
                          paste0(tx_df$ensembl_transcript_id, " (", tx_df$external_gene_name, ")"))
      updateSelectInput(session, "selected_transcript",
                        choices = choices,
                        selected = canonical_tx())
    })
    
    custom_domain_df <- reactive({
      if (!isTRUE(input$use_custom_domain)) return(NULL)
      table <- input$custom_domain_table
      if (is.null(table)) return(NULL)
      df <- tryCatch({
        df <- as.data.frame(rhandsontable::hot_to_r(table))
        domain_table_data(df)
        df
      }, error = function(e) {
        return(NULL)
      })
      if (!all(c("description", "start", "end") %in% colnames(df))) return(NULL)
      df <- df %>% dplyr::filter(!is.na(start) & !is.na(end) & nzchar(description))
      return(df)
    })
    
    domain_df <- reactive({
      custom <- custom_domain_df()
      if (!is.null(custom) && nrow(custom) > 0) {
        return(custom)
      }
      
      if (input$seq_source == "api") {
        data <- protein_data()
        if (!is.null(data$features)) {
          domains <- data$features[sapply(data$features, function(x) x$type == "Domain")]
          if (length(domains) > 0) {
            return(data.frame(
              start = as.integer(sapply(domains, function(x) x$location$start$value)),
              end   = as.integer(sapply(domains, function(x) x$location$end$value)),
              description = sapply(domains, function(x) x$description),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      return(NULL)
    })
    
    return(list(
      pid = reactive({ input$pid }),
      seq_source = reactive({ input$seq_source }),
      fasta_upload = reactive({ input$fasta_upload }),
      gnomad_upload = reactive({ input$gnomad_upload }),
      selected_transcript = reactive({ input$selected_transcript }),
      canonical_tx = canonical_tx,
      protein_data = protein_data,
      current_sequence = current_sequence,
      transcript_table = transcript_table,
      use_custom_domain = reactive({ input$use_custom_domain }),
      custom_domain_df = custom_domain_df,
      domain_df = domain_df
    ))
  })
}