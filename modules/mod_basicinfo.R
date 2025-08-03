# modules/mod_basicinfo.R

mod_basicinfo_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("seq_source"), "Sequence Source",
                 choices = c("API" = "api", "Upload File" = "upload"),
                 selected = "api", inline = TRUE),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'api'", ns("seq_source")),
      textInput(ns("pid"), "Enter UniProt ID")
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("seq_source")),
      fileInput(ns("fasta_upload"), "Upload FASTA File", accept = c(".fasta", ".fa"),
                placeholder = "Upload a FASTA file")
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("seq_source")),
      fileInput(ns("gnomad_upload"), "Upload gnomAD-like CSV",
                accept = c(".csv"),
                placeholder = "Optional variant CSV with HGVSp column"),
      helpText("Optional. Must contain HGVSp column like 'p.Arg12Cys'.")
    ),
    actionButton(ns("fetch"), "Fetch info"),
    tags$div(style = "margin-top: 10px;"), 
    conditionalPanel(
      condition = sprintf("input['%s'] == 'api'", ns("seq_source")),
      selectInput(ns("selected_transcript"), "Select Transcript", choices = NULL)
    )
  )
}

mod_basicinfo_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    canonical_tx <- reactiveVal(NULL)
    
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
    
    return(list(
      pid = reactive({ input$pid }),
      seq_source = reactive({ input$seq_source }),
      fasta_upload = reactive({ input$fasta_upload }),
      gnomad_upload = reactive({ input$gnomad_upload }),
      selected_transcript = reactive({ input$selected_transcript }),
      canonical_tx = canonical_tx,
      protein_data = protein_data,
      current_sequence = current_sequence,
      transcript_table = transcript_table
    ))
  })
}