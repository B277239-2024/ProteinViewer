# modules/mod_gnomad.R

mod_gnomad_server <- function(id, seq_source, selected_transcript, gnomad_upload) {
  moduleServer(id, function(input, output, session) {
    
    # --- 数据获取逻辑 ---
    
    fetch_gnomad_by_transcript <- function(transcript_id) {
      url <- "https://gnomad.broadinstitute.org/api"
      query <- '
      query TranscriptVariants($transcriptID: String!) {
        transcript(transcript_id: $transcriptID, reference_genome: GRCh38) {
          variants(dataset: gnomad_r4) {
              variant_id
              pos
              consequence
              transcript_id
              hgvsp
              exome { af ac an }
          }
        }
      }'
      body <- list(
        query = query,
        variables = list(transcriptID = transcript_id)
      )
      res <- httr::POST(url, body = body, encode = "json", httr::content_type_json())
      if (httr::status_code(res) == 200) {
        httr::content(res, as = "parsed", simplifyVector = TRUE)
      } else {
        warning(paste("gnomAD API request failed:", httr::status_code(res)))
        NULL
      }
    }
    
    uploaded_gnomad_df <- reactive({
      req(gnomad_upload())
      tryCatch({
        df_raw <- read.csv(gnomad_upload()$datapath, stringsAsFactors = FALSE, check.names = FALSE)
        col_map <- list(
          Variant_ID = "gnomAD ID", 
          Position = "Position",
          Consequence = "VEP Annotation",
          HGVSp = "HGVS Consequence",
          AF = "Allele Frequency",
          AC = "Allele Count",
          AN = "Allele Number"
        )
        df_cols <- lapply(names(col_map), function(newname) {
          oldname <- col_map[[newname]]
          if (oldname %in% colnames(df_raw)) {
            return(df_raw[[oldname]])
          } else {
            return(rep(NA, nrow(df_raw)))
          }
        })
        names(df_cols) <- names(col_map)
        df <- as.data.frame(df_cols, stringsAsFactors = FALSE)
        
        df$AA_Position <- sapply(df$HGVSp, function(hgvsp) {
          if (!is.na(hgvsp) && grepl("^p\\.", hgvsp)) {
            matches <- regmatches(hgvsp, regexec("\\d+", hgvsp))[[1]]
            if (length(matches) > 0) as.integer(matches[1]) else NA
          } else {
            NA
          }
        })
        df
      }, error = function(e) {
        showNotification("Failed to read gnomAD CSV file", type = "error")
        return(NULL)
      })
    })
    
    gnomad_df <- reactive({
      req(seq_source())
      if (seq_source() == "upload") {
        uploaded_gnomad_df()
      } else {
        req(selected_transcript())
        gnomad_data <- fetch_gnomad_by_transcript(selected_transcript())
        vars <- tryCatch(gnomad_data$data$transcript$variants, error = function(e) NULL)
        if (is.null(vars)) return(NULL)
        
        df <- data.frame(
          Variant_ID = vars$variant_id,
          Position = vars$pos,
          Consequence = vars$consequence,
          HGVSp = vars$hgvsp,
          Transcript_ID = vars$transcript_id,
          AF = vars$exome$af,
          AC = vars$exome$ac,
          AN = vars$exome$an
        )
        df$AA_Position <- sapply(df$HGVSp, function(hgvsp) {
          if (is.na(hgvsp) || hgvsp == "") return(NA)
          matches <- regmatches(hgvsp, regexec("\\d+", hgvsp))[[1]]
          if (length(matches) > 0) as.integer(matches[1]) else NA
        })
        df
      }
    })
    
    # --- summary UI ---
    summary_ui <- renderUI({
      df <- gnomad_df()
      if (is.null(df) || nrow(df) == 0) return("No variants found")
      consequence_counts <- sort(table(df$Consequence), decreasing = TRUE)
      consequence_list <- lapply(names(consequence_counts), function(name) {
        tags$li(paste(name, ":", consequence_counts[[name]]))
      })
      total_variants <- nrow(df)
      tags$div(
        tags$p(paste("Total variants: ", total_variants)),
        tags$p("Consequence types:"),
        tags$ul(consequence_list)
      )
    })
    
    # --- table output ---
    table_ui <- DT::renderDT({
      df <- gnomad_df()
      if (is.null(df) || nrow(df) == 0) return(NULL)
      DT::datatable(
        df,
        rownames = FALSE,
        filter = "top",
        extensions = 'Buttons',
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50),
          dom = 'lBfrtip',
          buttons = c('excel', 'pdf', 'csv'),
          scrollX = TRUE,
          columnDefs = list(list(width = '120px', targets = "_all"))
        ),
        colnames = c(
          "Variant ID" = "Variant_ID",
          "Genomic Position" = "Position",
          "Consequence" = "Consequence",
          "Protein Change" = "HGVSp",
          "Residue Position" = "AA_Position",
          "Allele Freq" = "AF",
          "Allele Count" = "AC",
          "Allele Number" = "AN"
        ),
        escape = FALSE
      )
    }, server = FALSE)
    
    return(list(
      gnomad_df = gnomad_df,
      summary_ui = summary_ui,
      table_ui = table_ui
    ))
  })
}