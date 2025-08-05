# modules/mod_alphamissense.R

mod_alphamissense_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("am_source"), "AlphaMissense Source",
                 choices = c("API" = "api", "Upload File" = "upload"),
                 selected = "upload", inline = TRUE),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("am_source")),
      fileInput(ns("alphamissense_upload"), "Upload AlphaMissense CSV", accept = ".csv")
    ),
    checkboxInput(ns("add_alphamissense"), "Add AlphaMissense Layer to 1D plot", value = FALSE),
    checkboxGroupInput(ns("am_prediction_filter"), "Prediction Category",
                       choices = c("likely_benign", "ambiguous", "likely_pathogenic"),
                       selected = NULL),
    helpText("Choose prediction types to display on heatmap."),
    uiOutput(ns("download_ui")),
    br(),
    uiOutput(ns("download_data_ui"))
  )
}

mod_alphamissense_server <- function(id, pid) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$download_ui <- renderUI({
      df <- alphamissense_df()
      if (!is.null(df) && nrow(df) > 0) {
        downloadButton(ns("download_am_heatmap"), "Download Heatmap (PNG)", icon = icon("download"))
      } else {
        return(NULL)
      }
    })
    
    output$download_data_ui <- renderUI({
      is_api <- input$am_source == "api"
      df <- alphamissense_df()
      
      if (is_api && !is.null(df) && nrow(df) > 0) {
        downloadButton(ns("download_am_data"), "Download AlphaMissense CSV", icon = icon("file-download"))
      } else {
        return(NULL)
      }
    })
    
    fetch_alphafold_am_url <- function(uniprot_id) {
      url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot_id)
      res <- httr::GET(url)
      if (httr::status_code(res) == 200) {
        json <- httr::content(res, as = "parsed", type = "application/json")
        if (!is.null(json[[1]]$amAnnotationsUrl)) {
          return(json[[1]]$amAnnotationsUrl)
        }
      }
      return(NULL)
    }
    
    alphamissense_df <- reactive({
      req(input$am_source)
      file_path <- NULL
      
      if (input$am_source == "upload") {
        req(input$alphamissense_upload)
        file_path <- input$alphamissense_upload$datapath
      } else if (input$am_source == "api"){
        req(pid())
        url <- fetch_alphafold_am_url(pid())
        if (is.null(url)) {
          showNotification("AlphaMissense annotations not found for this UniProt ID from AFDB, please try to upload the file.", type = "warning")
          return(NULL)
        }
        file_path <- tempfile(fileext = ".csv")
        tryCatch({
          download.file(url, destfile = file_path, quiet = TRUE)
        }, error = function(e) {
          showNotification("Failed to download AlphaMissense CSV from API.", type = "error")
          return(NULL)
        })
      }
      
      df <- readr::read_csv(file_path, show_col_types = FALSE)
      
      df <- df %>%
        dplyr::mutate(
          aa_from = stringr::str_sub(protein_variant, 1, 1),
          position = as.integer(stringr::str_extract(protein_variant, "\\d+")),
          aa_to = stringr::str_sub(protein_variant, -1, -1),
          aa_to = factor(aa_to, levels = rev(c(
            "A", "C", "D", "E", "F", "G", "H", "I", "K", 
            "L", "M", "N", "P", "Q", "R", "S", "T", "V", 
            "W", "Y"
          ))),
          am_class = dplyr::recode(am_class,
                                   "Amb" = "ambiguous",
                                   "LBen" = "likely_benign",
                                   "LPath" = "likely_pathogenic")
        ) %>%
        dplyr::rename(score = am_pathogenicity)
      
      df
    })
    
    alphamissense_enabled <- reactive({
      input$add_alphamissense && !is.null(alphamissense_df()) && nrow(alphamissense_df()) > 0
    })
    
    output$download_am_heatmap <- downloadHandler(
      filename = function() {
        paste0("alphamissense_heatmap_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(am_res$alphamissense_df())
        df <- alphamissense_df()
        
        filter_vals <- input$am_prediction_filter
        
        if (is.null(filter_vals) || length(filter_vals) == 0) {
          filter_vals <- unique(df$am_class)
        }
        
        df <- df[df$am_class %in% filter_vals, ]
        req(nrow(df) > 0)
        
        ref_df <- df %>%
          dplyr::select(position, aa_from) %>%
          dplyr::distinct() %>%
          dplyr::mutate(aa_from = factor(aa_from, levels = levels(df$aa_to)))
        
        p <- ggplot() +
          geom_tile(
            data = df,
            aes(x = position, y = aa_to, fill = score),
            color = NA
          ) +
          geom_tile(
            data = ref_df,
            aes(x = position, y = aa_from),
            fill = "black", width = 1, height = 1, alpha = 0.9
          ) +
          scale_fill_gradientn(
            colors = c("blue", "white", "red"),
            values = scales::rescale(c(0, 0.5, 1)),
            limits = c(0, 1),
            name = "Score"
          ) +
          labs(
            title = "AlphaMissense Heatmap",
            x = "Residue Position",
            y = "Alternative Amino Acid"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 8),
            panel.grid = element_blank()
          )
        
        ggsave(file, plot = p, width = 12, height = 6, dpi = 300, bg = "white")
      }
    )
    
    output$download_am_data <- downloadHandler(
      filename = function() {
        paste0("AlphaMissense_API_", Sys.Date(), ".csv")
      },
      content = function(file) {
        df <- alphamissense_df()
        if (is.null(df) || nrow(df) == 0) {
          showNotification("No AlphaMissense API data to download.", type = "error")
          return(NULL)
        }
        
        readr::write_csv(df, file)
      }
    )
    
    return(list(
      alphamissense_df = alphamissense_df,
      enabled = alphamissense_enabled,
      prediction_filter = reactive({ input$am_prediction_filter }),
      am_source = reactive({ input$am_source }),
      upload_info = reactive({ input$alphamissense_upload })
    ))
  })
}