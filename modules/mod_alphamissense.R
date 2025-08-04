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
    actionButton(ns("add_alphamissense"), "Add AlphaMissense Layer", icon = icon("plus")),
    checkboxGroupInput(ns("am_prediction_filter"), "Prediction Category",
                       choices = c("likely_benign", "ambiguous", "likely_pathogenic"),
                       selected = NULL),
    helpText("Choose prediction types to display on heatmap.")
  )
}

mod_alphamissense_server <- function(id, pid) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    alphamissense_enabled <- reactiveVal(FALSE)
    
    observe({
      df <- tryCatch(alphamissense_df(), error = function(e) NULL)
      ready <- !is.null(df) && nrow(df) > 0
      shinyjs::toggleState("add_alphamissense", condition = ready)
    })
    
    observeEvent(input$add_alphamissense, {
      req(alphamissense_df())  
      alphamissense_enabled(TRUE)
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
    
    return(list(
      alphamissense_df = alphamissense_df,
      enabled = alphamissense_enabled,
      prediction_filter = reactive({ input$am_prediction_filter }),
      am_source = reactive({ input$am_source }),
      upload_info = reactive({ input$alphamissense_upload })
    ))
  })
}