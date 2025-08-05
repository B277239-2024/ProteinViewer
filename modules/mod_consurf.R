# modules/mod_consurf.R

mod_consurf_1d_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("consurf_txt"), "Upload ConSurf TXT File", accept = c(".txt")),
    checkboxInput(ns("add_consurf"), "Add ConSurf Layer to 1D plot", value = FALSE),
    downloadButton(ns("download_consurf_plot"), "Download ConSurf Plot (PNG)", icon = icon("download"))
  )
}

mod_consurf_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    read_consurf_txt <- function(file_path) {
      df <- read.delim(file_path, skip = 27, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      colnames(df) <- trimws(colnames(df))
      
      df_clean <- df %>%
        dplyr::select(Position = POS, Score = SCORE, Grade = COLOR) %>%
        dplyr::mutate(
          Position = as.integer(Position),
          Score = as.numeric(Score),
          Grade = as.integer(Grade)
        )
      return(df_clean)
    }
    
    consurf_txt_df <- reactive({
      req(input$consurf_txt)
      tryCatch({
        read_consurf_txt(input$consurf_txt$datapath)
      }, error = function(e) {
        showNotification("Failed to read ConSurf TXT file", type = "error")
        NULL
      })
    })
    
    consurf_enabled <- reactive({
      input$add_consurf && !is.null(consurf_txt_df()) && nrow(consurf_txt_df()) > 0
    })
    
    consurf_cols <- c(
      "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
      "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
    )
    
    output$plot_ui <- renderUI({
      req(consurf_txt_df())
      tagList(
        plotOutput(ns("txt_plot"), height = "400px")
      )
    })
    
    output$txt_plot <- renderPlot({
      req(consurf_txt_df())
      df <- consurf_txt_df()
      
      ggplot(df, aes(x = Position, y = Score, fill = factor(Grade))) +
        geom_col(width = 1) +
        scale_fill_manual(values = consurf_cols, name = "ConSurf Grade") +
        theme_minimal() +
        labs(title = "ConSurf Conservation Score", x = "Residue Position", y = "Score")
    })
    
    output$download_consurf_plot <- downloadHandler(
      filename = function() {
        paste0("consurf_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(consurf_txt_df())
        df <- consurf_txt_df()
        
        consurf_cols <- c(
          "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
          "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
        )
        
        p <- ggplot(df, aes(x = Position, y = Score, fill = factor(Grade))) +
          geom_col(width = 1) +
          scale_fill_manual(values = consurf_cols, name = "ConSurf Grade") +
          theme_minimal() +
          labs(title = "ConSurf Conservation Score", x = "Residue Position", y = "Score")
        
        ggplot2::ggsave(file, plot = p, width = 10, height = 4, dpi = 300, bg = "white")
      }
    )
    
    plot_ui_element <- tagList(
      plotOutput(ns("txt_plot"), height = "400px")
    )
    
    consurf_pdb_file <- reactive({
      req(input$consurf_pdb_upload)
      tryCatch({
        path <- input$consurf_pdb_upload$datapath
        list(
          bio3d = bio3d::read.pdb(path),
          text = readChar(path, file.info(path)$size)
        )
      }, error = function(e) {
        showNotification("Failed to read ConSurf PDB file", type = "error")
        NULL
      })
    })
    
    show_consurf_3d <- reactive({
      isTRUE(input$toggle_consurf_3d)
    })
    
    return(list(
      consurf_df = consurf_txt_df,
      enabled = consurf_enabled,
      plot_ui = plot_ui_element
    ))
  })
}