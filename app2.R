library(shiny)
library(httr)
library(httr2)
library(jsonlite)
library(DT)
library(r3dmol)
library(shinycssloaders)
library(bio3d)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(plotly)
library(patchwork)
library(shinyBS)
library(shinyjs)
library(readr)
library(stringr)
library(grid) 
library(gridExtra)
library(shinydashboard)

source("modules/mod_alphamissense.R")
source("modules/mod_consurf.R")
source("modules/mod_ptm.R")
source("modules/mod_gnomad.R")
source("modules/mod_basicinfo.R")

# Define UI for application
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Protein info Viewer"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "max-width: 350px;",
      # Basic Info
      conditionalPanel(
        condition = "input.main_tabs == 'Basic Info'",
        mod_basicinfo_ui("bi")
      ),
      
      # 1D Plot
      conditionalPanel(
        condition = "input.main_tabs == '1D Plot'",
        mod_ptm_ui("ptm1"),
        mod_consurf_1d_ui("consurf1"),
        tags$div(style = "margin-top: 15px;"), 
        uiOutput("fells_button"),
        tags$hr(),
        h4("AlphaMissense"),
        mod_alphamissense_ui("am1")
      ),
      
      # 3D Structure
      conditionalPanel(
        condition = "input.main_tabs == '3D Structure'",
        selectInput("structure_source", "Choose Structure Source",
                    choices = c("AlphaFold", "PDB", "Upload"),
                    selected = "AlphaFold"),
        conditionalPanel(
          condition = "input.structure_source == 'Upload'",
          fileInput("pdb_upload", "Upload PDB File", accept = c(".pdb"))
        ),
        uiOutput("pdb_selector"),
        selectInput("set_style", "Choose Structure Style",
                    choices = c("Cartoon", "Line", "Stick", "Sphere", "Cross"),
                    selected = "Cartoon"),
        checkboxInput("spin", "Spin Structure", value = FALSE),
        checkboxInput("surface", "Show Surface", value = FALSE),
        checkboxInput("labels", "Show Mutation Labels", value = FALSE),
        fileInput("consurf_pdb_upload", "Upload ConSurf PDB", accept = ".pdb"),
        checkboxInput("toggle_consurf_3d", "Show ConSurf on 3D", value = FALSE), 
        
        numericInput("first", "First Residue", value = 1, min = 1),
        numericInput("last", "Last Residue", value = NULL, min = 1),
        actionButton("selectSpheres", "Highlight Variants with Spheres"),
        tags$hr(),
        sliderInput(
          inputId = "set_slab",
          label = "Set slab viewing depth",
          min = -150,
          max = 150,
          value = c(-150, 150),
          step = 10,
          dragRange = TRUE,
          animate = TRUE
        ),
        helpText("Adjust to clip front/back layers of the structure")
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(id = "main_tabs",
                  tabPanel("Basic Info", uiOutput("tab_basic")),
                  tabPanel("1D Plot", uiOutput("tab_1dplot")),
                  tabPanel("3D Structure", uiOutput("tab_3d"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  bi_res <- mod_basicinfo_server("bi")
  am_res <- mod_alphamissense_server("am1", pid = reactive({ bi_res$pid()}))
  consurf_res <- mod_consurf_server("consurf1")
  ptm_res <- mod_ptm_server("ptm1", protein_data = bi_res$protein_data)
  gnomad_res <- mod_gnomad_server("gn1", 
                                  seq_source = reactive({ bi_res$seq_source() }),
                                  selected_transcript = reactive({ bi_res$selected_transcript() }),
                                  gnomad_upload = reactive({ bi_res$gnomad_upload() })
  )
  gnomad_df <- gnomad_res$gnomad_df
  
  # Basic Info Tab
  output$tab_basic <- renderUI({
    has_seq <- !is.null(bi_res$current_sequence())
    has_gnomad <- !is.null(gnomad_df()) && nrow(gnomad_df()) > 0
    ptm_data <- ptm_res$ptm_df()
    
    if (!has_seq && !has_gnomad) return(helpText("Please upload a FASTA file or a gnomAD CSV."))
    content <- list()
    
    if(has_seq){
      content <- c(content, list(
        h4("Protein Info"),
        verbatimTextOutput("info"),
        h4("Sequence"),
        verbatimTextOutput("sequence"),
        h4("Domains"),
        tableOutput("domain_table")
      ))
    }
    
    if (has_gnomad) {
      content <- c(content, list(
        h4("GnomAD Summary"),
        uiOutput("gnomad_summary"),
        h4("GnomAD Variant Table"),
        DT::dataTableOutput("gnomad_table")
      ))
    }
    
    tagList(content)
  })
  
  # 1D Plot Tab
  output$tab_1dplot <- renderUI({
    req(bi_res$current_sequence())
    tagList(
      numericInput("vdvp_window", "Window size:", value = 3, min = 0.01, step = 0.1),
      bsTooltip("vdvp_window", 
                title = "Window size: Integer = fixed length, < 1 = fraction of the protein length", 
                placement = "right", 
                trigger = "hover"),
      plotlyOutput("variant_1dplot", height="600px"),
      br(),
      uiOutput("consurf_plot_ui"),
      br(),
      conditionalPanel(
        condition = "output.fellsAvailable == true",
        h4("FELLS Result"),
        plotOutput("fells_plot", height = "500px")
      ),
      br(),
      h4("AlphaMissense Pathogenicity Heatmap"),
      plotlyOutput("alphamissense_heatmap", height = "500px")
    )
  })
  
  output$consurf_plot_ui <- renderUI({
    consurf_res$plot_ui
  })
  
  consurf_colors <- c(
    "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
    "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
  )
  
  # 3D Structure Tab
  output$tab_3d <- renderUI({
    req(input$structure_source)
    tagList(
      textOutput("structure_info"),
      withSpinner(r3dmolOutput("structure_view", height = "500px")),
      br(),
      uiOutput("sphere_legend"),
      uiOutput("consurf_legend_ui"), 
      br(),
      downloadButton("download_pdb", "Download PDB File")
    )
  })
  
  gene_name <- reactive({
    tx_df <- bi_res$transcript_table()
    tx <- bi_res$selected_transcript()
    matched_gene <- tx_df$external_gene_name[tx_df$ensembl_transcript_id == tx]
    if (length(matched_gene) > 0 && nzchar(matched_gene)) matched_gene else "Unknown Gene"
  })
  
  missense_df <- reactive({
    df <- gnomad_df()
    req(df)
    df <- subset(df, grepl("missense_variant", df$Consequence, ignore.case = TRUE))
    df <- df[!is.na(df$AA_Position), ]
    df
  })
  
  # Basic info output
  output$info <- renderPrint({
    if (bi_res$seq_source() == "upload") {
      seq <- bi_res$current_sequence()
      if (!is.null(seq)) {
        list(
          Protein_Name = "User uploaded sequence",
          Length = nchar(seq),
          Gene = "N/A",
          Organism = "N/A"
        )
      } else {
        "No sequence available."
      }
    } else {
      data <- bi_res$protein_data()
      if (!is.null(data)){
        list(
          ID = data$primaryAccession,
          Protein_Name = data$proteinDescription$recommendedName$fullName$value,
          Gene = paste(sapply(data$genes, function(g) g$geneName$value), collapse = "; "),
          Organism = data$organism$scientificName,
          Length = data$sequence$length
        )
      } else {
        "Failed to fetch, please check the UniProt ID"
      }
    }
  })
  
  # Protein sequence output
  output$sequence <- renderText({
    if (bi_res$seq_source() == "upload") {
      seq <- bi_res$current_sequence()
      if (!is.null(seq)) {
        return(seq)
      } else {
        return("Failed to load uploaded sequence.")
      }
    }
    
    data <- bi_res$protein_data()
    if (!is.null(data)){
      data$sequence$value
    } else {
      "Failed to fetch, please check the UniProt ID"
    }
  })
  
  # Domain info table output
  output$domain_table <- renderTable({
    data <- bi_res$protein_data()
    if (!is.null(data) && !is.null(data$features)) {
      domain_features <- data$features[sapply(data$features, function(x) x$type == "Domain")]
      
      if (length(domain_features) > 0) {
        data.frame(
          Description = sapply(domain_features, function(x) x$description),
          Start = as.integer(sapply(domain_features, function(x) x$location$start$value)),
          End = as.integer(sapply(domain_features, function(x) x$location$end$value))
        )
      } else {
        data.frame(Message = "No DOMAIN features found.")
      }
    } else {
      data.frame(Message = "No feature data available.")
    }
  })
  
  output$gnomad_summary <- gnomad_res$summary_ui
  output$gnomad_table <- gnomad_res$table_ui
  
  # 1D plot Output
  output$variant_1dplot <- renderPlotly({
    req(bi_res$current_sequence())
    seq_val <- bi_res$current_sequence()
    data <- if (bi_res$seq_source() == "api") bi_res$protein_data() else NULL
    
    protein_len_df <- data.frame(start = 0, end = nchar(seq_val),
                                 ymin = 0.45, ymax = 0.55,
                                 label = paste0("Protein length: ", nchar(seq_val)))
    ## get the domain info
    domain_df <- NULL
    if (bi_res$seq_source() == "api" && !is.null(data$features)) {
      domains <- data$features[sapply(data$features, function(x) x$type == "Domain")]
      if (length(domains) == 0) return(NULL)
      
      domain_df <- data.frame(
        start = as.integer(sapply(domains, function(x) x$location$start$value)),
        end   = as.integer(sapply(domains, function(x) x$location$end$value)),
        description = sapply(domains, function(x) x$description)
      )
    }
    ## get missense variant info
    missense_df <- missense_df()
    
    ## get PTM info
    ptm <- ptm_res$ptm_df()
    
    p1 <- ggplot() +
      geom_rect(data = protein_len_df, 
                aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, text = label), 
                fill = "grey90")
    
    if (!is.null(domain_df)) {
      p1 <- p1 +
        geom_rect(data = domain_df,
                  aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55,
                      text = paste0("Domain: ", description, "\n", start, " - ", end)),
                  fill = "#fca6a6")
    }
    
    if (nrow(missense_df) > 0) {
      p1 <- p1 +
        geom_linerange(data = missense_df,
                       aes(x = AA_Position, ymin = 0.6, ymax = 0.7),
                       color = "slateblue", size = 0.3, alpha = 0.5)
    }
    
    if (!is.null(ptm) && nrow(ptm) > 0  && isTRUE(ptm_res$show())) {
      ptm_plot <- if (!is.null(ptm_res$selected_types())) {
        ptm |> dplyr::filter(TypeCategory %in% ptm_res$selected_types())
      } else {
        ptm |> dplyr::filter(FALSE) 
      }
      
      ptm_plot <- ptm |> dplyr::filter(TypeCategory %in% ptm_res$selected_types())
      
      if (nrow(ptm_plot) > 0) {
        ptm_color_map <- c(
          "Phosphorylation" = "#1f77b4",
          "Acetylation" = "#ff7f0e",
          "Succinylation" = "#2ca02c",
          "Methylation" = "#d62728",
          "Other" = "#9467bd"
        )
        
        p1 <- p1 +
          geom_point(data = ptm_plot,
                     aes(x = Position, y = 0.75, fill = TypeCategory, text = tooltip),
                     shape = 21, size = 2, color = "black", stroke = 0.3, alpha = 0.8)+
          scale_fill_manual(values = ptm_color_map, name = "PTM Type", drop = FALSE)
      }
    }
    
    if (consurf_res$enabled() && !is.null(consurf_res$consurf_df())) {
      consurf_data <- consurf_res$consurf_df()
      consurf_data <- consurf_data %>% filter(!is.na(Position), !is.na(Grade))
      consurf_data$Color <- consurf_colors[as.character(consurf_data$Grade)]
      
      p1 <- p1 + 
        geom_linerange(data = consurf_data,
                       mapping = aes(x = Position, ymin = 0.3, ymax = 0.4, 
                                     text = paste0("ConSurf\n", "Position: ", Position, "\n","Score: ", round(Score, 3), "\n","Grade: ", Grade)),
                       inherit.aes = FALSE,
                       color = consurf_data$Color,
                       size = 1.0, alpha = 0.9)
    }
    
    if (am_res$enabled() && !is.null(am_res$alphamissense_df())) {
      alpha_bar_df <- am_res$alphamissense_df() %>%
        dplyr::group_by(position) %>%
        dplyr::summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop")
      
      p1 <- p1 +
        geom_linerange(
          data = alpha_bar_df,
          aes(
            x = position,
            ymin = 0.15,
            ymax = avg_score * 0.1 + 0.15,
            text = paste0("AlphaMissense\nPosition: ", position, "\nAvg Score: ", round(avg_score, 3))
          ),
          color = "darkgreen",
          size = 1.0,
          alpha = 0.8)
    }
    
    p1 <- p1 +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
      ) +
      labs(x = NULL)
    
    mut_index <- missense_df$AA_Position |> unique()
    prot_len <- nchar(seq_val)
    # window <- 3
    user_input <- input$vdvp_window
    if (user_input < 1) {
      window <- max(1, floor(prot_len * user_input))  # fraction
    } else {
      window <- as.integer(user_input)  # fixed
    }
    
    vp <- length(mut_index) / prot_len
    
    index <- 1:prot_len
    vdvp_vals <- sapply(index, function(x) {
      vd <- sum(mut_index %in% x:(x + window))
      vd / window / vp
    })
    
    vdvp_df <- data.frame(x = index + window / 2, y = vdvp_vals)
    vdvp_df <- as.data.frame(spline(vdvp_df$x, vdvp_df$y))
    
    p2 <- ggplot(vdvp_df, aes(x = x, y = y)) +
      geom_line(color = "steelblue", size = 0.6) +
      theme_minimal() +
      labs(x = "Residue", y = "Vd/Vp") +
      theme(
        plot.title = element_blank()
      )
    
    p1_plotly <- ggplotly(p1, tooltip = "text") %>% layout(margin = list(b = 0), showlegend = TRUE, legend = list(orientation = "h", x = 0, xanchor = "left", y = 1.1))
    p2_plotly <- ggplotly(p2) %>% layout(margin = list(t = 0))
    
    title_text <- if (bi_res$seq_source() == "api" && !is.null(data)) {
      paste0("Protein: ", data$proteinDescription$recommendedName$fullName$value,
             "; Gene: ", gene_name())
    } else {
      "User uploaded sequence"
    }
    
    has_extra_layer <- consurf_res$enabled() || am_res$enabled()
    subplot_heights <- if (has_extra_layer) c(0.7, 0.3) else c(0.5, 0.5)
    
    subplot(p1_plotly, p2_plotly, nrows = 2, shareX = TRUE, titleY = TRUE,heights = subplot_heights) %>%
      layout(title = list(text = title_text, x = 0, xanchor = "left"),
             margin = list(t = 60))
  })
  
  # Alphamissense Heatmap Plot
  output$alphamissense_heatmap <- renderPlotly({
    df <- am_res$alphamissense_df()
    req(nrow(df) > 0)
    
    # filter prediction
    if (!is.null(am_res$prediction_filter()) && length(am_res$prediction_filter()) > 0) {
      if ("am_class" %in% colnames(df)) {
        df <- df[df$am_class %in% am_res$prediction_filter(), ]
      }
    }
    req(nrow(df) > 0)
    
    # Dynamic title
    file_name <- am_res$upload_info()$name
    dynamic_title <- if (!is.null(bi_res$pid()) && nzchar(bi_res$pid())) {
      paste0("AlphaMissense Heatmap - ", bi_res$pid())
    } else {
      paste0("AlphaMissense Heatmap - ", file_name)
    }
    
    # Ref aa
    ref_df <- df %>%
      select(position, aa_from) %>%
      distinct() %>%
      mutate(aa_from = factor(aa_from, levels = levels(df$aa_to)))
    
    p <- ggplot()+
      geom_tile(
        data = df, aes(x = position, y = aa_to, fill = score,
                       text = paste0(
                         "Position: ", position, "<br>",
                         "From: ", aa_from, " → ", aa_to, "<br>",
                         "Score: ", round(score, 3), "<br>",
                         "Prediction: ", am_class
                       )), color = NA) +
      geom_tile(
        data = ref_df %>%
          mutate(text = paste0("Variant: ", aa_from, position, " (reference)")), 
        aes(x = position, y = aa_from, text = text),
        inherit.aes = FALSE,
        fill = "black", width = 1, height = 1, alpha = 0.9) +
      scale_fill_gradientn(
        colors = c("blue", "white", "red"),
        values = scales::rescale(c(0, 0.5, 1)),
        limits = c(0, 1),
        name = "Score"
      ) +
      labs(
        title = dynamic_title,
        x = "Residue Position",
        y = "Alternative Amino Acid"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank()
      )
    ggplotly(p, tooltip = "text")
  })
  
  
  # Read Uploaded PDB File
  uploaded_pdb <- reactive({
    req(input$pdb_upload)
    tryCatch({
      pdb <- bio3d::read.pdb(input$pdb_upload$datapath)
      m_bio3d(pdb)
    }, error = function(e) {
      showNotification("Failed to parse uploaded PDB file", type = "error")
      return(NULL)
    })
  })
  
  # control residue length dynamically
  observe({
    seq <- bi_res$current_sequence()
    if (!is.null(seq)) {
      len <- nchar(seq)
      updateNumericInput(session, "first", max = len, value = 1)
      updateNumericInput(session, "last", max = len, value = len)
    }
  })
  
  # Alphafold API function
  fetch_alphafold_prediction <- function(uniprot_id) {
    url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot_id)
    res <- GET(url)
    if (status_code(res) == 200) {
      json <- content(res, as = "text", encoding = "UTF-8")
      parsed <- fromJSON(json, simplifyVector = FALSE)
      return(parsed)
    } else {
      warning(paste("AlphaFold API failed with status:", status_code(res)))
      NULL
    }
  }
  
  # Function to get PDB id from Uniprot id
  get_pdb_ids_from_uniprot <- function(uniprot_id) {
    url <- paste0("https://www.ebi.ac.uk/proteins/api/proteins/", uniprot_id)
    res <- httr::GET(url, httr::add_headers(Accept = "application/json"))
    if (httr::status_code(res) == 200) {
      data <- httr::content(res, as = "parsed", type = "application/json")
      pdbs <- data$dbReferences[sapply(data$dbReferences, function(x) x$type == "PDB")]
      if (length(pdbs) > 0) {
        return(sapply(pdbs, function(x) x$id))
      }
    }
    return(character(0))
  }
  
  # Get all the PDB id in the list
  pdb_id_list <- reactive({
    req(bi_res$pid())
    get_pdb_ids_from_uniprot(bi_res$pid())
  })
  
  output$pdb_selector <- renderUI({
    req(input$structure_source == "PDB")
    ids <- pdb_id_list()
    if (length(ids) > 0) {
      selectInput("selected_pdb_id", "Select PDB Structure", choices = ids, selected = ids[1])
    } else {
      helpText("No PDB structures available.")
    }
  })
  
  structure_source_real <- reactiveVal(NULL)
  
  observeEvent(bi_res$selected_transcript(), {
    req(bi_res$canonical_tx(), bi_res$selected_transcript())
    
    if (bi_res$canonical_tx() != bi_res$selected_transcript()) {
      showNotification(
        paste0("Current 3D structure is based on canonical transcript: ", bi_res$canonical_tx(),
               ", but you selected: ", bi_res$selected_transcript(), 
               ". Mapping may be inaccurate."),
        type = "warning",
        duration = 8
      )
    }
  })
  
  # 3D structure with r3dmol
  output$structure_view <- renderR3dmol({
    req(input$structure_source)
    
    if (input$structure_source == "Upload") {
      shiny::validate(
        need(!is.null(input$pdb_upload), "Please upload a PDB file.")
      )
    } else {
      shiny::validate(
        need(input$structure_source != "Upload" && (!isTruthy(input$fetch)), "Please click Fetch after entering UniProt ID.")
      )
    }
    
    uniprot_id <- bi_res$pid()
    pdb_data <- NULL
    
    if (input$structure_source == "AlphaFold") {
      structure_source_real("AlphaFold")
      pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",uniprot_id,"-F1-model_v4.pdb")))
    } else if (input$structure_source == "PDB") {
      req(input$selected_pdb_id)
      pdb_id <- input$selected_pdb_id
      
      tryCatch({
        pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://files.rcsb.org/download/", pdb_id, ".pdb")))
        structure_source_real("PDB")
      }, error = function(e){
        showNotification("No PDB structure available for this protein, fallback to AlphaFold.", type = "error")
        pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",bi_res$pid(),"-F1-model_v4.pdb")))
        structure_source_real("AlphaFold")
      })
    } else if (input$structure_source == "Upload") {
      pdb_data <- uploaded_pdb()
      structure_source_real("Upload")
    }
    
    r3dmol(
      viewer_spec = m_viewer_spec(backgroundColor = "white", 
                                  cartoonQuality = 25,
                                  lowerZoomLimit = 5,
                                  upperZoomLimit = 1000)
    ) %>%
      m_add_model(data = pdb_data, format = "pdb") %>%
      m_set_style(style = m_style_cartoon(color = "spectrum")) %>%
      m_zoom_to()
  })
  
  # Text about 3D source
  output$structure_info <- renderText({
    req(bi_res$pid())
    src <- structure_source_real()
    
    if (src == "AlphaFold") {
      paste("Structure from AlphaFold for Uniprot ID:", bi_res$pid())
    } else if(src == "PDB"){
      paste("Structure from PDB:", input$selected_pdb_id)
    } else if(src == "Upload"){
      "Structure from uploaded PDB file."
    } else{
      "No PDB structure found for this UniProt ID."
    }
  })
  
  # Download PDB button
  output$download_pdb <- downloadHandler(
    filename = function() {
      src <- structure_source_real()
      uniprot_id <- bi_res$pid()
      
      if(src == "AlphaFold"){
        af_data <- fetch_alphafold_prediction(uniprot_id)
        if (!is.null(af_data) && !is.null(af_data[[1]]$entryId)) {
          paste0(af_data[[1]]$entryId, ".pdb")
        } else {
          paste0("AF-", uniprot_id, "-model.pdb")
        }
      } else if(src == "PDB"){
        pdb_id <- input$selected_pdb_id
        if (!is.null(pdb_id)) {
          return(paste0(pdb_id, ".pdb"))
        } else {
          return("pdb_structure.pdb")
        }
      } else if (src == "Upload") {
        if (!is.null(input$pdb_upload)) {
          return(input$pdb_upload$name)
        } else {
          return("uploaded_structure.pdb")
        }
      }
    },
    
    content = function(file) {
      src <- structure_source_real()
      uniprot_id <- bi_res$pid()
      if (src == "AlphaFold") {
        af_data <- fetch_alphafold_prediction(uniprot_id)
        if (!is.null(af_data) && !is.null(af_data[[1]]$pdbUrl)) {
          pdb_url <- af_data[[1]]$pdbUrl
          download.file(pdb_url, destfile = file, mode = "wb")
        } else {
          showNotification("AlphaFold model not available for this protein.", type = "error")
        }
      } else if(src == "PDB"){
        pdb_id <- input$selected_pdb_id
        if (!is.null(pdb_id)) {
          pdb_url <- paste0("https://files.rcsb.org/download/", pdb_id, ".pdb")
          download.file(pdb_url, destfile = file, mode = "wb")
        } else {
          showNotification("No PDB ID selected for download.", type = "error")
        }
      } else if (src == "Upload") {
        req(input$pdb_upload)
        file.copy(input$pdb_upload$datapath, file)
      }
    }
  )
  
  observeEvent(input$set_style, {
    style <- switch(input$set_style,
                    "Line" = list(line = list()),
                    "Cartoon" = list(cartoon = list(color = "spectrum")),
                    "Stick" = list(stick = list()),
                    "Cross" = list(cross = list()),
                    "Sphere" = list(sphere = list()))
    m_set_style(id = "structure_view", style = style)
  })
  
  observeEvent(input$spin, {
    m_spin(id = "structure_view", speed = ifelse(input$spin, 0.3, 0))
  })
  
  observeEvent(input$surface, {
    if (input$surface) {
      m_add_surface(id = "structure_view", style = m_style_surface(opacity = 0.4))
    } else {
      m_remove_all_surfaces(id = "structure_view")
    }
  })
  
  missense_df_3d <- reactive({
    df <- missense_df()
    req(df, input$first, input$last)
    
    df <- df[df$AA_Position >= input$first & df$AA_Position <= input$last, ]
    df
  })
  
  extract_protein_variant <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "") return(NA)
    matches <- stringr::str_match(hgvsp, "^p\\.([A-Za-z]{3})(\\d+)([A-Za-z]{3})$")
    if (any(is.na(matches))) return(NA)
    
    aa3to1 <- c(
      Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
      Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
      Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
      Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
      Ter="*", Sec="U", Pyl="O", Asx="B", Glx="Z", Xaa="X"
    )
    paste0(aa3to1[[matches[2]]], matches[3], aa3to1[[matches[4]]])
  }
  
  alphamissense_color_map <- colorRampPalette(c("blue", "white", "red"))
  score_to_color <- function(scores) {
    ramp <- alphamissense_color_map(100)
    index <- pmin(100, pmax(1, round(scores * 100)))
    ramp[index]
  }
  
  # Labels about missense position
  observeEvent(input$labels, {
    req(input$labels %in% c(TRUE, FALSE))
    resis <- missense_df_3d()$AA_Position |> unique()
    req(length(resis) > 0)
    
    if (input$labels) {
      m_add_res_labels(
        id = "structure_view",
        sel = m_sel(resi = resis),
        style = m_style_label(
          backgroundColor = "#FF6F61",
          inFront = TRUE,
          fontSize = 12,
          showBackground = TRUE
        )
      )
    } else {
      m_remove_all_labels(id = "structure_view")
    }
  })
  
  # Highlight Variants with Spheres
  cols_3d <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd", "#8c564b")
  radii_3d <- c(1.15, 1.5, 1.85, 2.15, 2.5, 2.85)
  render_legend <- reactiveVal(FALSE)
  
  observeEvent(input$selectSpheres, {
    render_legend(TRUE)
    df <- missense_df_3d()
    req(nrow(df) > 0)
    # AF groups
    df <- df %>%
      filter(!is.na(AF)) %>%
      mutate(
        LogAF = ifelse(AF > 0, log10(AF * 1e6), -6),
        AF_Group = case_when(
          LogAF <= 1 ~ 1,
          LogAF <= 2 ~ 2,
          LogAF <= 3 ~ 3,
          LogAF <= 4 ~ 4,
          LogAF <= 5 ~ 5,
          TRUE       ~ 6)
      )
    
    # with alphamissense
    if (am_res$am_source() == "upload" && !is.null(am_res$upload_info())) {
      alpha_df <- am_res$alphamissense_df()
      df$protein_variant <- sapply(df$HGVSp, extract_protein_variant)
      
      df <- df %>%
        left_join(alpha_df %>% select(protein_variant, score), by = "protein_variant")
      
      df$Color <- ifelse(is.na(df$score), "#AAAAAA", score_to_color(df$score))
      
      # keep the variants with highest score
      df_grouped <- df %>%
        group_by(AA_Position) %>%
        summarise(
          radius = radii_3d[max(AF_Group, na.rm = TRUE)],
          color = Color[which.max(score)],
          .groups = "drop")
    } else if(am_res$am_source() == "api"){
      alpha_df <- am_res$alphamissense_df()
      if (!is.null(alpha_df) && nrow(alpha_df) > 0) {
        df$protein_variant <- sapply(df$HGVSp, extract_protein_variant)
        df <- df %>%
          left_join(alpha_df %>% select(protein_variant, score), by = "protein_variant")
        df$Color <- ifelse(is.na(df$score), "#AAAAAA", score_to_color(df$score))
        df_grouped <- df %>%
          group_by(AA_Position) %>%
          summarise(
            radius = radii_3d[max(AF_Group, na.rm = TRUE)],
            color = Color[which.max(score)],
            .groups = "drop")
      } else {
        df_grouped <- df %>%
          group_by(AA_Position) %>%
          summarise(
            af_group = max(AF_Group, na.rm = TRUE),
            radius = radii_3d[af_group],
            color = cols_3d[af_group],
            .groups = "drop")
      }
    } else {
      # no AlphaMissense (upload empty or API fail)
      df_grouped <- df %>%
        group_by(AA_Position) %>%
        summarise(
          af_group = max(AF_Group, na.rm = TRUE),
          radius = radii_3d[af_group],
          color = cols_3d[af_group],
          .groups = "drop")
    }
    
    # Spheres
    for (i in seq_len(nrow(df_grouped))) {
      m_add_style(
        id = "structure_view",
        sel = m_sel(resi = df_grouped$AA_Position[i], atom = "CA"),
        style = m_style_sphere(
          color = df_grouped$color[i],
          radius = df_grouped$radius[i],
          colorScheme = "none"))}
  })
  
  # Spheres legend
  output$am_score_bar <- renderImage({
    outfile <- tempfile(fileext = ".png")
    scores <- seq(0, 1, length.out = 100)
    df <- data.frame(x = scores, y = 1, score = scores)
    
    p <- ggplot(df, aes(x = x, y = y, fill = score)) +
      geom_tile() +
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "white", "red"))(100), limits = c(0, 1)) +
      scale_x_continuous(
        breaks = seq(0, 1, 0.25),
        labels = sprintf("%.2f", seq(0, 1, 0.25)),
        expand = c(0, 0)
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 10, 5, 10)
      )
    
    label_grob <- gridExtra::grid.arrange(
      textGrob("Benign", gp = gpar(fontsize = 10), just = "left"),
      textGrob("Pathogenic", gp = gpar(fontsize = 10), just = "right"),
      ncol = 2,
      widths = unit(c(0.5, 0.5), "npc")
    )
    
    png(outfile, width = 500, height = 100)
    gridExtra::grid.arrange(label_grob, p, ncol = 1, heights = c(0.3, 1))
    dev.off()
    
    list(src = outfile, contentType = "image/png", width = 500, height = 100)
  }, deleteFile = TRUE)
  
  output$sphere_legend <- renderUI({
    if (!render_legend()) return(NULL)
    show_am <- FALSE
    if (am_res$am_source() == "upload" && !is.null(am_res$upload_info())) {
      show_am <- TRUE
    } else if (am_res$am_source() == "api") {
      df <- am_res$alphamissense_df()
      show_am <- !is.null(df) && nrow(df) > 0
    }
    if (show_am) {
      tagList(
        tags$h5("Legend: AlphaMissense + AF-value Mode"),
        tags$p("Color → AlphaMissense pathogenicity score"),
        imageOutput("am_score_bar", height = "100px"),
        tags$p("Blue → benign   |   Red → pathogenic"),
        tags$p("Sphere size → AF Group (mutation frequency)"),
        tags$ul(
          lapply(1:6, function(i) {
            af_labels <- c("Extremely rare","Very rare","Rare","Low freq.","Common","Very common")
            tags$li(paste("Group", i, "→ radius =", radii_3d[i], "→", af_labels[i]))})))
    } else {
      tagList(
        tags$h5("Legend: AF-value-only Mode"),
        tags$p("Color → AF Group"),
        tags$table(
          lapply(1:6, function(i) {
            af_labels <- c("Extremely rare","Very rare","Rare","Low freq.","Common","Very common")
            tags$tr(
              tags$td(style = paste0("background-color:", cols_3d[i], ";width:20px;height:20px; border:1px solid #000;"), ""),
              tags$td(paste("Group", i, "(radius =", radii_3d[i], ")", "→", af_labels[i]))
            )})))}
  })
  
  output$consurf_legend_ui <- renderUI({
    req(input$toggle_consurf_3d)
    
    if (isTRUE(input$toggle_consurf_3d)) {
      consurf_colors <- c(
        "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
        "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
      )
      
      tags$div(
        style = "margin-top: 10px;",
        strong("ConSurf Grade Legend:"),
        tags$div(
          style = "display: flex; gap: 8px; flex-wrap: wrap; margin-top: 5px;",
          lapply(names(consurf_colors), function(g) {
            tags$div(
              style = paste0(
                "background-color:", consurf_colors[[g]], ";",
                "width: 28px; height: 22px; text-align: center;",
                "border: 1px solid #999; font-size: 12px; line-height: 22px;"
              ),
              g
            )
          })
        )
      )
    } else {
      NULL
    }
  })
  
  uploaded_consurf_pdb <- reactive({
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
  
  observeEvent(input$toggle_consurf_3d, {
    req(input$toggle_consurf_3d %in% c(TRUE, FALSE))
    req(consurf_res$consurf_df(), uploaded_consurf_pdb()) 
    
    consurf_data <- consurf_res$consurf_df() %>%
      dplyr::filter(!is.na(Position), !is.na(Grade))
    
    pdb_resnos <- unique(uploaded_consurf_pdb()$bio3d$atom$resno)
    consurf_data <- consurf_data %>% dplyr::filter(Position %in% pdb_resnos)
    
    consurf_colors <- c(
      "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
      "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
    )
    
    if (isTRUE(input$toggle_consurf_3d)) {
      for (i in seq_len(nrow(consurf_data))) {
        resi <- consurf_data$Position[i]
        grade <- as.character(consurf_data$Grade[i])
        color <- consurf_colors[[grade]]
        if (!is.null(color)) {
          m_add_style(
            id = "structure_view",
            sel = m_sel(resi = resi),
            style = m_style_cartoon(color = color)
          )
        }
      }
      showNotification("ConSurf coloring applied to 3D structure.", type = "message")
    } else {
      m_set_style(
        id = "structure_view",
        style = m_style_cartoon(color = "spectrum")
      )
      showNotification("ConSurf coloring removed from 3D structure.", type = "default")
    }
  }) 
  
  observeEvent(input$set_slab, {
    m_set_slab(
      id = "structure_view",  # same as the r3dmol Output id
      near = input$set_slab[1],
      far = input$set_slab[2]
    )
  })
  
  fells_job <- reactiveVal()
  fells_result <- reactiveVal()
  fells_class <- reactiveVal("btn-primary")
  
  output$fells_button <- renderUI({
    req(bi_res$current_sequence())
    actionButton("fells", label = "Submit to FELLS", class = fells_class())
  })
  
  output$fellsAvailable <- reactive({
    fells_class() == "btn-success"
  })
  outputOptions(output, "fellsAvailable", suspendWhenHidden = FALSE)
  
  observeEvent(input$fells, {
    data <- bi_res$protein_data()
    req(data)
    
    seq <- bi_res$current_sequence()
    pseq <- paste0(">", bi_res$pid(), "\n", seq)
    
    req <- request("http://protein.bio.unipd.it/fellsws/submit") |> 
      req_headers("Content-Type" = "multipart/form-data") |>
      req_body_multipart(sequence = pseq)
    
    resp <- req_perform(req)
    rb <- resp |> resp_body_json()
    fells_job(rb$jobid)
    
    showModal(modalDialog(title = "FELLS job submitted"))
    fells_class("btn-warning")
    updateActionButton(inputId = "fells", label = "FELLS submitted")
    
    job <- fells_job()
    status <- "submitted"
    
    while (status != "done") {
      Sys.sleep(10)
      check <- request(paste0("http://protein.bio.unipd.it/fellsws/status/", job))
      resp <- req_perform(check)
      rb <- resp |> resp_body_json()
      status <- rb$status
    }
    
    if (status == "done") {
      result_id <- rb$names[[1]][[2]]
      req2 <- request(paste0("http://protein.bio.unipd.it/fellsws/result/", result_id))
      res2 <- req_perform(req2)
      rb2 <- res2 |> resp_body_json()
      fells_result(rb2)
      
      fells_class("btn-success")
      updateActionButton(inputId = "fells", label = "FELLS complete")
      
      showModal(modalDialog(title = "FELLS job complete"))
    } else {
      fells_class("btn-danger")
      updateActionButton(inputId = "fells", label = "FELLS failed")
      showModal(modalDialog(title = "FELLS job failed"))
    }
  })
  
  fells_hsc <- reactive({
    fr <- fells_result()
    req(fr)
    
    phelix <- fr$p_h |> unlist()
    pstrand <- fr$p_e |> unlist()
    pcoil <- fr$p_c |> unlist()
    
    df <- data.frame(
      index = seq_along(phelix),
      Helix = as.numeric(phelix),
      Strand = as.numeric(pstrand),
      Coil = as.numeric(pcoil)
    ) |> 
      mutate(Coil = Coil * -1) |>
      tidyr::pivot_longer(-index) |>
      mutate(alpha = ifelse(name == "Coil", -value, value),
             alphag = cut(alpha, breaks = 5))
    
    df
  })
  
  fells_hd <- reactive({
    fr <- fells_result()
    req(fr)
    
    hca <- fr$hca |> unlist()
    dis <- fr$p_dis |> unlist()
    
    df <- data.frame(
      index = seq_along(hca),
      HCA = as.numeric(hca),
      Disorder = as.numeric(dis)
    ) |> 
      mutate(Disorder = Disorder * -1) |>
      tidyr::pivot_longer(-index) |>
      mutate(alpha = ifelse(name == "Disorder", -value, value),
             alphag = cut(alpha, breaks = 5))
    
    df
  })
  
  output$fells_plot <- renderPlot({
    req(fells_class() == "btn-success")
    
    hsc_cols <- c(Helix = "#91288c", Strand = "#ffa500", Coil = "gray50")
    dh_cols <- c(Disorder = "red", HCA = "black")
    
    p1 <- ggplot(fells_hsc(), aes(index, value, fill = name, alpha = alpha)) +
      geom_col(position = "identity") +
      scale_fill_manual(values = hsc_cols) +
      guides(alpha = "none") +
      theme_minimal() +
      ggtitle("Secondary Structure Prediction")
    
    p2 <- ggplot(fells_hd(), aes(index, value, fill = name, alpha = alpha)) +
      geom_col(position = "identity") +
      scale_fill_manual(values = dh_cols) +
      guides(alpha = "none") +
      theme_minimal() +
      ggtitle("Disorder & HCA")
    
    patchwork::wrap_plots(p1, p2, ncol = 1)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)