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

# Define UI for application
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Protein info Viewer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          # Basic Info
          conditionalPanel(
            condition = "input.main_tabs == 'Basic Info'",
            radioButtons("seq_source", "Sequence Source",
                         choices = c("API" = "api", "Upload File" = "upload"),
                         selected = "api",
                         inline = TRUE),
            conditionalPanel(
              condition = "input.seq_source == 'upload'",
              fileInput("fasta_upload", "Upload FASTA File", 
                        accept = c(".fasta", ".fa"),
                        placeholder = "Upload a FASTA file")
            ),
            conditionalPanel(
              condition = "input.seq_source == 'upload'",
              fileInput("gnomad_upload", "Upload gnomAD-like CSV",
                        accept = c(".csv"),
                        placeholder = "Optional variant CSV with HGVSp column"),
              helpText("Optional. Must contain HGVSp column like 'p.Arg12Cys'.")
            ),
            textInput("pid", "Enter UniProt ID"),
            actionButton("fetch", "Fetch info"),
            tags$div(style = "margin-top: 10px;"), 
            conditionalPanel(
              condition = "input.seq_source == 'api'",
              uiOutput("transcript_selector")
            )
          ),
          
          # 1D Plot
          conditionalPanel(
            condition = "input.main_tabs == '1D Plot'",
            uiOutput("ptm_control_ui"),
            uiOutput("ptm_filter_ui"),
            fileInput("consurf_txt", "Upload ConSurf TXT File", accept = c(".txt")),
            actionButton("add_consurf", "Add ConSurf Layer", icon = icon("plus")),
            tags$div(style = "margin-top: 15px;"), 
            uiOutput("fells_button")
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
            actionButton("selectSpheres", "Highlight Variants with Spheres")
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

# Function to query GnomAD GraphQL API
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
                    exome {
                        af
                        ac
                        an
                    }
                        
            }
        }
    }'
    body <- list(
        query = query,
        variables = list(
            transcriptID = transcript_id
        )
    )
    
    res <- POST(url, body = body, encode = "json", content_type_json())
    if (status_code(res) == 200) {
        content(res, as = "parsed", simplifyVector = TRUE)
    } else {
        warning(paste("gnomAD API request failed with status:", status_code(res)))
        NULL
    }
}

# Get transcript ID by Uniprot ID
get_transcripts_by_uniprot <- function(uniprot_id) {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    
    result <- getBM(
        attributes = c("uniprotswissprot", "ensembl_transcript_id", "external_gene_name", "transcript_biotype"),
        filters = "uniprotswissprot",
        values = uniprot_id,
        mart = ensembl
    )
    
    # option: only keep transcript about protein coding
    protein_tx <- subset(result, transcript_biotype == "protein_coding")
    protein_tx <- unique(protein_tx)
    return(protein_tx)
}

# Define server logic
server <- function(input, output, session) {
    # Basic Info Tab
    output$tab_basic <- renderUI({
        has_seq <- !is.null(current_sequence())
        has_gnomad <- !is.null(gnomad_df()) && nrow(gnomad_df()) > 0
        ptm_data <- ptm_df()
        
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
        req(current_sequence())
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
            )
        )
    })
    
    # 3D Structure Tab
    output$tab_3d <- renderUI({
        req(input$structure_source)
        tagList(
            withSpinner(r3dmolOutput("structure_view", height = "500px")),
            br(),
            textOutput("structure_info"),
            downloadButton("download_pdb", "Download PDB File")
        )
    })

    observe({
      req(input$seq_source)
        if (input$seq_source == "upload") {
            shinyjs::disable("fetch")
        } else {
            shinyjs::enable("fetch")
        }
    })
    
    # Get Uniprot Json Data
    protein_data <- eventReactive(input$fetch,{
        url <- paste0("https://rest.uniprot.org/uniprotkb/", input$pid, ".json")
        res <- GET(url)
        if (status_code(res) == 200) {
            content(res, as = "parsed", type = "application/json")
        } else {
            NULL
        }
    })
    
    # Read FASTA File
    uploaded_fasta_seq <- reactive({
        req(input$fasta_upload)
        tryCatch({
            fasta_lines <- readLines(input$fasta_upload$datapath)
            seq_lines <- fasta_lines[!grepl("^>", fasta_lines)]
            paste(seq_lines, collapse = "")
        }, error = function(e) {
            showNotification("Failed to read FASTA file. Please check the format.", type = "error")
            return(NULL)
        })
    })
    
    current_sequence <- reactive({
      req(input$seq_source)
        if (input$seq_source == "upload") {
            if (!is.null(input$fasta_upload)) {
                return(uploaded_fasta_seq())
            } else {
                return(NULL)
            }
        } else {
            data <- protein_data()
            if (!is.null(data)) data$sequence$value else NULL
        }
    })
    
    # Get the transcripts ID
    output$transcript_selector <- renderUI({
        tx_df <- transcript_table()
        if (nrow(tx_df) == 0) return(helpText("No transcript IDs found."))
        
        choices <- setNames(tx_df$ensembl_transcript_id, 
                            paste0(tx_df$ensembl_transcript_id, " (", tx_df$external_gene_name, ")"))
        
        selectInput("selected_transcript", "Select Transcript", choices = choices)
    })
    
    # Get Gnomad Data
    uploaded_gnomad_df <- reactive({
        req(input$gnomad_upload)
        
        tryCatch({
            df_raw <- read.csv(input$gnomad_upload$datapath, stringsAsFactors = FALSE, check.names = FALSE)
            
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
            
            # print(colnames(df_raw))
            # print(colnames(df))
            
            # Extract amino acid position
            df$AA_Position <- sapply(df$HGVSp, function(hgvsp) {
                if (!is.na(hgvsp) && grepl("^p\\.", hgvsp)) {
                    matches <- regmatches(hgvsp, regexec("\\d+", hgvsp))[[1]]
                    if (length(matches) > 0) as.integer(matches[1]) else NA
                } else {
                    NA
                }
            })
            
            df
            # print(head(df))
        }, error = function(e) {
            showNotification("Failed to read gnomAD CSV file.", type = "error")
            return(NULL)
        })
    })
    
    gnomad_df <- reactive({
      req(input$seq_source)
        if (input$seq_source == "upload") {
            uploaded_gnomad_df()
        } else {
            req(input$selected_transcript)
        
            gnomad_data <- fetch_gnomad_by_transcript(input$selected_transcript)
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
    
    transcript_table <- reactive({
        req(input$pid)
        get_transcripts_by_uniprot(input$pid)
    })
    
    gene_name <- reactive({
        tx_df <- transcript_table()
        tx <- input$selected_transcript
        matched_gene <- tx_df$external_gene_name[tx_df$ensembl_transcript_id == tx]
        if (length(matched_gene) > 0 && nzchar(matched_gene)) matched_gene else "Unknown Gene"
    })
    
    # Read txt File from Consurf
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
      }, error = function(e){
        showNotification("Failed to read ConSurf TXT file", type = "error")
        NULL
      })
    })
    
    consurf_enabled <- reactiveVal(FALSE)
    observe({
      if (!is.null(input$consurf_txt)) {
        shinyjs::enable("add_consurf")
      } else {
        shinyjs::disable("add_consurf")
      }
    })
    observeEvent(input$add_consurf, {
      req(consurf_txt_df())  
      consurf_enabled(TRUE)  
    })
    
    # Basic info output
    output$info <- renderPrint({
        if (input$seq_source == "upload") {
            seq <- uploaded_fasta_seq()
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
            data <- protein_data()
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
        if (input$seq_source == "upload") {
            seq <- uploaded_fasta_seq()
            if (!is.null(seq)) {
                return(seq)
            } else {
                return("Failed to load uploaded sequence.")
            }
        }
        
        data <- protein_data()
        if (!is.null(data)){
            data$sequence$value
        } else {
            "Failed to fetch, please check the UniProt ID"
        }
    })
    
    # Domain info table output
    output$domain_table <- renderTable({
        data <- protein_data()
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

    
    # GnomAD Summary info Output
    output$gnomad_summary <- renderUI({
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
    
    # Gnomad data table display by DT
    output$gnomad_table <- DT::renderDT({
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
    
    # Function about PTM (Post-Translational Modifications)
    get_ptm_info <- function(protein_json) {
        if (is.null(protein_json$features)) return(NULL)
        
        ptm_features <- protein_json$features[sapply(protein_json$features, function(x) x$type == "Modified residue")]
        
        if (length(ptm_features) == 0) return(NULL)
        
        df <- data.frame(
            Position = as.integer(sapply(ptm_features, function(x) x$location$start$value)),
            Type = sapply(ptm_features, function(x) x$description)
        )
        
        df
    }
    
    ptm_df <- reactive({
        data <- protein_data()
        raw <- get_ptm_info(data)
        if (is.null(raw)) return(NULL)
        
        raw$TypeCategory <- dplyr::case_when(
            grepl("phospho", raw$Type, ignore.case = TRUE) ~ "Phosphorylation",
            grepl("acetyl", raw$Type, ignore.case = TRUE) ~ "Acetylation",
            grepl("succinyl", raw$Type, ignore.case = TRUE) ~ "Succinylation",
            grepl("methyl", raw$Type, ignore.case = TRUE) ~ "Methylation",
            TRUE ~ "Other"
        )
        
        raw <- raw |>
            dplyr::group_by(Position, TypeCategory) |>
            dplyr::summarise(tooltip = paste(paste0(Position, "：", Type), collapse = "\n"), .groups = "drop")
        
        return(raw)
    })
    
    output$ptm_control_ui <- renderUI({
        req(ptm_df())
        checkboxInput("show_ptm", "Show PTM Sites", value = FALSE)
    })
    
    output$ptm_filter_ui <- renderUI({
        df <- ptm_df()
        req(df)
        ptm_color_map <- c(
            "Phosphorylation" = "#1f77b4",
            "Acetylation"     = "#ff7f0e",
            "Succinylation"   = "#2ca02c",
            "Methylation"     = "#d62728",
            "Other"           = "#9467bd"
        )
        
        types_count <- df %>%
            count(TypeCategory) %>%
            arrange(desc(n))
        
        checkbox_tags <- lapply(seq_len(nrow(types_count)), function(i) {
            row <- types_count[i, ]
            color <- ptm_color_map[[row$TypeCategory]]
            id <- paste0("ptm_", gsub(" ", "_", row$TypeCategory))
            
            tags$div(style = "margin-bottom:1.5px;",
                     tags$label(
                         tags$input(type = "checkbox", class = "ptm_check", name = "ptm_group", value = row$TypeCategory, checked = "checked"),
                         tags$span(style = sprintf("display:inline-block;width:10px;height:10px;border-radius:50%%;background:%s;margin-right:6px;", color)),
                         paste0(row$TypeCategory, " (", row$n, ")")
                     )
            )
        })
        
        tags$div(
            tags$strong("Select PTM types to show:"),
            tags$div(id = "ptm_checkboxes", checkbox_tags),
            # JavaScript to update selected checkboxes
            tags$script(HTML("
      Shiny.onInputChange('ptm_types', Array.from(document.querySelectorAll('input.ptm_check:checked')).map(x => x.value));
      document.querySelectorAll('input.ptm_check').forEach(el => {
        el.addEventListener('change', () => {
          const values = Array.from(document.querySelectorAll('input.ptm_check:checked')).map(x => x.value);
          Shiny.setInputValue('ptm_types', values);
        });
      });
    "))
        )
    })
    
    # Consurf colors
    consurf_colors <- c(
      "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
      "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664")
    
    # 1D plot Output
    output$variant_1dplot <- renderPlotly({
        req(current_sequence())
        seq_val <- current_sequence()
        data <- if (input$seq_source == "api") protein_data() else NULL
        
        protein_len_df <- data.frame(start = 0, end = nchar(seq_val),
                                     ymin = 0.4, ymax = 0.6,
                                     label = paste0("Protein length: ", nchar(seq_val)))
        ## get the domain info
        domain_df <- NULL
        if (input$seq_source == "api" && !is.null(data$features)) {
            domains <- data$features[sapply(data$features, function(x) x$type == "Domain")]
            if (length(domains) == 0) return(NULL)
        
            domain_df <- data.frame(
                start = as.integer(sapply(domains, function(x) x$location$start$value)),
                end   = as.integer(sapply(domains, function(x) x$location$end$value)),
                description = sapply(domains, function(x) x$description)
            )
        }
        ## get missense variant info
        df <- gnomad_df()
        missense_df <- subset(df, grepl("missense_variant", Consequence, ignore.case = TRUE))
        # print(missense_df)
        missense_df <- missense_df[!is.na(missense_df$AA_Position), ]
        
        ## get PTM info
        ptm <- ptm_df()
        
        p1 <- ggplot() +
            geom_rect(data = protein_len_df, 
                      aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, text = label), 
                      fill = "grey90")
        
        if (!is.null(domain_df)) {
            p1 <- p1 +
            geom_rect(data = domain_df,
                      aes(xmin = start, xmax = end, ymin = 0.4, ymax = 0.6,
                          text = paste0("Domain: ", description, "\n", start, " - ", end)),
                      fill = "#fca6a6")
        }
        
        if (nrow(missense_df) > 0) {
            p1 <- p1 +
            geom_linerange(data = missense_df,
                           aes(x = AA_Position, ymin = 0.65, ymax = 0.95),
                           color = "slateblue", size = 0.3, alpha = 0.5)
        }
        
        if (!is.null(ptm) && nrow(ptm) > 0  && isTRUE(input$show_ptm)) {
            ptm_plot <- if (!is.null(input$ptm_types)) {
                ptm |> dplyr::filter(TypeCategory %in% input$ptm_types)
            } else {
                ptm |> dplyr::filter(FALSE) 
            }
            
            ptm_plot <- ptm |> dplyr::filter(TypeCategory %in% input$ptm_types)
            
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
                           aes(x = Position, y = 1.0, fill = TypeCategory, text = tooltip),
                           shape = 21, size = 2, color = "black", stroke = 0.3, alpha = 0.8)+
                    scale_fill_manual(values = ptm_color_map, name = "PTM Type", drop = FALSE)
            }
        }
        
        if (consurf_enabled() && !is.null(consurf_txt_df())) {
          consurf_data <- consurf_txt_df()
          consurf_data <- consurf_data %>% filter(!is.na(Position), !is.na(Grade))
          consurf_data$Color <- consurf_colors[as.character(consurf_data$Grade)]
          
          p1 <- p1 + 
            geom_linerange(data = consurf_data,
                           mapping = aes(x = Position, ymin = 0.15, ymax = 0.25),
                           inherit.aes = FALSE,
                           color = consurf_data$Color,
                           size = 1.0, alpha = 0.9)
        }
        
        p1 <- p1 +
            theme_minimal() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
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
        
        title_text <- if (input$seq_source == "api" && !is.null(data)) {
            paste0("Protein: ", data$proteinDescription$recommendedName$fullName$value,
                   "; Gene: ", gene_name())
        } else {
            "User uploaded sequence"
        }
        
        subplot(p1_plotly, p2_plotly, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
            layout(title = list(text = title_text, x = 0, xanchor = "left"),
                   margin = list(t = 60))
    })
    
    # Consurf Plot
    output$consurf_plot_ui <- renderUI({
      req(consurf_txt_df())
      tagList(
        plotOutput("consurf_txt_plot", height = "400px")
      )
    })
    
    output$consurf_txt_plot <- renderPlot({
      req(consurf_txt_df())
      df <- consurf_txt_df()
      
      ggplot(df, aes(x = Position, y = Score, fill = factor(Grade))) +
        geom_col(width = 1) +
        scale_fill_manual(values = consurf_cols, name = "ConSurf Grade") +
        theme_minimal() +
        labs(title = "ConSurf Conservation Score", x = "Residue Position", y = "Score")
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
    # Get the first PDB id
    get_first_pdb_id <- function(uniprot_id) {
        pdb_ids <- get_pdb_ids_from_uniprot(uniprot_id)
        if (length(pdb_ids) > 0) return(pdb_ids[1])
        return(NULL)
    }
    
    # Get all the PDB id in the list
    pdb_id_list <- reactive({
        req(input$pid)
        get_pdb_ids_from_uniprot(input$pid)
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
    
    # 3D structure with r3dmol
    output$structure_view <- renderR3dmol({
        req(input$structure_source)
        
        if (input$structure_source == "Upload") {
            shiny::validate(
                need(!is.null(input$pdb_upload), "Please upload a PDB file.")
            )
        } else {
            shiny::validate(
                need(input$fetch > 0, "Please click Fetch after entering UniProt ID.")
            )
        }
        
        uniprot_id <- input$pid
        pdb_data <- NULL
        
        if (input$structure_source == "AlphaFold") {
            structure_source_real("AlphaFold")
            pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",input$pid,"-F1-model_v4.pdb")))
        } else if (input$structure_source == "PDB") {
            req(input$selected_pdb_id)
            pdb_id <- input$selected_pdb_id
            
            tryCatch({
                pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://files.rcsb.org/download/", pdb_id, ".pdb")))
                structure_source_real("PDB")
            }, error = function(e){
                showNotification("No PDB structure available for this protein, fallback to AlphaFold.", type = "error")
                pdb_data <- m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",input$pid,"-F1-model_v4.pdb")))
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
        req(input$pid)
        src <- structure_source_real()
        
        if (src == "AlphaFold") {
            paste("Structure from AlphaFold for Uniprot ID:", input$pid)
        } else if(src == "PDB"){
            paste("Structure from PDB:", input$selected_pdb_id)
        } else {
                "No PDB structure found for this UniProt ID."
            }
    })
    
    # Download PDB button
    output$download_pdb <- downloadHandler(
        filename = function() {
            src <- structure_source_real()
            uniprot_id <- input$pid
            
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
            uniprot_id <- input$pid
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
    
    fells_job <- reactiveVal()
    fells_result <- reactiveVal()
    fells_class <- reactiveVal("btn-primary")
    
    output$fells_button <- renderUI({
        req(current_sequence())
        actionButton("fells", label = "Submit to FELLS", class = fells_class())
    })
    
    output$fellsAvailable <- reactive({
        fells_class() == "btn-success"
    })
    outputOptions(output, "fellsAvailable", suspendWhenHidden = FALSE)
    
    observeEvent(input$fells, {
        data <- protein_data()
        req(data)
        
        seq <- current_sequence()
        pseq <- paste0(">", input$pid, "\n", seq)
        
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
