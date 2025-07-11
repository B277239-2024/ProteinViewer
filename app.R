library(shiny)
library(httr)
library(jsonlite)
library(DT)
library(r3dmol)
library(shinycssloaders)
library(bio3d)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(plotly)

# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Protein info Viewer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            textInput("pid", "Enter UniProt ID"),
            actionButton("fetch", "Fetch info"), 
            hr(),
            
            uiOutput("transcript_selector"),
            
            # Select the style of 3D structure
            selectInput("set_style", "Choose Structure Style",
                        choices = c("Cartoon", "Line", "Stick", "Sphere", "Cross"),
                        selected = "Cartoon"),
            
            # 3D structure spin
            checkboxInput("spin", "Spin Structure", value = FALSE),
            checkboxInput("surface", "Show Surface", value = FALSE),
            # show labels
            # checkboxInput("labels", "Show Labels on Selected", value = FALSE),
            # highlight selected residues
            # actionButton("selectSpheres", "Highlight Selected Residues")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           #h4("Protein Info"),
           #verbatimTextOutput("info"),
           #h4("Sequence"),
           #verbatimTextOutput("sequence"),
           #h4("Domains"),
           #tableOutput("domain_table"), 
           #h4("GnomAD Summary"),
           #uiOutput("gnomad_summary"),
           #h4("AlphaFold 3D Structure"),
           #downloadButton("download_pdb", "Download PDB File")
            uiOutput("main_content")
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
    output$main_content <- renderUI({
        data <- protein_data()
        req(data)
        
        tagList(
            h4("Protein Info"),
            verbatimTextOutput("info"),
            h4("Sequence"),
            verbatimTextOutput("sequence"),
            h4("Domains"),
            tableOutput("domain_table"),
            h4("GnomAD Summary"),
            uiOutput("gnomad_summary"),
            h4("GnomAD Variant Table"),
            DT::dataTableOutput("gnomad_table"),
            h4("1D Plot"),
            plotlyOutput("variant_1dplot", height="600px"),
            h4("AlphaFold 3D Structure"),
            withSpinner(r3dmolOutput("structure_view", height = "500px")),
            br(),
            downloadButton("download_pdb", "Download PDB File")
        )
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
    
    # Get the transcripts ID
    output$transcript_selector <- renderUI({
        tx_df <- transcript_table()
        if (nrow(tx_df) == 0) return(helpText("No transcript IDs found."))
        
        choices <- setNames(tx_df$ensembl_transcript_id, 
                            paste0(tx_df$ensembl_transcript_id, " (", tx_df$external_gene_name, ")"))
        
        selectInput("selected_transcript", "Select Transcript", choices = choices)
    })
    
    # 响应式封装GnomAD data获取
    gnomad_df <- reactive({
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
    
    # Basic info output
    output$info <- renderPrint({
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
    })
    
    # Protein sequence output
    output$sequence <- renderText({
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
    
    # 1D plot Output
    output$variant_1dplot <- renderPlotly({
        data <- protein_data()
        req(data)
        
        protein_len_df <- data.frame(start = 0, end = data$sequence$length,
                                     ymin = 0.4, ymax = 0.6,
                                     label = paste0("Protein length: ", data$sequence$length))
        ## get the domain info
        domains <- data$features[sapply(data$features, function(x) x$type == "Domain")]
        if (length(domains) == 0) return(NULL)
        
        domain_df <- data.frame(
            start = as.integer(sapply(domains, function(x) x$location$start$value)),
            end   = as.integer(sapply(domains, function(x) x$location$end$value)),
            description = sapply(domains, function(x) x$description)
        )
        ## get missense variant info
        df <- gnomad_df()
        missense_df <- subset(df, grepl("missense_variant", Consequence, ignore.case = TRUE))
        # print(missense_df)
        missense_df <- missense_df[!is.na(missense_df$AA_Position), ]
        
        p1 <- ggplot() +
            geom_rect(data = protein_len_df, 
                      aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, text = label), 
                      fill = "grey90") +
            geom_rect(data = domain_df,
                      aes(xmin = start, xmax = end, ymin = 0.4, ymax = 0.6,
                          text = paste0("Domain: ", description, "\n", start, " - ", end)),
                      fill = "#fca6a6") +
            geom_linerange(data = missense_df,
                           aes(x = AA_Position, ymin = 0.65, ymax = 0.95),
                           color = "slateblue", size = 0.3, alpha = 0.5) +
            theme_minimal() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none"
            ) +
            labs(x = NULL)
        
        #bin_width <- 3
        #missense_df$bin <- cut(missense_df$AA_Position,
         #                      breaks = seq(0, max(missense_df$AA_Position, na.rm = TRUE) + bin_width, by = bin_width),
          #                     right = FALSE)
        #bin_counts <- missense_df |>
         #   group_by(bin) |>
          #  summarise(
           #     start = as.numeric(gsub("\\[|,.*", "", bin)),  # 提取起点
            #    Vd = n()
            #)
        #bin_counts$Vp <- bin_width
        #bin_counts$y <- bin_counts$Vd / bin_counts$Vp
        #vd_df <- bin_counts
        mut_index <- missense_df$AA_Position |> unique()
        prot_len <- data$sequence$length
        window <- 3
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
        
        p1_plotly <- ggplotly(p1, tooltip = "text") %>% layout(margin = list(b = 0))
        p2_plotly <- ggplotly(p2) %>% layout(margin = list(t = 0))
        
        protein_name <- data$proteinDescription$recommendedName$fullName$value
        subplot(p1_plotly, p2_plotly, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
            layout(title = list(text = paste0("Protein: ", protein_name, "; Gene: ", gene_name()), x = 0, xanchor = "left"),
                   margin = list(t = 60))
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
    
    # 3D structure with r3dmol
    output$structure_view <- renderR3dmol({
        req(input$pid)
        
        shiny::validate(
            need(nzchar(input$pid), "Please enter a UniProt ID.")
        )
        
        pdb_url <- m_bio3d(bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",input$pid,"-F1-model_v4.pdb")))
        
        r3dmol(
            viewer_spec = m_viewer_spec(backgroundColor = "white", 
                                        cartoonQuality = 25,
                                        lowerZoomLimit = 5,
                                        upperZoomLimit = 1000)
        ) %>%
            m_add_model(data = pdb_url, format = "pdb") %>%
            m_set_style(style = m_style_cartoon(color = "spectrum")) %>%
            m_zoom_to()
    })
    
    
    
    # Download PDB button
    output$download_pdb <- downloadHandler(
        filename = function() {
            uniprot_id <- input$pid
            af_data <- fetch_alphafold_prediction(uniprot_id)
            
            if (!is.null(af_data) && !is.null(af_data[[1]]$entryId)) {
                paste0(af_data[[1]]$entryId, ".pdb")
            } else {
                paste0("AF-", uniprot_id, "-model.pdb")
            }
        },
        
        content = function(file) {
            uniprot_id <- input$pid
            af_data <- fetch_alphafold_prediction(uniprot_id)
            
            if (is.null(af_data) || length(af_data) == 0 || is.null(af_data[[1]]$pdbUrl)) {
                warning("AlphaFold model not available for this protein.")
                return(NULL)
            }
            
            pdb_url <- af_data[[1]]$pdbUrl
            download.file(pdb_url, destfile = file, mode = "wb")
        }
    )
    
    observeEvent(input$set_style, {
        style <- switch(input$set_style,
                        "Line" = list(line = list()),
                        "Cartoon" = list(cartoon = list()),
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
}

# Run the application 
shinyApp(ui = ui, server = server)
