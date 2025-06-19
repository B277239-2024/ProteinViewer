library(shiny)
library(httr)
library(jsonlite)
library(DT)
library(r3dmol)
library(shinycssloaders)
library(bio3d)

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
fetch_gnomad_by_gene <- function(gene_symbol) {
    url <- "https://gnomad.broadinstitute.org/api"
    query <- '
    query VariantsInGene($geneSymbol: String!) {
        gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            variants(dataset: gnomad_r4) {
                    variant_id
                    pos
                    consequence
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
            geneSymbol = gene_symbol
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
    
    # Basic info output
    output$info <- renderPrint({
        data <- protein_data()
        if (!is.null(data)){
            list(
                ID = data$primaryAccession,
                Protein_Name = data$proteinDescription$recommendedName$fullName$value,
                Gene = data$genes[[1]]$geneName$value,
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
    
    # Extract GnomAD variant dataframe
    get_gnomad_df <- function(data) {
        gene_symbol <- tryCatch(data$genes[[1]]$geneName$value, error = function(e) NULL)
        if (is.null(gene_symbol)) return(NULL)
        
        gnomad_data <- fetch_gnomad_by_gene(gene_symbol)
        if (is.null(gnomad_data)) return(NULL)
        
        vars <- tryCatch(gnomad_data$data$gene$variants, error = function(e) NULL)
        if (is.null(vars)) return(NULL)
        
        df <- data.frame(
            Variant_ID  = vars$variant_id,
            Position    = vars$pos,
            Consequence = vars$consequence,
            AF = sapply(vars$exome, function(x) {
                if (is.null(x) || !is.list(x)) return(NA_real_)
                if (!is.null(x$af)) x$af else NA_real_
            }),
            AC = sapply(vars$exome, function(x) {
                if (is.null(x) || !is.list(x)) return(NA_integer_)
                if (!is.null(x$ac)) x$ac else NA_integer_
            }),
            AN = sapply(vars$exome, function(x) {
                if (is.null(x) || !is.list(x)) return(NA_integer_)
                if (!is.null(x$an)) x$an else NA_integer_
            })
        )
        
        df
    }
    
    # GnomAD Summary info Output
    output$gnomad_summary <- renderUI({
        data <- protein_data()
        df <- get_gnomad_df(data)
        
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
