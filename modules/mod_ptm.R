# modules/mod_ptm.R

mod_ptm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("ptm_control_ui")),
    uiOutput(ns("ptm_filter_ui"))
  )
}

mod_ptm_server <- function(id, protein_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 提取 PTM 信息
    get_ptm_info <- function(protein_json) {
      if (is.null(protein_json$features)) return(NULL)
      ptm_features <- protein_json$features[sapply(protein_json$features, function(x) x$type == "Modified residue")]
      if (length(ptm_features) == 0) return(NULL)
      
      df <- data.frame(
        Position = as.integer(sapply(ptm_features, function(x) x$location$start$value)),
        Type = sapply(ptm_features, function(x) x$description)
        #stringsAsFactors = FALSE
      )
      return(df)
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
        dplyr::summarise(
          tooltip = paste(paste0(Position, "：", Type), collapse = "\n"), 
          .groups = "drop"
        )
      
      return(raw)
    })
    
    output$ptm_control_ui <- renderUI({
      req(ptm_df())
      checkboxInput(ns("show_ptm"), "Show PTM Sites", value = FALSE)
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
        dplyr::count(TypeCategory) %>%
        dplyr::arrange(desc(n))
      
      checkbox_tags <- lapply(seq_len(nrow(types_count)), function(i) {
        row <- types_count[i, ]
        color <- ptm_color_map[[row$TypeCategory]]
        id <- paste0("ptm_", gsub(" ", "_", row$TypeCategory))
        tags$div(style = "margin-bottom:1.5px;",
                 tags$label(
                   tags$input(
                     type = "checkbox", 
                     class = "ptm_check", 
                     name = "ptm_group", 
                     value = row$TypeCategory, 
                     checked = "checked"
                   ),
                   tags$span(style = sprintf("display:inline-block;width:10px;height:10px;border-radius:50%%;background:%s;margin-right:6px;", color)),
                   paste0(row$TypeCategory, " (", row$n, ")")
                 )
        )
      })
      
      tags$div(
        tags$strong("Select PTM types to show:"),
        tags$div(id = ns("ptm_checkboxes"), checkbox_tags),
        tags$script(HTML(sprintf("
          Shiny.onInputChange('%s', Array.from(document.querySelectorAll('#%s input.ptm_check:checked')).map(x => x.value));
          document.querySelectorAll('#%s input.ptm_check').forEach(el => {
            el.addEventListener('change', () => {
              const values = Array.from(document.querySelectorAll('#%s input.ptm_check:checked')).map(x => x.value);
              Shiny.setInputValue('%s', values);
            });
          });
        ", ns("ptm_types"), ns("ptm_checkboxes"), ns("ptm_checkboxes"), ns("ptm_checkboxes"), ns("ptm_types"))))
      )
    })
    
    return(list(
      ptm_df = ptm_df,
      show = reactive({ input$show_ptm }),
      selected_types = reactive({ input$ptm_types })
    ))
  })
}