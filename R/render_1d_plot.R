render_1d_plot <- function(
    seq_val,
    protein_data = NULL,
    missense_df,
    domain_df = NULL,
    ptm_df = NULL,
    ptm_show = FALSE,
    ptm_types = NULL,
    consurf_df = NULL,
    consurf_enabled = FALSE,
    am_df = NULL,
    am_enabled = FALSE,
    vdvp_window = 3,
    gene_name = NULL
) {
  consurf_colors <- c(
    "1" = "#10C8D2", "2" = "#89FDFD", "3" = "#D8FDFE", "4" = "#EAFFFF",
    "5" = "#FFFFFF", "6" = "#FBECF1", "7" = "#FAC9DE", "8" = "#F27EAB", "9" = "#A22664"
  )
  
  protein_len_df <- data.frame(start = 0, end = nchar(seq_val),
                               ymin = 0.45, ymax = 0.55,
                               label = paste0("Protein length: ", nchar(seq_val)))
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = protein_len_df, 
                       aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, text = label), 
                       fill = "grey90")
  
  if (!is.null(domain_df) && nrow(domain_df) > 0){
    p1 <- p1 +
      ggplot2::geom_rect(data = domain_df,
                         aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55,
                             text = paste0("Domain: ", description, "\n", start, " - ", end)),
                         fill = "#fca6a6")
  }
  
  if (nrow(missense_df) > 0) {
    p1 <- p1 +
      ggplot2::geom_linerange(data = missense_df,
                              aes(x = AA_Position, ymin = 0.6, ymax = 0.7),
                              color = "slateblue", size = 0.3, alpha = 0.5)
  }
  
  if (!is.null(ptm_df) && nrow(ptm_df) > 0 && isTRUE(ptm_show)) {
    ptm_plot <- if (!is.null(ptm_types)) {
      dplyr::filter(ptm_df, TypeCategory %in% ptm_types)
    } else {
      ptm_df[FALSE, ]
    }
    
    if (nrow(ptm_plot) > 0) {
      ptm_color_map <- c(
        "Phosphorylation" = "#1f77b4",
        "Acetylation" = "#ff7f0e",
        "Succinylation" = "#2ca02c",
        "Methylation" = "#d62728",
        "Other" = "#9467bd"
      )
      
      p1 <- p1 +
        ggplot2::geom_point(data = ptm_plot,
                            aes(x = Position, y = 0.75, fill = TypeCategory, text = tooltip),
                            shape = 21, size = 2, color = "black", stroke = 0.3, alpha = 0.8) +
        ggplot2::scale_fill_manual(values = ptm_color_map, name = "PTM Type", drop = FALSE)
    }
  }
  
  if (consurf_enabled && !is.null(consurf_df)) {
    consurf_data <- consurf_df %>% dplyr::filter(!is.na(Position), !is.na(Grade))
    consurf_data$Color <- consurf_colors[as.character(consurf_data$Grade)]
    
    p1 <- p1 + 
      ggplot2::geom_linerange(data = consurf_data,
                              mapping = aes(x = Position, ymin = 0.3, ymax = 0.4, 
                                            text = paste0("ConSurf\\n", "Position: ", Position, "\\n",
                                                          "Score: ", round(Score, 3), "\\n", "Grade: ", Grade)),
                              inherit.aes = FALSE,
                              color = consurf_data$Color,
                              size = 1.0, alpha = 0.9)
  }
  
  if (am_enabled && !is.null(am_df)) {
    alpha_bar_df <- am_df %>%
      dplyr::group_by(position) %>%
      dplyr::summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop")
    
    p1 <- p1 +
      ggplot2::geom_linerange(
        data = alpha_bar_df,
        aes(x = position,
            ymin = 0.15,
            ymax = avg_score * 0.1 + 0.15,
            text = paste0("AlphaMissense\\nPosition: ", position, "\\nAvg Score: ", round(avg_score, 3))
        ),
        color = "darkgreen",
        size = 1.0,
        alpha = 0.8)
  }
  
  p1 <- p1 +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10)
    ) +
    ggplot2::labs(x = NULL)
  
  mut_index <- missense_df$AA_Position |> unique()
  prot_len <- nchar(seq_val)
  window <- if (vdvp_window < 1) max(1, floor(prot_len * vdvp_window)) else as.integer(vdvp_window)
  vp <- length(mut_index) / prot_len
  index <- 1:prot_len
  vdvp_vals <- sapply(index, function(x) {
    vd <- sum(mut_index %in% x:(x + window))
    vd / window / vp
  })
  vdvp_df <- data.frame(x = index + window / 2, y = vdvp_vals)
  vdvp_df <- as.data.frame(spline(vdvp_df$x, vdvp_df$y))
  
  p2 <- ggplot2::ggplot(vdvp_df, aes(x = x, y = y)) +
    ggplot2::geom_line(color = "steelblue", size = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Residue", y = "Vd/Vp") +
    ggplot2::theme(plot.title = ggplot2::element_blank())
  
  title_text <- if (!is.null(protein_data)) {
    paste0("Protein: ", protein_data$proteinDescription$recommendedName$fullName$value,
           "; Gene: ", gene_name)
  } else {
    "User uploaded sequence"
  }
  
  has_extra_layer <- consurf_enabled || am_enabled
  subplot_heights <- if (has_extra_layer) c(0.7, 0.3) else c(0.5, 0.5)
  
  plotly::subplot(
    plotly::ggplotly(p1, tooltip = "text") %>% plotly::layout(margin = list(b = 0)),
    plotly::ggplotly(p2) %>% plotly::layout(margin = list(t = 0)),
    nrows = 2, shareX = TRUE, titleY = TRUE, heights = subplot_heights
  ) %>%
    plotly::layout(title = list(text = title_text, x = 0, xanchor = "left"),
                   margin = list(t = 60))
}