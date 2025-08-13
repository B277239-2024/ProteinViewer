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
    gene_name = NULL,
    fells_enabled = FALSE,
    fells_result = NULL,
    return_plot_only = FALSE
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
                       fill = "#d3d3d3")
  
  if (!is.null(domain_df) && nrow(domain_df) > 0){
    p1 <- p1 +
      ggplot2::geom_rect(data = domain_df,
                         aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55,
                             text = paste0("Domain: ", description, "\n", start, " - ", end)),
                         fill = "#fca6a6")
    domain_df$label_x <- (domain_df$start + domain_df$end) / 2
    domain_df$label_y <- 0.56 
    p1 <- p1 +
      ggplot2::geom_text(
        data = domain_df,
        aes(x = label_x, y = label_y, label = description),
        size = 3, vjust = 0, inherit.aes = FALSE
      )
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
                                            text = paste0("ConSurf\n", "Position: ", Position, "\n",
                                                          "Score: ", round(Score, 3), "\n", "Grade: ", Grade)),
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
            text = paste0("AlphaMissense\nPosition: ", position, "\nAvg Score: ", round(avg_score, 3))
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
  
  # FELLS subplots
  p3 <- p4 <- NULL
  if (fells_enabled && !is.null(fells_result)) {
    # Secondary structure
    hsc <- fells_result$p_h |> unlist()
    est <- fells_result$p_e |> unlist()
    coil <- fells_result$p_c |> unlist()
    df3 <- data.frame(index = 1:length(hsc),
                      Helix = hsc,
                      Strand = est,
                      Coil = coil) |>
      tidyr::pivot_longer(cols = c("Helix", "Strand", "Coil")) |>
      dplyr::mutate(
        value = as.numeric(value),
        value = dplyr::if_else(name == "Coil", -value, value))
    
    p3 <- ggplot2::ggplot(df3, aes(x = index, y = value, fill = name)) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::scale_fill_manual(values = c("Helix" = "#91288c", "Strand" = "#ffa500", "Coil" = "gray50")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title.x = element_blank()) +
      ggplot2::ggtitle("Secondary Structure Prediction") +
      ggplot2::labs(y = "Prob.")
    
    # Disorder / HCA
    dis <- fells_result$p_dis |> unlist()
    hca <- fells_result$hca |> unlist()
    df4 <- data.frame(index = 1:length(dis),
                      Disorder = dis,
                      HCA = hca) |>
      tidyr::pivot_longer(cols = c("Disorder", "HCA")) |>
      dplyr::mutate(
        value = as.numeric(value),
        value = dplyr::if_else(name == "Disorder", -value, value))
    
    p4 <- ggplot2::ggplot(df4, aes(x = index, y = value, fill = name)) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::scale_fill_manual(values = c("Disorder" = "red", "HCA" = "black")) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Disorder & HCA") +
      ggplot2::labs(y = "Score (±)")
  }
  
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
  
  p1_plot <- plotly::ggplotly(p1, tooltip = "text") %>%
    plotly::layout(margin = list(b = 0))
  
  p2_plot <- plotly::ggplotly(p2) %>%
    plotly::layout(margin = list(t = 0))
  
  p3_plot <- if (!is.null(p3)) {
    plotly::ggplotly(p3) %>% plotly::layout(margin = list(t = 0))
  } else NULL
  
  p4_plot <- if (!is.null(p4)) {
    plotly::ggplotly(p4) %>% plotly::layout(margin = list(t = 0))
  } else NULL
  
  plots <- list(p1_plot)
  if (!is.null(p3_plot)) plots <- append(plots, list(p3_plot))
  if (!is.null(p4_plot)) plots <- append(plots, list(p4_plot))
  plots <- append(plots, list(p2_plot))
  
  has_p3 <- !is.null(p3_plot)
  has_p4 <- !is.null(p4_plot)
  
  subplot_heights <- c(2)
  
  if (has_p3) subplot_heights <- c(subplot_heights, 1)
  if (has_p4) subplot_heights <- c(subplot_heights, 1)
  
  subplot_heights <- c(subplot_heights, 1)
  subplot_heights <- subplot_heights / sum(subplot_heights)
  
  title_text <- if (!is.null(protein_data)) {
    paste0("Protein: ", protein_data$proteinDescription$recommendedName$fullName$value,
           "; Gene: ", gene_name)
  } else {
    "User uploaded sequence"
  }
  
  if (return_plot_only) {
    combined_plot <- p1
    if (!is.null(p3)) combined_plot <- combined_plot / p3
    if (!is.null(p4)) combined_plot <- combined_plot / p4
    combined_plot <- combined_plot / p2
    return(combined_plot)
  }
  
  plotly::subplot(plots,
                  nrows = length(plots),
                  shareX = TRUE,
                  titleY = TRUE,
                  heights = subplot_heights) %>%
    plotly::layout(
      title = list(text = title_text, x = 0, xanchor = "left"),
      margin = list(t = 60)
    )
}