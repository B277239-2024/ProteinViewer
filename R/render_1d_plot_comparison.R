render_1d_plot_comparison <- function(seq_val,
                                      domain_df = NULL,
                                      missense_df = NULL,
                                      ptm_df = NULL,
                                      am_df = NULL,
                                      title = NULL,
                                      fixed_x_max = NULL,
                                      domain_color_map = NULL){
  protein_len <- nchar(seq_val)
  
  p <- ggplot() +
    # protein block
    geom_rect(aes(xmin = 1, xmax = protein_len, ymin = 0.56, ymax = 0.60),
              fill = "#d3d3d3", color = NA)
  
  # Domain block
  if (!is.null(domain_df) && nrow(domain_df) > 0 && !is.null(domain_color_map)){
    domain_df$label_x <- (domain_df$start + domain_df$end) / 2
    domain_df$label_y <- 0.605
    p <- p + geom_rect(data = domain_df,
                       aes(xmin = start, xmax = end, ymin = 0.56, ymax = 0.60, fill = description),
                       color = "black", size = 0.1) +
      geom_text(data = domain_df,
                aes(x = label_x, y = label_y, label = description),
                size = 5, angle = 0, hjust = 0.5, vjust = 0, inherit.aes = FALSE) +
      scale_fill_manual(values = domain_color_map, guide = "none", aesthetics = "fill")
  }
  
  # Missense lines
  if (!is.null(missense_df) && nrow(missense_df) > 0) {
    p <- p + geom_linerange(data = missense_df,
                            aes(x = AA_Position, ymin = 0.5, ymax = 0.54),
                            color = "slateblue", alpha = 0.6, size = 0.3)
  }
  
  # PTM points
  if (!is.null(ptm_df) && nrow(ptm_df) > 0) {
    ptm_color_map <- c(
      "Phosphorylation" = "#1f77b4",
      "Acetylation"     = "#ff7f0e",
      "Succinylation"   = "#2ca02c",
      "Methylation"     = "#d62728",
      "Other"           = "#9467bd"
    )
    ptm_df$TypeCategory <- factor(ptm_df$TypeCategory, levels = names(ptm_color_map))
    
    p <- p + geom_point(data = ptm_df,
                        mapping = aes(x = Position, y = 0.65, color = TypeCategory),
                        shape = 16, size = 2, alpha = 0.8, inherit.aes = FALSE) +
      scale_color_manual(values = ptm_color_map, drop = FALSE, name = "PTM Type")
  }
  
  # AlphaMissense average scores
  if (!is.null(am_df) && nrow(am_df) > 0) {
    am_mean <- am_df %>%
      dplyr::group_by(position) %>%
      dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
    
    p <- p + geom_linerange(data = am_mean,
                            aes(x = position, ymin = 0.42, ymax = 0.42 + score * 0.06),
                            color = "darkgreen", alpha = 0.8, size = 1)
  }
  
  if (!is.null(fixed_x_max)) {
    p <- p + xlim(0, fixed_x_max)
  }
  
  p <- p + 
    coord_cartesian(ylim = c(0.38, 0.65)) +
    theme_minimal() +
    labs(title = title, x = "Residue", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 11, face = "bold", hjust = 0)
    ) +
    theme(legend.position = "none")
  
  return(p)
}