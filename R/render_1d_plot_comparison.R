render_1d_plot_comparison <- function(seq_val,
                                      domain_df = NULL,
                                      missense_df = NULL,
                                      ptm_df = NULL,
                                      am_df = NULL,
                                      title = NULL) {
  protein_len <- nchar(seq_val)
  
  p <- ggplot() +
    # protein block
    geom_rect(aes(xmin = 1, xmax = protein_len, ymin = 0.45, ymax = 0.55),
              fill = "#d3d3d3", color = NA)
  
  # Domain block
  if (!is.null(domain_df) && nrow(domain_df) > 0) {
    p <- p + geom_rect(data = domain_df,
                       aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55),
                       fill = "#fca6a6", color = "black", size = 0.1)
  }
  
  # Missense lines
  if (!is.null(missense_df) && nrow(missense_df) > 0) {
    p <- p + geom_linerange(data = missense_df,
                            aes(x = AA_Position, ymin = 0.6, ymax = 0.7),
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
    p <- p + geom_point(data = ptm_df,
                        aes(x = Position, y = 0.75, fill = TypeCategory),
                        shape = 21, size = 2, stroke = 0.3, color = "black", alpha = 0.8) +
      scale_fill_manual(values = ptm_color_map, drop = FALSE)
  }
  
  # AlphaMissense average scores
  if (!is.null(am_df) && nrow(am_df) > 0) {
    am_mean <- am_df %>%
      dplyr::group_by(position) %>%
      dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
    
    p <- p + geom_linerange(data = am_mean,
                            aes(x = position, ymin = 0.15, ymax = 0.15 + score * 0.1),
                            color = "darkgreen", alpha = 0.8, size = 1)
  }
  
  p + theme_minimal() +
    labs(title = title, x = "Residue", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0)
    )
}