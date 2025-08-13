render_fells_plot <- function(hsc_df, dh_df, title_1 = "Secondary Structure Prediction", title_2 = "Disorder & HCA", use_alpha = TRUE) {
  hsc_cols <- c(Helix = "#91288c", Strand = "#ffa500", Coil = "gray50")
  dh_cols  <- c(Disorder = "red", HCA = "black")
  
  p1 <- ggplot(hsc_df, aes(index, value, fill = name, alpha = if (use_alpha) alpha else NULL)) +
    geom_col(position = "identity") +
    scale_fill_manual(values = hsc_cols) +
    guides(alpha = "none") +
    theme_minimal() +
    ggtitle(title_1) +
    labs(y = "Secondary structure probability")
  
  p2 <- ggplot(dh_df, aes(index, value, fill = name, alpha = if (use_alpha) alpha else NULL)) +
    geom_col(position = "identity") +
    scale_fill_manual(values = dh_cols) +
    guides(alpha = "none") +
    theme_minimal() +
    ggtitle(title_2) +
    labs(y = "Disorder/HCA (scaled)")
  
  patchwork::wrap_plots(p1, p2, ncol = 1)
}