
makeVolcano_EL <- function(df, title = "Volcano Plot") {
  
  ## Taking from DEG values with extra columns added (Function_Add_DE_Columns.R)
  
  my_plot_themes <- theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none",legend.text=element_text(size=14),
          legend.title = element_text(size = 14),
          plot.title = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=14), 
          plot.subtitle = element_text(size=10))
  
  my_volcano <- df %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE1, label = DE1_labels, text = gene)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = title) + 
    geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
  
  plot_build <- ggplot_build(my_volcano)
  
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  text_up <- df %>% filter(DE1 == "significant up") %>% nrow()
  text_down <- df %>% filter(DE1 == "significant down") %>% nrow()
  
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
  
  final_volcano <- my_volcano_annotated + my_plot_themes
  
  return(final_volcano)
}





