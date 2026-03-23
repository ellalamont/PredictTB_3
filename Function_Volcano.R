
makeVolcano_EL <- function(my_df, graph_title, DE_limit) {
  
  ## Make a volcano plot using output from Bob's pipeline
  ## Original p-values
  ## DE_limit can be 1, 2, 2.5
  
  Current_plot_themes <- theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none",legend.text=element_text(size=10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=14), 
          plot.subtitle = element_text(size=10), 
          # panel.background = element_rect(fill='transparent'),
          # plot.background = element_rect(fill='transparent', color=NA),
          # legend.background = element_rect(fill='transparent'),
          legend.box.background = element_blank())
  
  if (DE_limit == 1) {
    my_DE_col <- "DE1"
    my_DE_label <- "DE1_labels"
  } else if (DE_limit == 2) {
    my_DE_col <- "DE2"
    my_DE_label <- "DE2_labels"
  } else if (DE_limit == 2.5) {
    my_DE_col <- "DE2.5"
    my_DE_label <- "DE2.5_labels"
  }
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = .data[[my_DE_col]], label = .data[[my_DE_label]], text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = paste0(graph_title, " Log2Fold=", DE_limit, " AVG p-value")) + 
    geom_vline(xintercept = c(-DE_limit,DE_limit), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel # Changed from 3 to 4
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) # +
  
  
  # Determine the max and min axes values for labeling 
  plot_build <- ggplot_build(my_volcano)
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  # Add the gene number annotations
  text_up <- my_df %>% filter(.data[[my_DE_col]] == "significant up") %>% nrow()
  text_down <- my_df %>% filter(.data[[my_DE_col]] == "significant down") %>% nrow()
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "white") + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "white")
  

  final_volcano <- my_volcano_annotated + 
    scale_x_continuous(breaks = seq(floor(x_min), ceiling(x_max), by = 1)) +
    Current_plot_themes
  
  return(final_volcano)
}











# OLD BELOW
# makeVolcano_EL <- function(df, title = "Volcano Plot") {
#   
#   ## Taking from DEG values with extra columns added (Function_Add_DE_Columns.R)
#   
#   my_plot_themes <- theme_bw() +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     theme(legend.position = "none",legend.text=element_text(size=14),
#           legend.title = element_text(size = 14),
#           plot.title = element_text(size=10), 
#           axis.title.x = element_text(size=14), 
#           axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
#           axis.title.y = element_text(size=14),
#           axis.text.y = element_text(size=14), 
#           plot.subtitle = element_text(size=10))
#   
#   my_volcano <- df %>%
#     ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE1, label = DE1_labels, text = gene)) + # text is for plotly, could be GENE_ID
#     geom_point(alpha = 0.7) + 
#     labs(title = title) + 
#     geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
#     geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
#     geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
#     scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
#   
#   plot_build <- ggplot_build(my_volcano)
#   
#   y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
#   x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
#   x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
#   
#   text_up <- df %>% filter(DE1 == "significant up") %>% nrow()
#   text_down <- df %>% filter(DE1 == "significant down") %>% nrow()
#   
#   my_volcano_annotated <- my_volcano +
#     annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
#     annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
#   
#   final_volcano <- my_volcano_annotated + my_plot_themes
#   
#   return(final_volcano)
# }
# 
# 
# 
# 
# 
