makeBubble_EL <- function(df, myTitle = "Bubble Plot") {
  
  ## Taking a modified MetaGeneSet file from Bob's pipeline and making a bubble from it
  
  my_plot_themes <- theme_bw() +
    theme(legend.position = "right",legend.text=element_text(size=10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size=10), 
          axis.title.x = element_text(size=10), 
          axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
          # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
          axis.title.y = element_text(size=10),
          axis.text.y = element_text(size=10), 
          plot.subtitle = element_text(size=10))
  
  facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                        strip.text = element_text(size = 7))
  
  my_bubblePlot <- myBubble_df %>% 
    ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
    geom_point(aes(stroke = ifelse(Significance == "significant", 0.8, 0),
                   fill = FillColor),
               size = 4, shape = 21, alpha = 0.9) + 
    scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
    facet_grid(rows = vars(Group_wrapped), scales = "free_y", space = "free") + 
    guides(shape = "none") + 
    geom_vline(xintercept = 0) + 
    labs(title = myTitle, 
         y = NULL, x = "Log2Fold change") + 
    my_plot_themes + facet_themes + theme(legend.position = "none")
  
  return(my_bubblePlot)
  
  
  
  
  
}