# Look at the P_Genomic and N_Genomic for the biological samples
# E. Lamont
# 1/8/26

source("Import_data.R") # To get my_pipeSummary (duplicates are removed)

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=8),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        # axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2)#
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

# Labelled Colors
my_fav_colors <- c(`W0 sputum (cure)` = "#0072B2", `W0 sputum (relapse)` = "#bc5300", `Broth`= "black", `W2 sputum (cure)` = "#0072B2", `W2 sputum (relapse)` = "#bc5300", `THP1 spiked` = "#999999")
# Labelled Shapes
my_fav_shapes <- c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21, `Broth`= 23, `W2 sputum (cure)` = 24, `W2 sputum (relapse)` = 24, `THP1 spiked` = 23)


###########################################################
################ GRAB JUST SPUTUM SAMPLES #################

sputum_pipeSummary <- my_pipeSummary %>% 
  filter(!Type2 %in% c("THP1 spiked", "Broth", "W2 sputum (failure)", "W0 sputum (failure)"))

## NOTE!: There is a treatment failure that is not being plotted here. Causes a warning message. Removed it above

###########################################################
################ N_Genomic vs SAMPLE TYPE #################

### BOXPLOT ###
N_Genomic_box1 <- sputum_pipeSummary %>% 
  ggplot(aes(x = Type, y = N_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type2, shape = Type2), alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_fill_manual(values=my_fav_colors) +  
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,12000000), breaks = seq(0, 12000000, 2000000)) +
  scale_y_continuous(trans='log10') + 
  labs(title = "N_Genomic for all sputum (Run1-4)",
    x = "Sample type", 
    y = "# reads aligning to Mtb") + 
  my_plot_themes
N_Genomic_box1
# ggsave(N_Genomic_box1,
#        file = paste0("Sputum_N_Genomic_Box_v1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 6, height = 5, units = "in")


###########################################################
################ P_Genomic vs SAMPLE TYPE #################

P_Genomic_box1 <- sputum_pipeSummary %>% 
  ggplot(aes(x = Type, y = P_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type2, shape = Type2), alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_fill_manual(values=my_fav_colors) +  
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(title = "P_Genomic for all sputum (Run1-4)",
    x = "Sample type", 
    y = "% reads aligning to Mtb") + 
  my_plot_themes
P_Genomic_box1
# ggsave(P_Genomic_box1,
#        file = paste0("Sputum_P_Genomic_Box1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 6, height = 5, units = "in")


###########################################################
############ AtLeast.10.Reads vs SAMPLE TYPE ##############

TenReads_box2 <- sputum_pipeSummary %>% 
  ggplot(aes(x = Type, y = Txn_Coverage_f)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type2, shape = Type2), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  scale_fill_manual(values=my_fav_colors) +  
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100.2), breaks = seq(0, 100.2, 10)) + 
  geom_hline(yintercept = 60, linetype = "dashed", alpha = 0.5) + 
  # geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(title = "Genes with >= 10 reads aligning for all sputum (Run1-4)",
    # subtitle = "", 
    x = "Sample type", 
    y = "% transcriptional coverage") + 
  my_plot_themes
TenReads_box2
# ggsave(TenReads_box2,
#        file = paste0("Sputum_TenReads_box2.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 6, height = 5, units = "in")




