# Boxplot of TPM for single genes
# E. Lamont
# 3/12/26

# Hassan is interest in Rv2557 so start with that one

# source("Import_data.R")


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=8),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        # axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9))

facet_themes <- theme(strip.background=element_rect(fill="lightgrey", linewidth = 0.9),
                      strip.text = element_text(size = 10))

# Labelled Colors
my_fav_colors <- c(`Cure` = "#0072B2", `Relapse` = "#bc5300")
# Labelled Shapes
# my_fav_shapes <- c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21, `Broth`= 23, `W2 sputum (cure)` = 24, `W2 sputum (relapse)` = 24, `THP1 spiked` = 23)

###########################################################
######################  CHOOSE GENE  ######################

# myGene <- "Rv2557"
myGene <- "Rv2703"

###########################################################
######################  PROCESS DATA  #####################

my_tpmf <- GoodSamples60_tpmf %>% 
  dplyr::select(contains("W")) %>% 
  rownames_to_column("gene") %>%
  filter(gene == myGene) %>%
  pivot_longer(cols = contains("W"), names_to = "SampleID2", values_to = "TPM") %>%
  left_join(GoodSputum60_pipeSummary %>% dplyr::select(SampleID2, Outcome, Week), 
            by = "SampleID2")

###########################################################
########################  BOXPLOT  ########################

my_TPM_Boxplot <- my_tpmf %>%
  ggplot(aes(x = Outcome, y = TPM)) +
  # ggplot(aes(x = Outcome, y = log2(TPM+1))) +
  geom_boxplot(aes(fill = Outcome), width = 0.6, outlier.size = 0.9, alpha = 0.4) +
  geom_point(alpha = 0.8, size = 1.7, position = position_jitter(0.2)) +
  facet_grid(~ Week, scales = "free") +
  stat_compare_means(aes(group = Outcome), method = "t.test", paired = F, label = "p.format", size = 3, label.x = 1.4) +
  scale_fill_manual(values=my_fav_colors) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  labs(title = myGene,
       subtitle = "Two-sample t-test. Mean is red diamond",
       x = "Oucome",
       y = "TPM") +
  my_plot_themes + facet_themes
my_TPM_Boxplot
# ggsave(my_TPM_Boxplot,
#        file = paste0(myGene, "_v1.pdf"),
#        path = "Figures/Gene_Boxplots",
#        width = 6, height = 4, units = "in")

###########################################################
#####################  BOXPLOT LOG2  ######################

my_log2.TPM_Boxplot <- my_tpmf %>% 
  ggplot(aes(x = Outcome, y = log2(TPM+1))) + 
  geom_boxplot(aes(fill = Outcome), width = 0.6, outlier.size = 0.9, alpha = 0.4) + 
  geom_point(alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  facet_grid(~ Week, scales = "free") +
  stat_compare_means(aes(group = Outcome), method = "t.test", paired = F, label = "p.format", size = 3, label.x = 1.4) + 
  scale_fill_manual(values=my_fav_colors) +  
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") + 
  labs(title = myGene,
       subtitle = "Two-sample t-test. Mean is red diamond", 
       x = "Oucome", 
       y = "log2(TPM+1)") + 
  my_plot_themes + facet_themes
my_log2.TPM_Boxplot
ggsave(my_log2.TPM_Boxplot,
       file = paste0(myGene, "_log2.pdf"),
       path = "Figures/Gene_Boxplots",
       width = 6, height = 4, units = "in")




















