# Boxplot of TPM for Gene Sets
# E. Lamont
# 3/20/26

# source("Import_data.R")


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=8),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9))

facet_themes <- theme(strip.background=element_rect(fill="lightgrey", linewidth = 0.9),
                      strip.text = element_text(size = 10))

# Labelled Colors
my_fav_colors <- c(`Cure` = "#0072B2", `Relapse` = "#bc5300")

###########################################################
####################  CHOOSE GENESET  #####################

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))
Current_GeneSet_df %>% select(GeneSet) %>% unique()

# Needs to match the GeneSets names
myGeneSet <- "Enduring hypoxic response"
# myGeneSet <- "dosR regulon"
# myGeneSet <- "Non-specific stress responses"
myGenes <- Current_GeneSet_df %>% filter(GeneSet == myGeneSet) %>% pull(Gene)

###########################################################
######################  PROCESS DATA  #####################

my_tpmf <- GoodSamples60_tpmf %>% 
  dplyr::select(contains("W")) %>% 
  rownames_to_column("gene") %>%
  filter(gene %in% myGenes) %>%
  pivot_longer(cols = contains("W"), names_to = "SampleID2", values_to = "TPM") %>%
  left_join(GoodSputum60_pipeSummary %>% dplyr::select(SampleID2, Outcome, Week), 
            by = "SampleID2")

my_tpmf_avg <- my_tpmf %>%
  group_by(gene, Week, Outcome) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), 
            SD_TPM = sd(TPM, na.rm = TRUE),
            N_Samples = n(),
            .groups = "drop")

my_Log2_tpmf_avg <- my_tpmf %>%
  mutate(Log2_TPM = log2(TPM+1)) %>%
  group_by(gene, Week, Outcome) %>%
  summarise(mean_Log2_TPM = mean(Log2_TPM, na.rm = TRUE), 
            SD_Log2_TPM = sd(Log2_TPM, na.rm = TRUE),
            N_Samples = n(),
            .groups = "drop")


###########################################################
########################  BOXPLOT  ########################

my_TPM_Boxplot <- my_tpmf_avg %>%
  ggplot(aes(x = Week, y = mean_TPM)) +
  geom_boxplot(aes(fill = Outcome), width = 0.6, outlier.size = 0.9, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1, position = position_jitter(0.2)) +
  facet_grid(~ Outcome, scales = "free") +
  stat_compare_means(aes(group = Week), method = "t.test", paired = F, label = "p.format", size = 3, label.x = 1.4) +
  scale_fill_manual(values=my_fav_colors) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  labs(title = paste0(myGeneSet, ". Each dot is the average value of one gene"),
       subtitle = "Two-sample t-test. Mean is red diamond",
       y = "TPM") +
  my_plot_themes + facet_themes
my_TPM_Boxplot
ggsave(my_TPM_Boxplot,
       file = paste0(myGeneSet, "_tpmf.pdf"),
       path = "Figures/GeneSet_Boxplots",
       width = 7, height = 5, units = "in")

###########################################################
#####################  BOXPLOT LOG2  ######################

my_Log2.TPM_Boxplot <- my_Log2_tpmf_avg %>%
  ggplot(aes(x = Week, y = mean_Log2_TPM)) +
  geom_boxplot(aes(fill = Outcome), width = 0.6, outlier.size = 0.9, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1, position = position_jitter(0.2)) +
  facet_grid(~ Outcome, scales = "free") +
  stat_compare_means(aes(group = Week), method = "t.test", paired = F, label = "p.format", size = 3, label.x = 1.4) +
  scale_fill_manual(values=my_fav_colors) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  labs(title = paste0(myGeneSet, ". Each dot is the average value of one gene"),
       subtitle = "Two-sample t-test. Mean is red diamond",
       x = "Oucome",
       y = "Log2(TPM+1)") +
  my_plot_themes + facet_themes
my_Log2.TPM_Boxplot
ggsave(my_Log2.TPM_Boxplot,
       file = paste0(myGeneSet, "_log2.pdf"),
       path = "Figures/GeneSet_Boxplots",
       width = 7, height = 5, units = "in")



###########################################################
###################  BOXPLOT LOG2 v2  #####################
# Order the facet differently

my_Log2.TPM_Boxplot_v2 <- my_Log2_tpmf_avg %>%
  ggplot(aes(x = Outcome, y = mean_Log2_TPM)) +
  geom_boxplot(aes(fill = Outcome), width = 0.6, outlier.size = 0.9, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1, position = position_jitter(0.2)) +
  facet_grid(~ Week, scales = "free") +
  stat_compare_means(aes(group = Outcome), method = "t.test", paired = F, label = "p.format", size = 3, label.x = 1.4) +
  scale_fill_manual(values=my_fav_colors) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  labs(title = paste0(myGeneSet, ". Each dot is the average value of one gene"),
       subtitle = "Paired Two-sample t-test. Mean is red diamond",
       x = "Oucome",
       y = "Log2(TPM+1)") +
  my_plot_themes + facet_themes
my_Log2.TPM_Boxplot_v2
ggsave(my_Log2.TPM_Boxplot_v2,
       file = paste0(myGeneSet, "_log2_v2.pdf"),
       path = "Figures/GeneSet_Boxplots",
       width = 7, height = 5, units = "in")

