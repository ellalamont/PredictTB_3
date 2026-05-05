# Compare genes DEG for each W2 relapse vs all cures
# E. Lamont
# 5/5/26

source("Import_DEG_sets.R")

# https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html
# install.packages("ggVennDiagram")
library(ggVennDiagram)


# Plot basics
my_plot_themes <- theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9))

###########################################################
#################### PROCESS THE DATA #####################

# These are all with p<0.05 and log2fold>1

genelist_DE1_UP <- list(
  "W2_11031" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_11031.vs.W2_cures"]] %>%
    filter(LOG2FOLD >= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12020" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12020.vs.W2_cures"]] %>%
    filter(LOG2FOLD >= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12073" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12073.vs.W2_cures"]] %>%
    filter(LOG2FOLD >= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12084" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12084.vs.W2_cures"]] %>%
    filter(LOG2FOLD >= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_13041" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_13041.vs.W2_cures"]] %>%
    filter(LOG2FOLD >= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID))

genelist_DE1_DOWN <- list(
  "W2_11031" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_11031.vs.W2_cures"]] %>%
    filter(LOG2FOLD <= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12020" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12020.vs.W2_cures"]] %>%
    filter(LOG2FOLD <= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12073" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12073.vs.W2_cures"]] %>%
    filter(LOG2FOLD <= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_12084" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_12084.vs.W2_cures"]] %>%
    filter(LOG2FOLD <= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID),
  "W2_13041" = DEG.list_W2Relapses.vs.AllCures_f2[["MTb.W2_13041.vs.W2_cures"]] %>%
    filter(LOG2FOLD <= 1 & AVG_PVALUE < 0.05) %>% pull(GENE_ID))


###########################################################
################### SAVE LISTS AS CSV #####################

# UP list
max_len <- max(lengths(genelist_DE1_UP))
genelist_padded <- lapply(genelist_DE1_UP, function(x) {
  length(x) <- max_len
  x
})
df_out <- bind_cols(genelist_padded)
colnames(df_out) <- names(genelist_DE1_UP)
write.csv(df_out, "Data/VennDiagram_Genes/W2_EachRelapse.Vs.AllCures_DE1_UP.csv", row.names = FALSE)

# DOWN list
max_len <- max(lengths(genelist_DE1_DOWN))
genelist_padded <- lapply(genelist_DE1_DOWN, function(x) {
  length(x) <- max_len
  x
})
df_out <- bind_cols(genelist_padded)
colnames(df_out) <- names(genelist_DE1_DOWN)
write.csv(df_out, "Data/VennDiagram_Genes/W2_EachRelapse.Vs.AllCures_DE1_DOWN.csv", row.names = FALSE)



###########################################################
########### VENN DIAGRAM - UP REGULATED GENES #############

Venn_DEG_UP <- ggVennDiagram(genelist_DE1_UP,
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 650),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes each W2 relapse compared to cure(n=8)",
       subtitle = "P<0.05, Log2fold>1",
       fill = "Number of genes") 
Venn_DEG_UP
# ggsave(Venn_DEG_UP,
#        file = "W2.EachRelapse.Vs.AllCure_DE1.UP.pdf",
#        path = "Figures/VennDiagrams",
#        width = 9, height = 6, units = "in")

# Make it plotly! 
ggVennDiagram(genelist_DE1_UP,
              show_intersect = T)


###########################################################
########### VENN DIAGRAM - DOWN REGULATED GENES #############

Venn_DEG_DOWN <- ggVennDiagram(genelist_DE1_DOWN,
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 350),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes each W2 relapse compared to cure(n=8)",
       subtitle = "P<0.05, Log2fold>1",
       fill = "Number of genes") 
Venn_DEG_DOWN
# ggsave(Venn_DEG_DOWN,
#        file = "W2.EachRelapse.Vs.AllCure_DE1.DOWN.pdf",
#        path = "Figures/VennDiagrams",
#        width = 9, height = 6, units = "in")

# Make it plotly! 
ggVennDiagram(genelist_DE1_DOWN,
              show_intersect = T)

