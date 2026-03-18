# W0->W2 Using lm
# E. Lamont
# 3/18/26


# Seeing if there are any differences between the delta change from W0->W2 between cure and relapse samples

library(broom)
library(fgsea)

source("Import_data.R")
source("Import_GeneSets.R")


################################################
################# PREPARE DATA #################

source("Delta_Weeks0to2.R") # BothWk_log2TPMf_delta


BothWk_log2TPMf_delta$Outcome <- factor(BothWk_log2TPMf_delta$Outcome, 
                                        levels = c("Cure", "Relapse"))


################################################
################### BASIC LM ###################

lm_results <- BothWk_log2TPMf_delta %>%
  group_by(Gene) %>%
  do({
    fit <- lm(delta ~ Outcome, data = .)
    tidy(fit)
  }) %>%
  ungroup()

lm_results_clean <- lm_results %>%
  filter(term == "OutcomeRelapse") %>%
  dplyr::select(Gene, estimate, p.value)

# The estimate is ∆(Relapse - Cure)

lm_results_clean <- lm_results_clean %>%
  mutate(padj_fdr = p.adjust(p.value, method = "fdr"))


group_means <- BothWk_log2TPMf_delta %>%
  group_by(Gene, Outcome) %>%
  summarise(mean_delta = mean(delta), .groups = "drop") %>%
  pivot_wider(names_from = Outcome, values_from = mean_delta) %>%
  dplyr::rename("Cure_meanDelta" = "Cure",
         "Relapse_meanDelta" = "Relapse")

lm_results_clean <- lm_results_clean %>%
  left_join(group_means, by = "Gene")


ggplot(lm_results_clean, aes(x = estimate, y = -log10(padj_fdr))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Delta difference (Relapse - Cure)",
       y = "-log10(FDR)")

################################################
##################### LIMMA ####################

# reshape to matrix
# delta_mat <- BothWk_log2TPMf_delta %>%
#   dplyr::select(Gene, Patient, delta) %>%
#   pivot_wider(names_from = Patient, values_from = delta) %>%
#   column_to_rownames("Gene") %>%
#   as.matrix()
# 
# delta_meta <- BothWk_log2TPMf_delta %>%
#   distinct(Patient, Outcome)
# 
# design <- model.matrix(~ Outcome, data = delta_meta)
# 
# fit <- lmFit(delta_mat, design)
# fit <- eBayes(fit)
# # Warning message:
# # Zero sample variances detected, have been offset away from zero 
# 
# res <- topTable(fit, coef = "OutcomeRelapse", number = Inf)

# Looks pretty similar to the basic lm()

################################################
##################### GSEA #####################

# Rank genes
ranked_genes <- lm_results_clean %>%
  mutate(rank_stat = -log10(p.value) * sign(estimate)) %>% # This is using the non-adjusted pvalue! 
  dplyr::select(Gene, rank_stat)

# Convert to named vector
ranked_genes_vec <- ranked_genes$rank_stat
names(ranked_genes_vec) <- ranked_genes$Gene

# Sort decreasing
ranked_genes_vec <- sort(ranked_genes_vec, decreasing = TRUE)

# Load custom gene sets. Needs to have columns Gene, GeneSet. Will take from allGeneSetList
custom_list <- allGeneSetList[["EllaGeneSets_2026.01.30"]]

fgsea_res <- fgsea(
  pathways = custom_list,
  stats = ranked_genes_vec,
  minSize = 3
)

# Warning message:
#   In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                There are ties in the preranked stats (0.02% of the list).
#                              The order of those tied genes will be arbitrary, which may produce unexpected results.

# NES>0 is enriched in Relapse??
# NES<0 is enriched in Cure??


# Pick one pathway
plotEnrichment(custom_list[["dosR regulon"]], ranked_genes_vec)

fgsea_res %>%
  filter(pval < 0.05) %>%
  ggplot(aes(x = reorder(pathway, NES), y = NES)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

################################################
################## BUBBLE PLOT #################

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
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

# Need this for facet grouping
EllaGeneSets_2026.01.30 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2026.01.30.csv")

# Add some useful columns
fgsea_res2 <- fgsea_res %>%
  mutate(FDR_Significance = ifelse(padj < 0.05, "significant", "not significant"),
         FillColor = case_when(FDR_Significance == "significant" & NES > 0 ~ "pos",
                               FDR_Significance == "significant" & NES < 0 ~ "neg",
                               TRUE ~ "ns"),
         PathName_2 = paste0(pathway, " (n=", size, ")")) %>%
  left_join(EllaGeneSets_2026.01.30 %>% # Add the Group names
              dplyr::rename(pathway = GeneSet) %>%
              dplyr::select(pathway, Group) %>% 
              distinct(pathway, .keep_all = TRUE),
            by = "pathway") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))
  
my_bubblePlot <- fgsea_res2 %>% 
  ggplot(aes(x = NES, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = FillColor),
             size = 4, shape = 21, alpha = 0.9) + 
  scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
  facet_grid(rows = vars(Group_wrapped), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  # scale_x_continuous(limits = c(-4, 3.5), breaks = seq(-4, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "GSEA from linear model W0->W2 Relapse (n=5) vs Cure (n=8)", 
       subtitle = "Log2(W2.TPM+1)-Log2(W0.TPM+1) -> lm(delta ~ Outcome) -> GSEA \nFDR-adjusted p-values",
       y = NULL, x = "Normalized Enrichment Score (NES)") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
ggsave(my_bubblePlot,
       file = paste0("Delta_lm_v1", ".pdf"),
       path = "Figures/TimeCourse/Bubbles",
       width = 8, height = 8, units = "in")







