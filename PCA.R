# PCA on all the biological samples passing filter
# E. Lamont
# 1/8/26

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R") # To get GoodSamples_tpmf and GoodSamples_pipeSummary 


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Labelled Colors
my_fav_colors <- c(`W0 sputum (cure)` = "#0072B2", `W0 sputum (relapse)` = "#bc5300", `Broth`= "black", `W2 sputum (cure)` = "#0072B2", `W2 sputum (relapse)` = "#bc5300", `THP1 spiked` = "#999999")
# Labelled Shapes
my_fav_shapes <- c(`W0 sputum (cure)` = 21, `W0 sputum (relapse)` = 21, `Broth`= 23, `W2 sputum (cure)` = 24, `W2 sputum (relapse)` = 24, `THP1 spiked` = 23)

###########################################################
############### PCA GOODSAMPLES TPM_F 60% #################
# >60% genes with at least 10 reads, already subsetted in Import_data.R
# Rv GENES ONLY INCLUDED IN TXN COVERAGE

# Convert gene column to rownames
my_tpm <- GoodSamples60_tpmf %>% # column_to_rownames(var = "X")
  select(!contains("THP1")) # Remove the spiked samples

# Keep genes with >5 tpm (should be using counts!) in at least 50% of samples 
# 1/16/26: Not doing this now but maybe should?? Doesn't really change what the graph looks like
# keep <- rowSums(my_tpm > 5) >= 0.5 * ncol(my_tpm)
# my_tpm2 <- my_tpm[keep, ] # now only 3948 genes (instead of 4030)

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm)) # or my_tpm2

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 14.6% of variance
summary_PCA[2,1] # PC2 explains 8% of variance
summary_PCA[3,1] # PC3 explains 6.9% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID2 = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, GoodSamples60_pipeSummary, by = "SampleID2", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type2, shape = Type2)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = Patient), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >60% genes with at least 10 reads (Run1-4)",
       subtitle = "TPM filtered (Rv genes only)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
# ggsave(PCA_fig,
#        file = paste0("GoodSamples_tpmf_txnCov60_v1.pdf"),
#        path = "Figures/PCA",
#        width = 10, height = 6, units = "in")


###########################################################
############ PCA GOODSAMPLES TPM_F 60% - W0 ONLY ##########
# >60% genes with at least 10 reads, already subsetted in Import_data.R
# Rv GENES ONLY INCLUDED IN TXN COVERAGE

# Convert gene column to rownames
my_tpm <- GoodSamples60_tpmf %>% # column_to_rownames(var = "X")
  select(contains("W0")) # Only keep W0 samples

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 9.9% of variance
summary_PCA[2,1] # PC2 explains 9.1% of variance
summary_PCA[3,1] # PC3 explains 7.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID2 = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, GoodSamples60_pipeSummary, by = "SampleID2", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type2, shape = Type2)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  geom_text_repel(aes(label = SampleID2), size= 2, box.padding = 0.4, segment.color = "grey", max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >60% genes with at least 10 reads (Run1-4) W0 only",
       subtitle = "TPM filtered (Rv genes only) Numbers are % txn coverage",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
# ggsave(PCA_fig,
#        file = paste0("W0_GoodSamples_tpmf_txnCov60_v2.pdf"),
#        path = "Figures/PCA",
#        width = 10, height = 6, units = "in")



###########################################################
############ PCA GOODSAMPLES TPM_F 60% - W2 ONLY ##########
# >60% genes with at least 10 reads, already subsetted in Import_data.R
# Rv GENES ONLY INCLUDED IN TXN COVERAGE

# Convert gene column to rownames
my_tpm <- GoodSamples60_tpmf %>% # column_to_rownames(var = "X")
  select(contains("W2")) # Only keep W0 samples

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 15.1% of variance
summary_PCA[2,1] # PC2 explains 13.6% of variance
summary_PCA[3,1] # PC3 explains 11.3% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID2 = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, GoodSamples60_pipeSummary, by = "SampleID2", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type2, shape = Type2)) + 
  geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = N_Genomic), size= 2, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >60% genes with at least 10 reads (Run1-4) W2 only",
       subtitle = "TPM filtered (Rv genes only)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
# ggsave(PCA_fig,
#        file = paste0("W2_GoodSamples_tpmf_txnCov60_v1.pdf"),
#        path = "Figures/PCA",
#        width = 10, height = 6, units = "in")


