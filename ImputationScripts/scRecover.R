# scRecover for Imputation
# E. Lamont
# 3/25/26


# https://rdrr.io/bioc/scRecover/f/vignettes/scRecover.Rmd
# https://onlinelibrary.wiley.com/doi/10.1002/sim.10334

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRecover")

library(scRecover)
library(BiocParallel)
suppressMessages(library(SingleCellExperiment))


###########################################################
##################### ORGANIZE DATA #######################

head(W0_RawReadsf)
W0_Metadata



###########################################################
######################## scRecover ########################

scRecover(counts = W0_RawReadsf, labels = W0_Metadata$Week, outputDir = "ImputationScripts/outDir_scRecover/")

W0_scRecover_RawReadsf <- read.csv("ImputationScripts/outDir_scRecover/scRecover+scImpute.csv")

###########################################################
###################### CONVERT TO TPM #####################

source("Function_CalculateTPM.R")
W0_scRecover_tpmf <- CalculateTPM_RvOnly(W0_scRecover_RawReadsf)

###########################################################
############################ PCA ##########################

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
        plot.subtitle = element_text(size=9))

# Labelled Colors
my_fav_colors <- c(`Cure` = "#0072B2", `Relapse` = "#bc5300")
# Labelled Shapes
my_fav_shapes <- c(`Cure` = 21, `Relapse` = 21)

# Convert gene column to rownames
my_tpm <- W0_scRecover_tpmf # %>% column_to_rownames(var = "X")

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm)) # or my_tpm2

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 12.3% of variance
summary_PCA[2,1] # PC2 explains 5.1% of variance
summary_PCA[3,1] # PC3 explains 4.7% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, W0_Metadata, by = "SampleID", )

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Outcome)) + 
  geom_point(aes(fill = Outcome), shape = 21, size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors) +  
  # scale_shape_manual(values = my_fav_shapes) + 
  geom_text_repel(aes(label = Patient), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
  labs(title = "scRecover all Week 0",
       subtitle = "TPM filtered (Rv genes only)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig
# ** The samples that are farther away on the PCA are still the samples that have low txn coverage

