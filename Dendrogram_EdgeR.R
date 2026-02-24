# Clustering with EdgeR
# E. Lamont
# 2/23/26

# Not sure the best way to cluster but let's try EdgeR methods
## This starts with raw reads and makes CPM

# https://www.r-bloggers.com/2022/08/downstream-bioinformatics-analysis-of-omics-data-with-edger/
# https://stats.stackexchange.com/questions/459063/best-practices-in-the-selection-of-distance-metric-and-clustering-methods-for-ge


source("Import_data.R")

library(edgeR)
# library(dendextend)
library(ggtree)
library(ape)
library(ggtreeExtra)

###########################################################
###################### PROCESS DATA #######################

# Clean up the raw reads and make a matrix
matrix_RawReads <- GoodSamples60_RawReadsf %>% 
  select(-contains("THP1")) %>%
  column_to_rownames("X") %>%
  as.matrix()

# Convert to an EdgeR dge format
dge <- DGEList(counts = matrix_RawReads)

# Normalize the data
# method="TMM" is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data. https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors
dge2 <- calcNormFactors(dge, method = "TMM")

# Convert to CPM
logCPM <- cpm(dge2, 
              log = T, prior.count = 1, # Log2 transformation (+1)
              normalized.lib.sizes = T) # Because we normalized with TMM

# Make a distance matrix (pairwise distances between each sample)
Euc_dist_matrix <- dist(t(logCPM)) # This is Euclidean distances, looks at differences
# Ignore the error....
Corr_dist_matrix <- as.dist(1-cor(logCPM)) # This is correlation distances which looks more at similarity

# look at the different types of distance calculations
par(mfrow=c(1,2))
plot(Euc_dist_matrix, main="Euclidean")
plot(Corr_dist_matrix, main="Correlation")


###########################################################
#################### EUCLIDEAN WARD.D2 ####################

# Generate the clusters
hc_euc <- hclust(Euc_dist_matrix, method = "ward.D2")
# plot(hc_euc)

# Convert hclust to phylo object
phylo_tree <- ape::as.phylo(hc_euc)

# Create clean metadata
tree_metadata <- data.frame(SampleID2 = phylo_tree$tip.label)

# Match the pipeSummary info
tree_metadata$Type2 <- my_pipeSummary$Type2[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Outcome <- my_pipeSummary$Outcome[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Week <- my_pipeSummary$Week[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]

# Determine colors
tree_metadata[is.na(tree_metadata)] <- "Broth"
Week_Colors <- c("Week 0" = "#0072B2",
                 "Week 2" = "#bc5300",
                 "Broth" = "#999999")
Outcome_Colors <- c("Cure" = "aquamarine2",
                    "Relapse" = "red3",
                    "Broth" = "#999999")

# Ensure metadata is in tip order
tree_metadata <- tree_metadata %>%
  filter(SampleID2 %in% phylo_tree$tip.label) %>%
  slice(match(phylo_tree$tip.label, SampleID2))

# Compute distance from root for all tips
tip_depths <- node.depth.edgelength(phylo_tree)[1:length(phylo_tree$tip.label)]
max_tip <- max(tip_depths)

tree1 <- ggtree(phylo_tree, layout = "rectangular") + 
  geom_tiplab(size = 2.5, align = T, offset = 10) + 
  geom_fruit(data = tree_metadata, geom = geom_tile, 
             mapping = aes(y = SampleID2, fill = Week), 
             width = 5, offset = 0) + 
  scale_fill_manual(name = "Week", values = Week_Colors) + 
  ggnewscale::new_scale_fill() + # Start a new fill scale for Outcome so they don't overwrite
  geom_fruit(data = tree_metadata, geom = geom_tile,
             mapping = aes(y = SampleID2, fill = Outcome), 
             width = 5, offset = 0.02) + 
  scale_fill_manual(name = "Outcome", values = Outcome_Colors) + 
  xlim(0, max_tip*1.2) +   # 20% extra space
  labs(title = "Euclidean distance, ward.D2 method")
tree1
ggsave(tree1,
       file = paste0("Euclidean_WardD2.pdf"),
       path = "Figures/Dendrograms",
       width = 10, height = 6, units = "in")


###########################################################
################### CORRELATION WARD.D2 ###################

# Generate the clusters
hc_corr <- hclust(Corr_dist_matrix, method = "ward.D2")
# plot(hc_corr)

# Convert hclust to phylo object
phylo_tree <- ape::as.phylo(hc_corr)

# Create clean metadata
tree_metadata <- data.frame(SampleID2 = phylo_tree$tip.label)

# Match the pipeSummary info
tree_metadata$Type2 <- my_pipeSummary$Type2[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Outcome <- my_pipeSummary$Outcome[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Week <- my_pipeSummary$Week[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]

# Determine colors
tree_metadata[is.na(tree_metadata)] <- "Broth"
Week_Colors <- c("Week 0" = "#0072B2",
                 "Week 2" = "#bc5300",
                 "Broth" = "#999999")
Outcome_Colors <- c("Cure" = "aquamarine2",
                    "Relapse" = "red3",
                    "Broth" = "#999999")

# Ensure metadata is in tip order
tree_metadata <- tree_metadata %>%
  filter(SampleID2 %in% phylo_tree$tip.label) %>%
  slice(match(phylo_tree$tip.label, SampleID2))

# Compute distance from root for all tips
tip_depths <- node.depth.edgelength(phylo_tree)[1:length(phylo_tree$tip.label)]
max_tip <- max(tip_depths)

tree1 <- ggtree(phylo_tree, layout = "rectangular") + 
  geom_tiplab(size = 2.5, align = T, offset = 0.03) + 
  geom_fruit(data = tree_metadata, geom = geom_tile, 
             mapping = aes(y = SampleID2, fill = Week), 
             width = 0.02, offset = 0) + 
  scale_fill_manual(name = "Week", values = Week_Colors) + 
  ggnewscale::new_scale_fill() + # Start a new fill scale for Outcome so they don't overwrite
  geom_fruit(data = tree_metadata, geom = geom_tile,
             mapping = aes(y = SampleID2, fill = Outcome),
             width = 0.02, offset = 0.02) +
  scale_fill_manual(name = "Outcome", values = Outcome_Colors) +
  xlim(0, max_tip*1.2) +   # 20% extra space
  labs(title = "Correlation distance, ward.D2 method")
tree1
ggsave(tree1,
       file = paste0("Correlation_WardD2.pdf"),
       path = "Figures/Dendrograms",
       width = 10, height = 6, units = "in")

###########################################################
################### EUCLIDEAN COMPLETE ####################

# Generate the cluster
hc_euc2 <- hclust(Euc_dist_matrix, method = "complete")
# plot(hc_euc2)

# Convert hclust to phylo object
phylo_tree <- ape::as.phylo(hc_euc2)

# Create clean metadata
tree_metadata <- data.frame(SampleID2 = phylo_tree$tip.label)

# Match the pipeSummary info
tree_metadata$Type2 <- my_pipeSummary$Type2[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Outcome <- my_pipeSummary$Outcome[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Week <- my_pipeSummary$Week[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]

# Determine colors
tree_metadata[is.na(tree_metadata)] <- "Broth"
Week_Colors <- c("Week 0" = "#0072B2",
                 "Week 2" = "#bc5300",
                 "Broth" = "#999999")
Outcome_Colors <- c("Cure" = "aquamarine2",
                    "Relapse" = "red3",
                    "Broth" = "#999999")

# Ensure metadata is in tip order
tree_metadata <- tree_metadata %>%
  filter(SampleID2 %in% phylo_tree$tip.label) %>%
  slice(match(phylo_tree$tip.label, SampleID2))

# Compute distance from root for all tips
tip_depths <- node.depth.edgelength(phylo_tree)[1:length(phylo_tree$tip.label)]
max_tip <- max(tip_depths)

tree1 <- ggtree(phylo_tree, layout = "rectangular") + 
  geom_tiplab(size = 2.5, align = T, offset = 10) + 
  geom_fruit(data = tree_metadata, geom = geom_tile, 
             mapping = aes(y = SampleID2, fill = Week), 
             width = 5, offset = 0) + 
  scale_fill_manual(name = "Week", values = Week_Colors) + 
  ggnewscale::new_scale_fill() + # Start a new fill scale for Outcome so they don't overwrite
  geom_fruit(data = tree_metadata, geom = geom_tile,
             mapping = aes(y = SampleID2, fill = Outcome), 
             width = 5, offset = 0.02) + 
  scale_fill_manual(name = "Outcome", values = Outcome_Colors) + 
  xlim(0, max_tip*1.2) +   # 20% extra space
  labs(title = "Euclidean distance, complete method")
tree1
ggsave(tree1,
       file = paste0("Euclidean_complete.pdf"),
       path = "Figures/Dendrograms",
       width = 10, height = 6, units = "in")

###########################################################
################## CORRELATION COMPLETE ###################

# Generate the clusters
hc_corr2 <- hclust(Corr_dist_matrix, method = "complete")
# plot(hc_corr2)

# Convert hclust to phylo object
phylo_tree <- ape::as.phylo(hc_corr2)

# Create clean metadata
tree_metadata <- data.frame(SampleID2 = phylo_tree$tip.label)

# Match the pipeSummary info
tree_metadata$Type2 <- my_pipeSummary$Type2[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Outcome <- my_pipeSummary$Outcome[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]
tree_metadata$Week <- my_pipeSummary$Week[match(tree_metadata$SampleID2, my_pipeSummary$SampleID2)]

# Determine colors
tree_metadata[is.na(tree_metadata)] <- "Broth"
Week_Colors <- c("Week 0" = "#0072B2",
                 "Week 2" = "#bc5300",
                 "Broth" = "#999999")
Outcome_Colors <- c("Cure" = "aquamarine2",
                    "Relapse" = "red3",
                    "Broth" = "#999999")

# Ensure metadata is in tip order
tree_metadata <- tree_metadata %>%
  filter(SampleID2 %in% phylo_tree$tip.label) %>%
  slice(match(phylo_tree$tip.label, SampleID2))

# Compute distance from root for all tips
tip_depths <- node.depth.edgelength(phylo_tree)[1:length(phylo_tree$tip.label)]
max_tip <- max(tip_depths)

tree1 <- ggtree(phylo_tree, layout = "rectangular") + 
  geom_tiplab(size = 2.5, align = T, offset = 0.015) + 
  geom_fruit(data = tree_metadata, geom = geom_tile, 
             mapping = aes(y = SampleID2, fill = Week), 
             width = 0.02, offset = 0) + 
  scale_fill_manual(name = "Week", values = Week_Colors) + 
  ggnewscale::new_scale_fill() + # Start a new fill scale for Outcome so they don't overwrite
  geom_fruit(data = tree_metadata, geom = geom_tile,
             mapping = aes(y = SampleID2, fill = Outcome),
             width = 0.01, offset = 0.02) +
  scale_fill_manual(name = "Outcome", values = Outcome_Colors) +
  xlim(0, max_tip*1.2) +   # 20% extra space
  labs(title = "Correlation distance, complete method")
tree1
ggsave(tree1,
       file = paste0("Correlation_complete.pdf"),
       path = "Figures/Dendrograms",
       width = 10, height = 6, units = "in")



