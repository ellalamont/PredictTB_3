# Unsupervised Clustering
# 1/15/26


# https://www.biostars.org/p/9495904/
# https://ngs101.com/how-to-cluster-rna-seq-data-to-uncover-gene-expression-patterns-hierarchical-and-k-means-methods-for-absolute-beginners/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8386505/

# Hierarchical clustering?

###########################################################
##################### ORGANIZE DATA #######################

# Keep only W0 sputum
W0SputumSamples60_pipeSummary <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(Outcome != "Failure")

P_metadata <- W0SputumSamples60_pipeSummary %>% select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)

P_TPM <- GoodSamples60_tpmf %>% select(all_of(W0SputumSamples60_pipeSummary$SampleID2))

# Keep genes with >5 counts in at least 50% of samples 
# (https://ngs101.com/how-to-cluster-rna-seq-data-to-uncover-gene-expression-patterns-hierarchical-and-k-means-methods-for-absolute-beginners/)
keep <- rowSums(P_TPM > 5) >= 0.5 * ncol(P_TPM)
P_TPM_2 <- P_TPM[keep, ] # now only 3942 genes (instead of 4030)

# Put everything in one dataframe
P_TPM_t <- P_TPM_2 %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")
# my_df2 <- my_df %>% select(-SampleID2)

# Remove na
# my_df2 <- na.omit(my_df2) # Didn't change anything

# Outcome must be factor with "Relapse" first
# my_df2$Outcome <- factor(my_df2$Outcome, levels = c("Relapse", "Cure"))

###########################################################
######################### HCLUST ##########################

# Using TPM... but should I be using DEG?








