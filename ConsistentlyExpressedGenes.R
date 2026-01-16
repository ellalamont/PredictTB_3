# Looking for consistently expressed genes
# 1/15/26
# E Lamont

# After meeting with Evan Johnson and Howard Fan about using imputation to fill in the gaps of the partial genomes
# Want to know if there is a set of genes that is consistently expressed across the samples to know if imputation if worth it
# W0 and W2 separately
# Genes with >=100 reads 
# Want ~200 genes in 80% of samples


###########################################################
############################ W0 ###########################

# Keep only W0 sputum (ALL DATA)
W0SputumSamplesAll_pipeSummary <- my_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(Outcome != "Failure") %>%
  filter(N_Genomic > 1e6) %>% 
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome))

# Get a df of the raw reads
W0_RawReads <- All_RawReads_f %>%
  column_to_rownames("X") %>% 
  select(all_of(W0SputumSamplesAll_pipeSummary$SampleID2))

# What is 80% of the samples
min_samples <- ceiling(0.8 * ncol(W0_RawReads)) 
# 52 is 80% for all
# 38 is 80% for >1e6 reads

# Filter to get just the genes that have >=100 reads in at least 80% of the samples
genes_passing <- W0_RawReads %>%
  mutate(nSamples_w100reads = rowSums(across(everything(), ~ .x >= 100))) %>%
  filter(nSamples_w100reads >= min_samples)
nrow(genes_passing) 
# 155 genes for all
# 2,351 genes for >1e6 read samples

# Filter to get just the genes that have >=10 reads in at least 80% of the samples
# genes_passing <- W0_RawReads %>%
#   mutate(nSamples_w10reads = rowSums(across(everything(), ~ .x >= 10))) %>%
#   filter(nSamples_w10reads >= min_samples)
# nrow(genes_passing) # 904 genes 


###########################################################
############################ W2 ###########################

# Keep only W2 sputum (ALL DATA)
W2SputumSamplesAll_pipeSummary <- my_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>%
  filter(Outcome != "Failure") %>%
  filter(N_Genomic > 1e6) %>%
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome))

# Get a df of the raw reads
W2_RawReads <- All_RawReads_f %>%
  column_to_rownames("X") %>% 
  select(all_of(W2SputumSamplesAll_pipeSummary$SampleID2))

# What is 80% of the samples
min_samples <- ceiling(0.8 * ncol(W2_RawReads)) 
# 45 is 80% for all
# 24 is 80% for >1e6 reads

# Filter to get just the genes that have >=100 reads in at least 80% of the samples
genes_passing <- W2_RawReads %>%
  mutate(nSamples_w100reads = rowSums(across(everything(), ~ .x >= 100))) %>%
  filter(nSamples_w100reads >= min_samples)
nrow(genes_passing) 
# 15 genes for all
# 103 genes for >1e6 read samples

# Filter to get just the genes that have >=10 reads in at least 80% of the samples
# genes_passing <- W2_RawReads %>%
#   mutate(nSamples_w10reads = rowSums(across(everything(), ~ .x >= 10))) %>%
#   filter(nSamples_w10reads >= min_samples)
# nrow(genes_passing) # 227 genes




