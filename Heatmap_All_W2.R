# Try a heatmap with pheatmap
# E. Lamont 
# 10/29/25

# Comparing TPM, log2(TPM+1) and top variable genes with vars()

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Import_data.R")

# Just want the good sputum and broth

# Define the colors
# my_annotation_colors <- list(
#   Type2 = c("W0 sputum (cure)" = "#0072B2",
#             "W0 sputum (relapse)" = "red", 
#             "W0 sputum (failure)" = "orange2",
#             "W2 sputum (cure)" = "green4",
#             "W2 sputum (relapse)" = "#6A3D9A", 
#             "W2 sputum (failure)" = "orange4",
#             "Broth" = "#999999"))

my_annotation_colors <- list(
  Type2 = c("W2 sputum (cure)" = "green4",
            "W2 sputum (relapse)" = "#6A3D9A"))


###########################################################
#################### TESTING PHEATMAP #####################

my_tpm <- GoodSamples60_tpmf %>% select(-contains("THP1"))


# pheatmap(All_tpm_matrix[1:10,], scale = "row")

testing <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]]) # Guess this doesn't need to be a matrix

# pheatmap(testing, scale = "row")

# Grab the columns needed to give colors
pheatmap_pipeSummary <- GoodSamples60_pipeSummary %>%
  filter(SampleID2 %in% colnames(testing)) %>%
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
pheatmap_pipeSummary <- pheatmap_pipeSummary[colnames(testing), , drop = FALSE]

pheatmap(testing, 
         annotation_col = pheatmap_pipeSummary, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 5)


###########################################################
#################### ALL W2 TPM TOP200 ####################
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

### Subset just the top 200 expressed genes ##
top200 <- Pheatmap_tpm %>%
  rowMeans() %>% # Find the mean value of each gene
  sort(decreasing = TRUE) %>%
  names() %>%
  head(200)

Pheatmap_tpm_Top200 <- Pheatmap_tpm[top200,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top200)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top200), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top200)))

pheatmap(Pheatmap_tpm_Top200, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)

###########################################################
################ ALL W2 Log2(TPM+1) TOP200 ################
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

### Subset just the top 200 expressed genes ##
top200 <- Pheatmap_log2TPM %>%
  rowMeans() %>% # Find the mean value of each gene
  sort(decreasing = TRUE) %>%
  names() %>%
  head(200)

Pheatmap_tpm_Top200 <- Pheatmap_log2TPM[top200,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top200)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top200), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top200)))

pheatmap(Pheatmap_tpm_Top200, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)

###########################################################
############ ALL W2 Log2(TPM+1) VARIABLE TOP200 ###########
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

# See the most variable genes (across all samples)
vars <- apply(Pheatmap_log2TPM, 1, var)
top200 <- names(sort(vars, decreasing = TRUE))[1:200]

Pheatmap_tpm_Top200 <- Pheatmap_log2TPM[top200,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top200)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top200), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top200)))

pheatmap(Pheatmap_tpm_Top200, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)


###########################################################
#################### ALL W2 TPM TOP50 #####################
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

### Subset just the top 200 expressed genes ##
top50 <- Pheatmap_tpm %>%
  rowMeans() %>% # Find the mean value of each gene
  sort(decreasing = TRUE) %>%
  names() %>%
  head(50)

Pheatmap_tpm_Top50 <- Pheatmap_tpm[top50,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top50)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top50), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top50)))

pheatmap(Pheatmap_tpm_Top50, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)


###########################################################
################ ALL W2 log2(TPM+1) TOP50 #################
# 4/7/26

# Grab the TPM
All_tpmf_log2 <- log2(All_tpm_f + 1)
Pheatmap_tpm <- All_tpmf_log2 %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

### Subset just the top 200 expressed genes ##
top50 <- Pheatmap_tpm %>%
  rowMeans() %>% # Find the mean value of each gene
  sort(decreasing = TRUE) %>%
  names() %>%
  head(50)

Pheatmap_tpm_Top50 <- Pheatmap_tpm[top50,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top50)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top50), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top50)))

pheatmap(Pheatmap_tpm_Top50, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)

###########################################################
############ ALL W2 Log2(TPM+1) VARIABLE TOP50 ############
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

# See the most variable genes (across all samples)
vars <- apply(Pheatmap_log2TPM, 1, var)
top50 <- names(sort(vars, decreasing = TRUE))[1:50]

Pheatmap_tpm_Top50 <- Pheatmap_log2TPM[top50,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top50)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top50), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top50)))

pheatmap(Pheatmap_tpm_Top50, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)


###########################################################
################ ALL W2 log2(TPM+1) TOP20 #################
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

### Subset just the top 20 expressed genes ##
top20 <- Pheatmap_log2TPM %>%
  rowMeans() %>% # Find the mean value of each gene
  sort(decreasing = TRUE) %>%
  names() %>%
  head(20)

Pheatmap_tpm_Top20 <- Pheatmap_tpm[top20,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top20)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top20), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top20)))

pheatmap(Pheatmap_tpm_Top20, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 1)


###########################################################
############ ALL W2 Log2(TPM+1) VARIABLE TOP20 ############
# 4/7/26
# BEST ONE

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 10) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

# See the most variable genes (across all samples)
vars <- apply(Pheatmap_log2TPM, 1, var)
top20 <- names(sort(vars, decreasing = TRUE))[1:20]

Pheatmap_tpm_Top20 <- Pheatmap_log2TPM[top20,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top20)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top20), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top20)))

pheatmap(Pheatmap_tpm_Top20, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 2,
         main = "All W2: Keep TPM>10 in >70% -> Log2(TPM+1) -> vars() -> top20Genes")

###########################################################
############ ALL W2 Log2(TPM+1) VARIABLE TOP30 ############
# 4/7/26

# Grab the TPM
Pheatmap_tpm <- All_tpm_f %>% 
  select(any_of(my_pipeSummary$SampleID2)) %>% # Only keep samples that are in my_pipeSummary
  select(contains("W2")) %>% # There are only 54 samples here.. should there be 56?
  filter(rowSums(.) != 0) # Remove genes that are 0

# Keep only genes that are TPM >= 30 in >70% of samples
keep <- rowMeans(Pheatmap_tpm >= 30) >= 0.7
Pheatmap_tpm <- Pheatmap_tpm[keep, ] # Now only 842 genes

# Log2 transform
Pheatmap_log2TPM <- log2(Pheatmap_tpm + 1)

# See the most variable genes (across all samples)
vars <- apply(Pheatmap_log2TPM, 1, var)
top30 <- names(sort(vars, decreasing = TRUE))[1:30]

Pheatmap_tpm_Top30 <- Pheatmap_log2TPM[top30,]

# Grab the metadata
Pheatmap_Metadata <- my_pipeSummary %>% 
  filter(SampleID2 %in% colnames(Pheatmap_tpm_Top30)) %>% 
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Pheatmap_Metadata <- Pheatmap_Metadata[colnames(Pheatmap_tpm_Top30), , drop = FALSE]
stopifnot(all(Pheatmap_Metadata$SampleID2 == colnames(Pheatmap_tpm_Top30)))

pheatmap(Pheatmap_tpm_Top30, 
         annotation_col = Pheatmap_Metadata, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 2)
