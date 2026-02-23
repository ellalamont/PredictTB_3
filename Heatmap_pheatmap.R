# Try a heatmap with pheatmap
# E. Lamont 
# 10/29/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Import_data.R")

# Just want the good sputum and broth


###########################################################
###################### PROCESS DATA #######################

my_tpm <- GoodSamples60_tpmf %>% select(-contains("THP1"))

###########################################################
#################### TESTING PHEATMAP #####################

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

# Define the colors
my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "W0 sputum (failure)" = "orange2",
            "W2 sputum (cure)" = "green4",
            "W2 sputum (relapse)" = "#6A3D9A", 
            "W2 sputum (failure)" = "orange4",
            "Broth" = "#999999")
)


pheatmap(testing, 
         annotation_col = pheatmap_pipeSummary, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_cols = 5)



######## NOTHING BELOW HAS BEEN DONE #########

# # Try to get good row annotations based on MTb functional group
# # Start with MTb.TB.Phenotypes.AllGeneSets
# # Convert to dataframe
# Gene_Category <- do.call(rbind, lapply(names(Walter2015GeneSets), function(category) {
#   data.frame(Gene = Walter2015GeneSets[[category]], Category = category, stringsAsFactors = FALSE)
# })) %>% 
#   filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
#   distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
#   column_to_rownames(var = "Gene")
# 
# 
# pheatmap(my_tpm_2_matrix, 
#          # annotation_row = Gene_Category, 
#          fontsize_row = 1,
#          annotation_col = my_pipeSummary["Week"],
#          annotation_colors = my_annotation_colors,
#          scale = "row",
#          cutree_rows = 8,
#          cutree_cols = 5)
# 
# 
# 
# 
# 
# # Pull out what the clusters are: 
# # use silent = TRUE to suppress the plot
# my_heatmap <- pheatmap(my_tpm_2_matrix, 
#                        annotation_row = Gene_Category, 
#                        annotation_col = my_pipeSummary["Week"],
#                        annotation_colors = my_annotation_colors,
#                        scale = "row",
#                        cutree_rows = 7,
#                        cutree_cols = 5,
#                        silent = TRUE)
# # Extract the row clustering information
# row_clusters <- cutree(my_heatmap$tree_row, k = 7)
# # Convert to a data frame for easier handling
# row_cluster_df <- data.frame(Gene = names(row_clusters), Cluster = row_clusters)
# # View the first few rows
# head(row_cluster_df)
# 
# 
# ###########################################################
# ############ W0 AND BROTH WITH CLUSTERING #################
# 
# my_tpm_3_matrix <- my_tpm_2 %>% select(-c("S_503937", "S_577208", "S_575533_MtbrRNA")) %>%
#   as.matrix()
# 
# Gene_Category <- do.call(rbind, lapply(names(MTb.TB.Phenotypes.AllGeneSets), function(category) {
#   data.frame(Gene = MTb.TB.Phenotypes.AllGeneSets[[category]], Category = category, stringsAsFactors = FALSE)
# })) %>% 
#   # filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
#   distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
#   column_to_rownames(var = "Gene")
# 
# pheatmap(my_tpm_3_matrix, 
#          annotation_row = Gene_Category, 
#          fontsize_row = 1,
#          annotation_col = my_pipeSummary["Week"],
#          annotation_colors = my_annotation_colors,
#          scale = "row",
#          cutree_rows = 6,
#          cutree_cols = 2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
