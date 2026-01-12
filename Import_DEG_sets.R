# Import Bob's differential gene expression data and the metagenesets data
# 1/12/26
# E. Lamont

source("Import_data.R")

# Data is coming from the Lenovo PredictTB_allRuns
# 1/12/26: Using the 60% Txn coverage cutoff runs 1-4


###########################################################
################### IMPORT BOB's DE DATA ##################

`W0.cure.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_Ra/W0_cure.MTb.Meta.JOINED.txt")
`W0.relapse.ComparedTo.W0.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W0.relapse/W0_relapse.MTb.Meta.JOINED.txt")
`W2.cure.ComparedTo.W0.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W2.cure/W2_cure.MTb.Meta.JOINED.txt")
`W0.relapse.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.relapse_vs_Ra/W0_relapse.MTb.Meta.JOINED.txt")
`W2.cure.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W2.cure_vs_Ra/W2_cure.MTb.Meta.JOINED.txt")
`W2.relapse.ComparedTo.W2.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W2.cure_vs_W2.relapse/W2_relapse.MTb.Meta.JOINED.txt")

# Haven't done this for this R project
# Comparing the same sample run twice
# W0_13051_Compare <- read.delim("Data/Differential_Expression_Run1to3/W0_13051/W0_13051_S42.MTb.Meta.JOINED.txt")
# W2_11058_Compare <- read.delim("Data/Differential_Expression_Run1to3/W2_11058/W2_11058_S31.MTb.Meta.JOINED.txt")
# W2_12010_Compare <- read.delim("Data/Differential_Expression_Run1to3/W2_12010/W2_12010_S33.MTb.Meta.JOINED.txt")
# W2_12043_Compare <- read.delim("Data/Differential_Expression_Run1to3/W2_12043/W2_12043_S25.MTb.Meta.JOINED.txt")
# W2_13026_Compare <- read.delim("Data/Differential_Expression_Run1to3/W2_13026/W2_13026_S25.MTb.Meta.JOINED.txt")
# THP1_1e6_Compare <- read.delim("Data/Differential_Expression_Run1to3/THP1_1e6/THP1_1e6_1_S67.MTb.Meta.JOINED.txt")


###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`W0.cure.ComparedTo.Ra`,
                 `W0.relapse.ComparedTo.W0.cure`, 
                 `W2.cure.ComparedTo.W0.cure`, 
                 `W0.relapse.ComparedTo.Ra`,
                 `W2.cure.ComparedTo.Ra`,
                 `W2.relapse.ComparedTo.W2.cure`)

# Make a list of the names
df_names <- c("W0.cure.ComparedTo.Ra",
              "W0.relapse.ComparedTo.W0.cure", 
              "W2.cure.ComparedTo.W0.cure", 
              "W0.relapse.ComparedTo.Ra",
              "W2.cure.ComparedTo.Ra",
              "W2.relapse.ComparedTo.W2.cure")

# Give the df list the correct df names
names(list_dfs) <- df_names


###########################################################
################# REMOVE NON-CODING GENES #################
# 8/15/25: After talking to DRS, decided to remove all the non-coding genes and all the MT genes, and leave just the coding genes starting with Rv. So need to remove these at the raw read level and calculate new TPM
# The Pathcap people also had issues with ncRNAs: https://www.nature.com/articles/s41598-019-55633-6#Sec8

list_dfs_f <- lapply(list_dfs, function(df) {
  df %>% filter(str_detect(GENE_ID, "^Rv\\d+.*"))
})


###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# 10/10/25: Adjusting this so it is all FDR p-values and I have Log2Fold >1 and >2

# Make a new list to hold dataframes with extra columns
list_dfs_f2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs_f)) {
  
  current_df <- list_dfs_f[[i]]
  current_df_name <- df_names[i]
  
  # Calculate the FDR p-value
  current_df$FDR_PVALUE <- p.adjust(current_df$AVG_PVALUE, method = "fdr")
  
  # Columns for Log2Fold > 1
  current_df$DE1 <- ifelse(current_df$LOG2FOLD < -1 & current_df$FDR_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 1 & current_df$FDR_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE1 <- factor(current_df$DE1, levels = ordered_DE)
  current_df$DE1_labels <- ifelse(current_df$DE1 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 2
  current_df$DE2 <- ifelse(current_df$LOG2FOLD < -2 & current_df$FDR_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 2 & current_df$FDR_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2 <- factor(current_df$DE2, levels = ordered_DE)
  current_df$DE2_labels <- ifelse(current_df$DE2 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 1 and AVG_PVALUE < 0.05
  current_df$DE1_ogP <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE1_ogP <- factor(current_df$DE1_ogP, levels = ordered_DE)
  current_df$DE1_ogP_labels <- ifelse(current_df$DE1_ogP != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 2 and AVG_PVALUE < 0.05
  current_df$DE2_ogP <- ifelse(current_df$LOG2FOLD < -2 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 2 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2_ogP <- factor(current_df$DE2_ogP, levels = ordered_DE)
  current_df$DE2_ogP_labels <- ifelse(current_df$DE2_ogP != "not significant", current_df$GENE_NAME, NA)

  list_dfs_f2[[current_df_name]] <- current_df
}




###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############

`MetaGeneSets_W0.cure.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_Ra/W0_cure.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W0.relapse.ComparedTo.W0.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W0.relapse/W0_relapse.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2.cure.ComparedTo.W0.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W2.cure/W2_cure.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W0.relapse.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W0.relapse_vs_Ra/W0_relapse.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2.cure.ComparedTo.Ra` <- read.delim("Data/DE_Run1to4_60TxnCov/W2.cure_vs_Ra/W2_cure.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2.relapse.ComparedTo.W2.cure` <- read.delim("Data/DE_Run1to4_60TxnCov/W2.cure_vs_W2.relapse/W2_relapse.MTb.MetaGeneSets.UP.txt")

# Specific gene sets I made
TAR_Poonawala2024_W2.cure.ComparedTo.W0.cure <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W2.cure_TAR/W2_cure.MTb.MetaGeneSets.UP.txt")
EllaGeneSets_2025.11.05_W2.cure.ComparedTo.W0.cure <- read.delim("Data/DE_Run1to4_60TxnCov/W0.cure_vs_W2.cure_EllaGeneSets_2025.11.05/W2_cure.MTb.MetaGeneSets.UP.txt")
# EllaGeneSets_2025.11.05_W2.relapse.ComparedTo.W2.cure <- read.delim("Data/DE_Run1to4_60TxnCov/W2.cure_vs_W2.relapse_EllaGeneSets_2025.11.05/W2_relapse.MTb.MetaGeneSets.UP.txt") # THIS ONE NEEDS TO BE RE-DONE! Still not working as of 1/12/26

