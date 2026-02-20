# Import Bob's differential gene expression data and the metagenesets data
# 1/12/26
# E. Lamont

source("Import_data.R")

# Data is coming from the Lenovo PredictTB_allRuns
# 1/12/26: Using the 60% Txn coverage cutoff runs 1-4


###########################################################
################### IMPORT BOB's DE DATA ##################
# 2/20/26: updated to be cleaner import

parent_dir <- "Data/DE_Run1to4_60TxnCov"

folders <- list.dirs(parent_dir, full.names = T, recursive = F)
folders <- folders[!grepl("GeneSets", basename(folders))] # Remove the folders containing GeneSets because these are duplicates for the DEG

# Function to get the JOINED.txt files from each folder
list_dfs <- lapply(folders, function(f) {
  
  # Find the JOINED file inside each folder
  current_file <- list.files(f, pattern = "\\.JOINED\\.txt$", full.names = T)
  
  # Read in the file
  if (length(current_file) == 1) {
    read.delim(current_file)
  } else {NULL}
}
  )

# Add the names to the list of lists
names(list_dfs) <- basename(folders)

# Also make a list of the names
df_names <- names(list_dfs)


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
# 2/20/26: Removing the FDR p-values because David says Bob's AVG_PVALUE is good. Doing Log2Fold >1, >2, >2.5

# Make a new list to hold dataframes with extra columns
list_dfs_f2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs_f)) {
  
  current_df <- list_dfs_f[[i]]
  current_df_name <- df_names[i]
  
  # Columns for Log2Fold > 1 and AVG_PVALUE < 0.05
  current_df$DE1 <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE1 <- factor(current_df$DE1, levels = ordered_DE)
  current_df$DE1_labels <- ifelse(current_df$DE1 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 2 and AVG_PVALUE < 0.05
  current_df$DE2 <- ifelse(current_df$LOG2FOLD < -2 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 2 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2 <- factor(current_df$DE2, levels = ordered_DE)
  current_df$DE2_labels <- ifelse(current_df$DE2 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 2.5 and AVG_PVALUE < 0.05
  current_df$DE2.5 <- ifelse(current_df$LOG2FOLD < -2.5 & current_df$AVG_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 2.5 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2.5 <- factor(current_df$DE2.5, levels = ordered_DE)
  current_df$DE2.5_labels <- ifelse(current_df$DE2.5 != "not significant", current_df$GENE_NAME, NA)

  list_dfs_f2[[current_df_name]] <- current_df
}




###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############
# 2/20/26: updated to be cleaner import

parent_dir <- "Data/DE_Run1to4_60TxnCov"
folders <- list.dirs(parent_dir, full.names = T, recursive = F)

# Function to get the UP.txt files from each folder
list_GeneSets <- lapply(folders, function(f) {
  
  # Find the UP file inside each folder
  current_file <- list.files(f, pattern = "\\.UP\\.txt$", full.names = T)
  
  # Read in the file
  if (length(current_file) == 1) {
    read.delim(current_file)
  } else {NULL}
}
)

# Add the names to the list of lists
names(list_GeneSets) <- paste0("GeneSets_", basename(folders))

# Also make a list of the names
list_GeneSets_names <- names(list_GeneSets)
