

add_DE_columns_function <- function(my_df, log2fold_cutoff = 1, padj_cutoff = 0.05) {
  
  ## Adds columns to DEG dataframe results for plotting with the volcano plot function
  ## Needs to have already run rownames_to_column("gene")
  
  # Determine column names
  DE_col_name <- paste0("DE", log2fold_cutoff)
  label_col_name <- paste0(DE_col_name, "_labels")
  
  # Make the df with new columns
  my_df[[DE_col_name]] <- ifelse(my_df$log2FoldChange < -log2fold_cutoff & my_df$padj < padj_cutoff, "significant down",
                                 ifelse(my_df$log2FoldChange > log2fold_cutoff & my_df$padj < padj_cutoff, "significant up", "not significant"))
  
  # Make the DE significance a factor
  my_df[[DE_col_name]] <- factor(my_df[[DE_col_name]], levels = c("significant down", "not significant", "significant up"))
  
  # Add gene name column
  my_df[[label_col_name]] <- ifelse(my_df[[DE_col_name]] != "not significant", my_df$gene, NA)
  
  return(my_df)
}