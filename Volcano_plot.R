# Make a volcano plot 
# 10/13/25

source("Import_DEG_sets.R") # list_dfs_f2
source("Function_Volcano.R")


# 3/20/26: Try with the volcano function I made
# 4/1/26 Trying again with new function
W2_CvsR_plot <- makeVolcano_EL(list_dfs_f2$W2.cure_vs_W2.relapse, "W2.cure_vs_W2.relapse", 
                               DE_limit = 2, x_max = 6.2, x_min = -6.2)
W2_CvsR_plot
ggsave(plot = W2_CvsR_plot,
       file = paste0("W2.cure_vs_W2.relapse", "_PacTB_v2.png"),
       # path = "Figures/Volcano_60TxnCov/PacTB2026",
       dpi = 600,
       width = 7, height = 5, units = "in")

W0_CvsR_plot <- makeVolcano_EL(list_dfs_f2$W0.cure_vs_W0.relapse, "W0.cure_vs_W0.relapse", 
                               DE_limit = 2, x_max = 6.2, x_min = -6.2)
W0_CvsR_plot
ggsave(W0_CvsR_plot,
       file = paste0("W0.cure_vs_W0.relapse", "_PacTB_v2.png"),
       # path = "Figures/Volcano_60TxnCov/PacTB2026",
       dpi = 600,
       width = 7, height = 5, units = "in")

W0_vs_Ra_plot <- makeVolcano_EL(list_dfs_f2$W0_vs_Ra, "W0_vs_Ra", 
                                DE_limit = 2, x_max = 6.2, x_min = -6.2)
W0_vs_Ra_plot
ggsave(W0_vs_Ra_plot,
       file = paste0("W0_vs_Ra", "_PacTB_v2.png"),
       # path = "Figures/Volcano_60TxnCov/PacTB2026",
       dpi = 600,
       width = 7, height = 5, units = "in")

W2_vs_Ra_plot <- makeVolcano_EL(list_dfs_f2$W2_vs_Ra, "W2_vs_Ra", 
                                DE_limit = 2, x_max = 6.2, x_min = -6.2)
W2_vs_Ra_plot
ggsave(W2_vs_Ra_plot,
       file = paste0("W2_vs_Ra", "_PacTB_v2.png"),
       # path = "Figures/Volcano_60TxnCov/PacTB2026",
       dpi = 600,
       width = 7, height = 5, units = "in")


W0_CvsR_subset_plot <- makeVolcano_EL(list_dfs_f2$W0.Relapse_vs_Cure_W2Subset, "W0.Relapse_vs_Cure_W2Subset", DE_limit = 2, x_max = 6.2, x_min = -6.2)
W0_CvsR_subset_plot
ggsave(W0_CvsR_subset_plot,
       file = paste0("W0.cure_vs_W0.relapse_W2Subset", "_v1.pdf"),
       path = "Figures/Volcano_60TxnCov",
       width = 7, height = 5, units = "in")





# my_path <- "Figures/Volcano_60TxnCov/Log2Fold1_AVG_PVALUE"
# for (i in 1:length(list_dfs_f2)) { ## USING FILTERED DATA ##
#   current_df_name <- df_names[i]
#   filename <- paste0(current_df_name, "_f_AVG_PVALUE_Run1to4.pdf")
#   my_plot <- make_volcano_function_ogP(list_dfs_f2[[i]], df_names[i], 1) ## USING FILTERED DATA ##
#   ggsave(my_plot,
#          file = filename,
#          path = my_path,
#          width = 7, height = 5, units = "in")
# }
