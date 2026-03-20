# Specific gene sets comparisons I have chosen
# E. Lamont
# 3/19/26


# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result
# I know the UP and DOWN files are a little different, just working with the UP files for now, should maybe ask Bob if that is okay at some point

source("Import_DEG_sets.R")
source("Import_GeneSets.R")
source("Function_Bubbles.R")

###########################################################
##################### TAR: W0 vs Ra #######################

# Need this for facet grouping
Current_GeneSet_df <- read.csv("Data/GeneSet_Data/TAR_Poonawala2024_GeneSets.csv")

# Add some useful columns
myBubble_df <- list_GeneSets$GeneSets_W2_vs_W0_TAR %>%
  mutate(Significance = ifelse(AVG_PVALUE < 0.05, "significant", "not significant"),
         FillColor = case_when(Significance == "significant" & LOG2FOLD > 0 ~ "pos",
                               Significance == "significant" & LOG2FOLD < 0 ~ "neg",
                               TRUE ~ "ns"),
         PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>% 
  filter(!str_detect(PathName, " 24hr")) # Don't look at the >24hr set

myBubble_plot <- makeBubble_EL(list_GeneSets$GeneSets_W2_vs_W0_TAR, myTitle = "GeneSets_W2_vs_W0_TAR")
myBubble_plot


###########################################################
################# EllaGeneSets: W0 vs Ra ##################

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W0_vs_Ra"

myData <- paste0(CurrentComparison, "_", CurrentGeneSet)

# Need this for facet grouping
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))

# Add the Grouping column
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))

myBubble_plot <- makeBubble_EL(myBubble_df, myTitle = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

###########################################################
################# EllaGeneSets: W2 vs Ra ##################

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W2_vs_Ra"

myData <- paste0(CurrentComparison, "_", CurrentGeneSet)

# Need this for facet grouping
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))

# Add the Grouping column
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))

myBubble_plot <- makeBubble_EL(myBubble_df, myTitle = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")


###########################################################
################# EllaGeneSets: W2 vs W0 ##################

# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W2_vs_W0"

myData <- paste0(CurrentComparison, "_", CurrentGeneSet)

# Need this for facet grouping
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))

# Add the Grouping column
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))

myBubble_plot <- makeBubble_EL(myBubble_df, myTitle = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")


###########################################################
################# EllaGeneSets: W2 vs W0 ##################

# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W0.Cure_vs_Relapse"

myData <- paste0(CurrentComparison, "_", CurrentGeneSet)

# Need this for facet grouping
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))

# Add the Grouping column
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))

myBubble_plot <- makeBubble_EL(myBubble_df, myTitle = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")



###########################################################
################# EllaGeneSets: W2 vs W0 ##################

# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W2.Cure_vs_Relapse"

myData <- paste0(CurrentComparison, "_", CurrentGeneSet)

# Need this for facet grouping
Current_GeneSet_df <- read.csv(paste0("Data/GeneSet_Data/", CurrentGeneSet, ".csv"))

# Add the Grouping column
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19))

myBubble_plot <- makeBubble_EL(myBubble_df, myTitle = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")


