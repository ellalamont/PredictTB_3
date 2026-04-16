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
# list_GeneSets_names

CurrentComparison <- "W0_vs_Ra_TAR"

myData <- paste0(CurrentComparison)

# Need this for facet grouping
Current_GeneSet_df <- read.csv("Data/GeneSet_Data/TAR_Poonawala2024_GeneSets.csv")

# Add some useful columns
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>% 
  filter(!str_detect(PathName, " 24hr")) # Don't look at the >24hr set

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 4, units = "in")

###########################################################
##################### TAR: W2 vs Ra #######################
# list_GeneSets_names

CurrentComparison <- "W2_vs_Ra_TAR"

myData <- paste0(CurrentComparison)

# Need this for facet grouping
Current_GeneSet_df <- read.csv("Data/GeneSet_Data/TAR_Poonawala2024_GeneSets.csv")

# Add some useful columns
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>% 
  filter(!str_detect(PathName, " 24hr")) # Don't look at the >24hr set

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 4, units = "in")

###########################################################
##################### TAR: W2 vs W0 #######################
# list_GeneSets_names

CurrentComparison <- "W2_vs_W0_TAR"

myData <- paste0(CurrentComparison)

# Need this for facet grouping
Current_GeneSet_df <- read.csv("Data/GeneSet_Data/TAR_Poonawala2024_GeneSets.csv")

# Add some useful columns
myBubble_df <- list_GeneSets2[[myData]] %>%
  left_join(Current_GeneSet_df %>% # Add the Group names
              dplyr::rename(PathName = GeneSet) %>%
              dplyr::select(PathName, Group) %>% 
              distinct(PathName, .keep_all = TRUE),
            by = "PathName") %>% 
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>% 
  filter(!str_detect(PathName, " 24hr")) # Don't look at the >24hr set

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 4, units = "in")



###########################################################
################# EllaGeneSets: W0 vs Ra ##################
# 4/15/26: Updated gene set

CurrentGeneSet <- "EllaGeneSets_2026.04.15"
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

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
# myBubble_df_PacTB <- myBubble_df %>%
#   filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
#                       "Growth", "Hypoxia related", # "Metal", 
#                       "Nutrient Starvation", "Stress responses"))
# myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData, 
#                                      x_min = -3.5, x_max = 3.5)
# myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")

### 4/15/26: Subset for TB Symposium ###
myBubble_df_TBSymp <- myBubble_df %>%
  filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
                      "Growth", "Respiration", "Hypoxia",
                      "Stress responses")) %>%
  filter(!PathName %in% c("Non-specific stress responses", "Oxidative stress response")) %>% 
  mutate(PathName_2 = ifelse(PathName_2 == "Non-specific stress responses Full (n=64)", "Non-specific stress responses (n=64)", PathName_2))
myBubble_plot_TBSymp <- makeBubble_EL(myBubble_df_TBSymp, title = myData,
                                     x_min = -3.5, x_max = 3.5)
myBubble_plot_TBSymp
# ggsave(myBubble_plot_TBSymp,
#        file = paste0(myData, "_TBSymp_v2.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 6.2, units = "in")


###########################################################
################# EllaGeneSets: W2 vs Ra ##################
# 4/15/26: Updated gene set

CurrentGeneSet <- "EllaGeneSets_2026.04.15"
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

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")


### Subset for PacTB ###
# myBubble_df_PacTB <- myBubble_df %>%
#   filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
#                       "Growth", "Hypoxia related", # "Metal", 
#                       "Nutrient Starvation", "Stress responses"))
# myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData, 
#                                      x_min = -3.5, x_max = 3.5)
# myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")

### 4/15/26: Subset for TB Symposium ###
myBubble_df_TBSymp <- myBubble_df %>%
  filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
                      "Growth", "Respiration", "Hypoxia",
                      "Stress responses")) %>%
  filter(!PathName %in% c("Non-specific stress responses", "Oxidative stress response")) %>% 
  mutate(PathName_2 = ifelse(PathName_2 == "Non-specific stress responses Full (n=64)", "Non-specific stress responses (n=64)", PathName_2))
myBubble_plot_TBSymp <- makeBubble_EL(myBubble_df_TBSymp, title = myData,
                                      x_min = -3.5, x_max = 3.5)
myBubble_plot_TBSymp
# ggsave(myBubble_plot_TBSymp,
#        file = paste0(myData, "_TBSymp_v2.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 6.2, units = "in")


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

# myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
# myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
myBubble_df_PacTB <- myBubble_df %>%
  filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
                      "Growth", "Hypoxia related", "Metal", 
                      "Nutrient Starvation", "Stress responses"))
myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData)
myBubble_plot_PacTB


###########################################################
############ EllaGeneSets: W0 CURE vs RELAPSE #############
# Updated 4/15/26
# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.04.15"
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

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
# myBubble_df_PacTB <- myBubble_df %>%
#   filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
#                       "Growth", "Hypoxia related", # "Metal", 
#                       "Nutrient Starvation", "Stress responses"))
# myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData, 
#                                      x_min = -3.5, x_max = 3.5)
# myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")



###########################################################
############# EllaGeneSets: W2 CURE vs RELAPSE ############
# Updated 4/15/26
# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.04.15"
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

myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
# myBubble_df_PacTB <- myBubble_df %>%
#   filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
#                       "Growth", "Hypoxia related", # "Metal", 
#                       "Nutrient Starvation", "Stress responses"))
# myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData, 
#                                      x_min = -3.5, x_max = 3.5)
# myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")


### 4/15/26: Subset for TB Symposium ###
myBubble_df_TBSymp <- myBubble_df %>%
  filter(Group %in% c("Toxin/Antitoxin", "ESX genes", "Growth",
                      "Hypoxia related", "Nucleic Acid",
                      "Stress responses")) %>%
  filter(!PathName %in% c("Non-specific stress responses", "Oxidative stress response")) %>% 
  mutate(PathName_2 = ifelse(PathName_2 == "Non-specific stress responses Full (n=64)", "Non-specific stress responses (n=64)", PathName_2))
myBubble_plot_TBSymp <- makeBubble_EL(myBubble_df_TBSymp, title = myData,
                                      x_min = -2, x_max = 1.5)
myBubble_plot_TBSymp
ggsave(myBubble_plot_TBSymp,
       file = paste0(myData, "_TBSymp_v1.png"),
       path = "Figures/Bubbles",
       dpi = 600,
       width = 6, height = 4.2, units = "in")


###########################################################
############### EllaGeneSets: W2 CURE VS Ra ###############

# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W2.Cure_vs_Ra"

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

# myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
# myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
myBubble_df_PacTB <- myBubble_df %>%
  filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
                      "Growth", "Hypoxia related", # "Metal", 
                      "Nutrient Starvation", "Stress responses"))
myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData)
myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")

###########################################################
############# EllaGeneSets: W2 RELAPSE VS Ra ##############

# list_GeneSets_names

CurrentGeneSet <- "EllaGeneSets_2026.03.19"
CurrentComparison <- "W2.Relapse_vs_Ra"

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

# myBubble_plot <- makeBubble_EL(myBubble_df, title = myData)
# myBubble_plot
# ggsave(myBubble_plot,
#        file = paste0(myData, "_v1.pdf"),
#        path = "Figures/Bubbles",
#        width = 6, height = 10, units = "in")

### Subset for PacTB ###
myBubble_df_PacTB <- myBubble_df %>%
  filter(Group %in% c("Cell wall synthesis and remodeling", "Fatty Acid and Cholesterol",
                      "Growth", "Hypoxia related", # "Metal", 
                      "Nutrient Starvation", "Stress responses"))
myBubble_plot_PacTB <- makeBubble_EL(myBubble_df_PacTB, title = myData)
myBubble_plot_PacTB
# ggsave(myBubble_plot_PacTB,
#        file = paste0(myData, "_PacTB_v1.png"),
#        path = "Figures/Bubbles",
#        dpi = 600,
#        width = 6.5, height = 7, units = "in")
