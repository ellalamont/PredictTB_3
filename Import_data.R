# Import data for Run1-4 of the PredictTB samples
# 9/18/25, 10/23/25, 1/8/26

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction
library(stringr)
library(readxl) # To import excel files as dataframes

# DuffyTools
# library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
# library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

## PREDICTTB_RUN1
# This has been edited to include more metadata!
Run1_pipeSummary <- read.csv("Data/PredictTB_Run1/Pipeline.Summary.Details.csv") 

# Just get the samples I want 
Run1_pipeSummary <- Run1_pipeSummary %>% 
  filter(Type %in% c("THP1 spiked", "Week 0 sputum", "Week 2 sputum")) %>% 
  mutate(Run = "PredictTB_Run1")

# BROTH DATA FROM PROBETEST 5
ProbeTest5_pipeSummary <- read.csv("Data/ProbeTest5/ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") 
Broth_pipeSummary <- ProbeTest5_pipeSummary %>% 
  filter(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9", "THP1_1e6_1a_S28"))

## PREDICTTB_RUN2
# This has been edited to include more metadata!
Run2_pipeSummary <- read.csv("Data/PredictTB_Run2/Pipeline.Summary.Details.csv") 
Run2_pipeSummary <- Run2_pipeSummary %>%
  mutate(Run = "PredictTB_Run2")

## PREDICTTB_RUN2.5
# This has been edited to include more metadata!
Run2.5_pipeSummary <- read.csv("Data/PredictTB_Run2_5/Pipeline.Summary.Details.csv") 
Run2.5_pipeSummary <- Run2.5_pipeSummary %>%
  mutate(Run = "PredictTB_Run2.5")

## PREDICTTB_RUN3
# This has been edited to include more metadata!
Run3_pipeSummary <- read.csv("Data/PredictTB_Run3/Pipeline.Summary.Details.csv") 
Run3_pipeSummary <- Run3_pipeSummary %>%
  mutate(Run = "PredictTB_Run3")

## PREDICTTB_RUN4
# This has been edited to include more metadata!
Run4_pipeSummary <- read.csv("Data/PredictTB_Run4/Pipeline.Summary.Details.csv") 
Run4_pipeSummary <- Run4_pipeSummary %>%
  mutate(Run = "PredictTB_Run4")

# Merge the pipeSummaries
All_pipeSummary <- merge(Run1_pipeSummary, Run2_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Broth_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Run2.5_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Run3_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Run4_pipeSummary, all = T)
# Merge two columns
All_pipeSummary <- All_pipeSummary %>% 
  mutate(Type = coalesce(Type, Sample_Type)) %>%
  mutate(Type2 = coalesce(Type2, Sample_Type)) %>% 
  mutate(Type2 = str_replace(Type2, "^THP1$", "THP1 spiked")) %>% 
  dplyr::select(-Sample_Type) %>%
  dplyr::select(-c(N_Splice, P_Splice, Replicate, RT, CFU_per_g, CFU_per_mL, Ra_cells, EukrRNADep, Hyb_Time, Probe, Probe_ng, Pooled_Set, X, mRNA_ng, ct, ttd)) %>%
  filter(SampleID != "Undetermined_S0")

# Make a second SampleID column
# Remove _S* From names
All_pipeSummary$SampleID2 <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
All_pipeSummary$Run2 <- gsub(x = All_pipeSummary$Run, pattern = "PredictTB_", replacement = "")
All_pipeSummary <- All_pipeSummary %>% mutate(SampleID2 = paste0(Run2, "_", SampleID2))

# Add a % of genes with >= 10 reads aligning column
All_pipeSummary <- All_pipeSummary %>%
  mutate(Txn_Coverage = round(AtLeast.10.Reads/4499*100))


###########################################################
############### IMPORT AND PROCESS TPM VALUES #############

Run1_tpm <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.TPM.csv")
Run1_tpm <- Run1_tpm %>% 
  dplyr::select(X, THP1_1e6_1_S67, contains("W")) %>%
  rename_with(~ paste0("Run1_", .), -X) # Add Run1 to the beginning of every column because some have the same name

# Just pull the tpm of the THP1 spiked from another run: THP1 1e6_1 (Predict rack 2 box 1 I04)
# Need THP1 1e6_1a from the Januaray run. Also need 
ProbeTest5_tpm <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") 
ProbeTest5_tpm_Broth <- ProbeTest5_tpm %>% 
  dplyr::select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28) %>%
  rename_with(~ paste0("ProbeTest5_", .), -X) 

Run2_tpm <- read.csv("Data/PredictTB_Run2/Mtb.Expression.Gene.Data.TPM.csv")
Run2_tpm <- Run2_tpm %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run2_", .), -X) 

Run2.5_tpm <- read.csv("Data/PredictTB_Run2_5/Mtb.Expression.Gene.Data.TPM.csv")
Run2.5_tpm <- Run2.5_tpm %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run2.5_", .), -X) 

Run3_tpm <- read.csv("Data/PredictTB_Run3/Mtb.Expression.Gene.Data.TPM.csv")
Run3_tpm <- Run3_tpm %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run3_", .), -X) 

Run4_tpm <- read.csv("Data/PredictTB_Run4/Mtb.Expression.Gene.Data.TPM.csv")
Run4_tpm <- Run4_tpm %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run4_", .), -X) 

All_tpm <- merge(Run1_tpm, ProbeTest5_tpm_Broth, all = T)
All_tpm <- merge(All_tpm, Run2_tpm, all = T)
All_tpm <- merge(All_tpm, Run2.5_tpm, all = T)
All_tpm <- merge(All_tpm, Run3_tpm, all = T)
All_tpm <- merge(All_tpm, Run4_tpm, all = T)

# Remove the _S at the end
names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)


###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% 
  dplyr::select(X, THP1_1e6_1_S67, contains("W")) %>%
  dplyr::select(-contains("NoRT")) %>%
  rename_with(~ paste0("Run1_", .), -X) # Add Run1 to the beginning of every column because some have the same name

ProbeTest5_RawReads <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.readsM_moreTrim.csv")
ProbeTest5_RawReads_Broth <- ProbeTest5_RawReads %>% 
  dplyr::select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28) %>%
  rename_with(~ paste0("ProbeTest5_", .), -X) 

Run2_RawReads <- read.csv("Data/PredictTB_Run2/Mtb.Expression.Gene.Data.readsM.csv")
Run2_RawReads <- Run2_RawReads %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run2_", .), -X) 

Run2.5_RawReads <- read.csv("Data/PredictTB_Run2_5/Mtb.Expression.Gene.Data.readsM.csv")
Run2.5_RawReads <- Run2.5_RawReads %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run2.5_", .), -X) 

Run3_RawReads <- read.csv("Data/PredictTB_Run3/Mtb.Expression.Gene.Data.readsM.csv")
Run3_RawReads <- Run3_RawReads %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run3_", .), -X) 

Run4_RawReads <- read.csv("Data/PredictTB_Run4/Mtb.Expression.Gene.Data.readsM.csv")
Run4_RawReads <- Run4_RawReads %>% 
  dplyr::select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run4_", .), -X) 
  

# Merge the RawReads I collected above
All_RawReads <- merge(Run1_RawReads, Run2_RawReads, all = T)
All_RawReads <- merge(All_RawReads, ProbeTest5_RawReads_Broth)
All_RawReads <- merge(All_RawReads, Run2.5_RawReads)
All_RawReads <- merge(All_RawReads, Run3_RawReads)
All_RawReads <- merge(All_RawReads, Run4_RawReads)

# Remove the _S at the end
names(All_RawReads) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_RawReads))

# # Just keep the samples passing filter
# All_RawReads <- All_RawReads %>% dplyr::select("X", all_of(GoodSampleList), "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
# 
# 
# 
###########################################################
################ MAKE TPM FROM Rv RAW READS ###############

# Bob's TPM includes all the non-coding RNAs, make a new TPM from just the protein-coding genes

# Keep only the protein coding Rv genes
All_RawReads_f <- All_RawReads %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))

source("Function_CalculateTPM.R")
All_tpm_f <- CalculateTPM_RvOnly(All_RawReads_f)

# Remove the _S at the end
# names(All_tpm_f) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_tpm_f))


###########################################################
#################### BATCH CORRECTION #####################
# Using CombatSeq

# count_matrix <- as.matrix(All_RawReads_f %>% column_to_rownames("X"))
# # Ensure integer counts
# mode(count_matrix) <- "integer"
# 
# # Reorder metadata to match column order in count_matrix
# meta <- All_pipeSummary[match(colnames(count_matrix), All_pipeSummary$SampleID2), ]
# 
# # Check alignment
# all(meta$SampleID2 == colnames(count_matrix))  # should be TRUE
# 
# # IF NOT TRUE RUN THESE:
# ## This will show the samples in count_matrix that don't match metadata
# # colnames(count_matrix)[!colnames(count_matrix) %in% All_pipeSummary$SampleID2]
# ## And the opposite: metadata samples not in count_matrix
# # All_pipeSummary$SampleID2[!All_pipeSummary$SampleID2 %in% colnames(count_matrix)]
# 
# 
# # Extract batch (run) and condition (for checking later)
# batch <- meta$Run
# condition <- meta$Type2
# 
# # Run ComBat-Seq
# combat_counts <- ComBat_seq(
#   count_matrix,
#   batch = batch,
#   group = condition # optional, helps preserve biological signal
# )
# 
# # Now convert to TPM
# All_RawReads_f.bc <- as.data.frame(combat_counts) %>% rownames_to_column("X")
# All_tpm_f.bc <- CalculateTPM_RvOnly(All_RawReads_f.bc)

###########################################################
######## CALCULATE TXN COVERAGE FROM Rv GENES ONLY ########

# Count, for each column (sample), how many genes have >= 10 reads
NumGoodReads <- colSums(All_RawReads_f >= 10)

# Add as new column in All_pipeSummary, matching by SampleID2
All_pipeSummary$AtLeast.10.Reads_f <- NumGoodReads[All_pipeSummary$SampleID2]

# Add transcriptional coverage
All_pipeSummary <- All_pipeSummary %>% mutate(Txn_Coverage_f = round(AtLeast.10.Reads_f/4030*100))

###########################################################
#################### ADD ARM INFORMATION ##################
# 12/15/25

source("Import_SampleMetadata.R")

All_pipeSummary <- All_pipeSummary %>% select(-Arm) %>%
  left_join(my_metadata %>% 
              filter(Visit == "Day 0") %>% 
              distinct(Patient, .keep_all = TRUE), 
            select(Patient, Arm), by = "Patient")

###########################################################
########### LIST OF SAMPLES PASSING INSPECTION ############

# Using the transcriptional coverage from Rv genes only! 

# Remove the duplicates (these were hand picked)
my_pipeSummary <- All_pipeSummary %>% 
  filter(!SampleID2 %in% c("Run1_W0_12024", "Run3_W0_12029", "Run1_W0_12043", "Run1_W0_12082", "Run2_W2_13016", "Run1_W0_13026", "Run3_W0_13051", "Run2.5_W0_13051", "Run2_W0_13094", "Run1_W0_14051", "Run1_W0_14136", "Run1_W0_15072", "Run2_W0_15083", "Run2_W2_11058", "Run3_W2_12008", "Run2_W2_12010", "Run3_W2_12012", "Run3_W2_12029", "Run2_W2_12032", "Run3_W2_12043", "Run3_W2_13026_A","Run3_W2_13026_B", "Run3_W2_13045", "Run2_W2_14113", "Run2_W2_14136", "Run3_W2_15017", "Run2_W2_15029", "Run2_W2_15065", "Run2_W2_15089", "Run2_W2_12007", "Run4_W2_12008", "Run4_W2_12012", "Run4_W2_12019", "Run3_W2_12020", "Run3_W2_12025", "Run2_W2_12028", "Run4_W2_12029", "Run4_W2_12052", "Run3_W2_13001", "Run4_W0_13001", "Run4_W2_13026", "Run4_W2_13051", "Run4_W2_14025", "Run4_W2_14051", "Run4_W2_14136")) %>% 
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome))

# Know broth is good
BrothSampleList <- my_pipeSummary %>% 
  filter(str_detect(SampleID, "Broth")) %>%
  pull(SampleID)

# With 60% Transcriptional Coverage (11/3/25)
# 1/8/26: Removing the >1 M reads requirement
GoodSampleList60 <- my_pipeSummary %>%
  filter(Txn_Coverage_f >= 60) %>% 
  pull(SampleID2)
SputumSampleList60 <- GoodSampleList60[grep("W", GoodSampleList60)] # 59 as of Run4
GoodSamples60_pipeSummary <- my_pipeSummary %>% 
  filter(SampleID2 %in% GoodSampleList60) 
GoodSamples60_tpmf <- All_tpm_f %>% dplyr::select(all_of(GoodSampleList60))
GoodSamples60_RawReadsf <- All_RawReads_f %>% dplyr::select(X, all_of(GoodSampleList60))



###########################################################
############# TPM WITH <10TPM GENES REMOVED ###############
# 11/4/25: Remove the genes with <10 TPM in any sample because they cannot be trusted!!
# Need to save a list of this so I can remove them from other datasets

# Not done yet!!

###########################################################
##################### SUMMARY NUMBERS #####################
my_pipeSummary %>% 
  filter(!is.na(Patient)) %>% # 123 total sputum samples
  # filter(N_Genomic >= 1000000) %>% 
  filter(Txn_Coverage_f >=60) %>%
  group_by(Week, Outcome) %>%
  summarize(N_samples = n())

my_pipeSummary %>% 
  filter(!is.na(Patient)) %>% # 123 total sputum samples
  # filter(N_Genomic >= 1000000) %>% 
  # filter(Txn_Coverage_f >=60) %>%
  group_by(Week, Outcome, main_lineage) %>%
  summarize(N_samples = n())

my_pipeSummary %>% 
  filter(!is.na(Patient)) %>% # 123 total sputum samples
  # filter(N_Genomic >= 1000000) %>% 
  filter(Txn_Coverage_f >=60) %>%
  group_by(Week, Outcome, Arm) %>%
  summarize(N_samples = n())


###########################################################
################### EXPORT METADATA INFO ##################
# 1/29/26 for Adding to the DEG code on the lenovo

GoodSputum60_pipeSummary <- GoodSamples60_pipeSummary %>% filter(SampleID2 %in% colnames(GoodSputum60_tpmf))


###########################################################
##################### EXPORT TPM INFO #####################
# 1/29/26 for Jason Yang

GoodSputum60_tpmf <- GoodSamples60_tpmf %>% select(contains("W"))
GoodSputum60_pipeSummary <- GoodSamples60_pipeSummary %>% filter(SampleID2 %in% colnames(GoodSputum60_tpmf))

# Clean up the names
names(GoodSputum60_tpmf) <- sub("^Run[0-9]+_", "", names(GoodSputum60_tpmf))
GoodSputum60_pipeSummary$SampleID <- gsub(pattern = "_S[0-9]+$", replacement = "", x = GoodSputum60_pipeSummary$SampleID)
GoodSputum60_Metadata <- GoodSputum60_pipeSummary %>% select(SampleID, Week, Patient, Outcome, Type2, Arm, main_lineage)

# Remove the failure sample (12025) and the cousin relapse (13001, 13026)
GoodSputum60_tpmf <- GoodSputum60_tpmf %>% select(-any_of(c("W0_12025", "W0_13001", "W2_13001", "W0_13026", "W2_13026")))
GoodSputum60_Metadata <- GoodSputum60_Metadata %>% filter(!SampleID %in% c("W0_12025", "W0_13001", "W2_13001", "W0_13026", "W2_13026"))

# Check the tpm and metadata have the same samples
stopifnot(setequal(GoodSputum60_Metadata$SampleID, colnames(GoodSputum60_tpmf)))

# Save as csv
write.csv(GoodSputum60_tpmf, "Data/ForJasonYang/SputumTPM_20260129.csv")
write.csv(GoodSputum60_Metadata, "Data/ForJasonYang/SputumMetadata_20260129.csv", row.names = F)
            