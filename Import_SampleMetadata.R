# Import all the data sheets I got via Email about the sputum samples
# 7/29/25

# source("Import_data.R")

library(readxl)

###########################################################
####################### IMPORT DATA #######################

# Import the first data specific to the samples we have
Sample_list <- read_excel("Data/Sample_Metadata/Predict TRIzol sutum sample list_UW.xlsx")
Sample_list <- Sample_list %>% select(`Patient ID`, Visit, `Aliquot Number`, Box, Position, Comments)
names(Sample_list)[names(Sample_list) == "Patient ID"] <- "PatientID"
names(Sample_list)[names(Sample_list) == "Aliquot Number"] <- "AliquotNumber"

# Make a SUBJID column to match with the other datashets
Sample_list$SUBJID <- Sample_list$PatientID
Sample_list$SUBJID <- gsub("PD-", "", Sample_list$SUBJID)
Sample_list$SUBJID <- gsub("-", "", Sample_list$SUBJID)

# Import the data with outcome and arm and merge
Outcome_Arm <- read.csv("Data/Sample_Metadata/Predict trizol sample set full PIDs with arm outcome.csv")
metadata <- merge(Sample_list, Outcome_Arm[,2:4], by = "SUBJID", all = TRUE)

# Import the datasheet with extra values (doesn't look like this is all the information for all the samples I have... just the samples used in the subMIC experiments)
Extra_values <- read_excel("Data/Sample_Metadata/Predict trizol sample subMIC set with outcome arm BL and PETCT.xlsx", sheet = "Predict trizol sample subMIC se")

metadata_2 <- merge(metadata, Extra_values, by = c("SUBJID", "outcome", "arm"), all = TRUE)

# Import the lineage information
Lineage_info <- read_excel("Data/Sample_Metadata/predict_148_17NOV2023_TBprofilerResults_TS.xlsx", sheet = "Summary")
Lineage_info$SUBJID <- gsub("_.*$", "", Lineage_info$sample) # New column with any characters after (.*) to end of sring ($) removed 
Lineage_info$Collection_TimePoint <- gsub("^.*?_", "", Lineage_info$sample)

metadata_3 <- merge(metadata_2, Lineage_info, by = c("SUBJID"), all = TRUE)

# Import the sub breakpoint MICs for RIF and INH
subMIC_info <- read_excel("Data/Sample_Metadata/Sub-breakpoint MIC_20260126_YP.xlsx")
subMIC_info$SUBJID <- gsub("PD-", "", subMIC_info$`Patient ID`)
subMIC_info$SUBJID <- gsub("-", "", subMIC_info$SUBJID)

metadata_4 <- merge(metadata_3, subMIC_info, by = c("SUBJID"), all = TRUE)

###########################################################
##################### EXPERIMENTAL SET ####################
# Want to look at the data that I understand and with the feasibility test set removed

Exp_metadata <- metadata_4 # 1516 rows

# Add a week column
Exp_metadata$Week <- Exp_metadata$Visit
Exp_metadata$Week <- gsub("Week ", "", Exp_metadata$Week)
Exp_metadata$Week <- gsub("Day ", "", Exp_metadata$Week)

# Remove the feasibility test test
Exp_metadata_2 <- Exp_metadata %>% filter(Comments != ("feasibility test set")) # 1311 rows

# Remove the day 0 and week 2 extras (we asked for ALL the week 0 and week 2 samples, but not all patients were processed through the study in the subMIC dataset)
Exp_metadata_3 <- Exp_metadata_2 %>% filter(Age != "NA")

# Add a patient column to match the sequencing
Exp_metadata_3 <- Exp_metadata_3 %>% mutate(Patient = paste0("P_", SUBJID)) # 1023 rows


###########################################################
##################### SUBSET COLUMNS ######################

my_metadata <- Exp_metadata_3 %>% select(outcome, arm, Visit, Age, SEX, TBprev, BMI, Weight, CurrentSmoker, PrevSmoker, SmokeDuration, TTD, XpertCT_wk0, Patient, main_lineage, sub_lineage) %>% rename(Arm = arm)






