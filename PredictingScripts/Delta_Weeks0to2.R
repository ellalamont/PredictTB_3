# Look at the change between Weeks 0 and 2 for patients with both
# E. Lamont
# 3/13/26


source("Import_data.R")


###########################################################
################ PATIENTS WITH BOTH WEEKS #################

BothWk_PIDs <- c("P_11011", "P_11012", "P_11031", "P_12008", "P_12010",
                 "P_12019", "P_12020", "P_12073", "P_12084", "P_13012",
                 "P_13016", "P_13041", "P_14005")

BothWk_pipeSummary <- GoodSputum60_pipeSummary %>% 
  filter(Patient %in% BothWk_PIDs)

BothWk_SampleList <- BothWk_pipeSummary %>% pull(SampleID)

BothWk_tpmf_log2 <- GoodSputum60_tpmf_log2 %>% dplyr::select(all_of(BothWk_SampleList))

BothWk_tpmf <- GoodSputum60_tpmf %>% dplyr::select(all_of(BothWk_SampleList))


###########################################################
######################### ∆ DELTA #########################

BothWk_log2TPMf_delta <- BothWk_tpmf_log2 %>%
  rownames_to_column("Gene") %>% 
  pivot_longer(-Gene, names_to = "SampleID", values_to = "log2TPMf") %>%
  left_join(BothWk_pipeSummary %>% dplyr::select(SampleID, Week, Patient, Outcome),
            by = "SampleID") %>%
  pivot_wider(id_cols = c(Gene, Patient, Outcome), 
              names_from = Week, values_from = log2TPMf) %>%
  dplyr::rename(W0_log2TPMf = "Week 0", W2_log2TPMf = "Week 2") %>%
  mutate(delta = W2_log2TPMf - W0_log2TPMf)

###########################################################
######################## TPM RATIO ########################

options(scipen = 999)
BothWk_TPMf_ratio <- BothWk_tpmf %>%
  rownames_to_column("Gene") %>% 
  pivot_longer(-Gene, names_to = "SampleID", values_to = "TPMf") %>%
  left_join(BothWk_pipeSummary %>% dplyr::select(SampleID, Week, Patient, Outcome),
            by = "SampleID") %>%
  pivot_wider(id_cols = c(Gene, Patient, Outcome), 
              names_from = Week, values_from = TPMf) %>%
  dplyr::rename(W0_TPMf = "Week 0", W2_TPMf = "Week 2") %>%
  mutate(Ratio = (W2_TPMf+0.1) / (W0_TPMf+0.1)) # Had to add 0.1 to avoid inf with dividing by 0


###########################################################
####################### COMBINE DFs #######################

BothWks_df <- BothWk_log2TPMf_delta %>% 
  left_join(BothWk_TPMf_ratio, by = c("Gene", "Patient", "Outcome"))
