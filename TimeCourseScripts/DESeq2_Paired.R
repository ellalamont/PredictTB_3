# W0->W2 Using DESeq2 paired
# E. Lamont
# 3/18/26



# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
# https://www.biostars.org/p/437279/

# Update design formula to be ~ subject + condition

source("Import_data.R")

library(DESeq2)

################################################
################# PREPARE DATA #################

# Keep only the Patients that are in W0 and W2
BothWk_PIDs <- c("P_11011", "P_11012", "P_11031", "P_12008", "P_12010",
                 "P_12019", "P_12020", "P_12073", "P_12084", "P_13012",
                 "P_13016", "P_13041", "P_14005")
BothWk_pipeSummary <- GoodSputum60_pipeSummary %>% 
  filter(Patient %in% BothWk_PIDs)
BothWk_SampleList2 <- BothWk_pipeSummary %>% pull(SampleID2)
BothWk_RawReadsf <- GoodSamples60_RawReadsf %>% dplyr::select(X, all_of(BothWk_SampleList2)) %>%
  column_to_rownames("X") %>%
  as.matrix()

mode(BothWk_RawReadsf) <- "integer"

# Remove genes with <5 reads
keep <- rowSums(BothWk_RawReadsf > 5) >= 0.5 * ncol(BothWk_RawReadsf)
BothWk_RawReadsf <- BothWk_RawReadsf[keep, ] # now only 3950 genes (instead of 4499)

# Generate the metadata
my_metadata <- my_pipeSummary %>%
  filter(SampleID2 %in% colnames(BothWk_RawReadsf)) %>%
  dplyr::select(SampleID2, Outcome, Patient, Week) %>%
  column_to_rownames("SampleID2")

# Make sure everything is in the same order
my_metadata <- my_metadata[colnames(BothWk_RawReadsf), ]

# Make variables factors
my_metadata$Outcome <- factor(my_metadata$Outcome)
my_metadata$Week <- factor(my_metadata$Week, levels = c("Week 0", "Week 2"))
my_metadata$Patient <- factor(my_metadata$Patient)

################################################
############ GENERATE DESEQ2 OBJECT ############

# Ensure sample names line up
stopifnot(all(colnames(BothWk_RawReadsf) == rownames(my_metadata)))
stopifnot(all(rownames(my_metadata) == colnames(BothWk_RawReadsf)))

# Make the DESeqDataSet object (using all the data)
dds <- DESeqDataSetFromMatrix(countData = BothWk_RawReadsf,
                              colData = my_metadata,
                              design = ~ Patient + Week + Week:Outcome)

# Error in checkFullRank(modelMatrix) : 
#   the model matrix is not full rank, so the model cannot be fit as specified.
# One or more variables or interaction terms in the design formula are linear
# combinations of the others and must be removed.
# 
# Please read the vignette section 'Model matrix not full rank':
#   
#   vignette('DESeq2')


# *** This method doesn't work because the Patients only have one outcome... I think...
