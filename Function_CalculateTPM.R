CalculateTPM <- function(Raw_reads) {
  ## This takes the raw reads from Bob's pipeline and turns them into TPM
  ## Requires the MTb.MapSet.rda to be in the R environment
  
  # 1. Extract and organize the gene lengths 
  load("Data/MTb.MapSet.rda")
  my_geneLengths <- mapSet[["geneMap"]] %>% dplyr::select(GENE_ID, NAME, N_EXON_BASES)
  
  my_geneLengths_ordered <- my_geneLengths[match(Raw_reads$X, my_geneLengths$GENE_ID), ]
  my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
  
  # 2. Convert column to rowname in the raw data
  Raw_reads_2 <- Raw_reads %>% column_to_rownames("X")
  
  # 3. Divide the raw reads by the gene lengths in KB
  All_RPK <- Raw_reads_2 / my_geneLengths_ordered$Kilobases
  
  # 4. Generate the "per million" scaling factor (sum of RPK values / 1e6 for each column)
  ScalingFactor <- colSums(All_RPK) / 1e6
  
  # 5. Divide the RPK by the Scaling Factor
  All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")
  
  return(All_tpm)
  
}

CalculateTPM_RvOnly <- function(Raw_reads) {
  ## This takes the raw reads from Bob's pipeline and turns them into TPM
  ## Requires the MTb.MapSet.rda to be in the R environment
  ## The genes are filtered here to only include protein coding Rv genes, input needs to be raw reads filtered in this way also!
  
  # 1. Extract and organize the gene lengths 
  load("Data/MTb.MapSet.rda")
  my_geneLengths <- mapSet[["geneMap"]] %>% dplyr::select(GENE_ID, NAME, N_EXON_BASES)
  
  my_geneLengths_f <- my_geneLengths %>%
    filter(grepl("^Rv[0-9]+[A-Za-z]?$", GENE_ID))
  
  my_geneLengths_ordered <- my_geneLengths_f[match(Raw_reads$X, my_geneLengths_f$GENE_ID), ]
  my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
  
  # 2. Convert column to rowname in the raw data
  Raw_reads_2 <- Raw_reads %>% column_to_rownames("X")
  
  # 3. Divide the raw reads by the gene lengths in KB
  All_RPK <- Raw_reads_2 / my_geneLengths_ordered$Kilobases
  
  # 4. Generate the "per million" scaling factor (sum of RPK values / 1e6 for each column)
  ScalingFactor <- colSums(All_RPK) / 1e6
  
  # 5. Divide the RPK by the Scaling Factor
  All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")
  
  return(All_tpm)
  
}


############# TESTING #############


# load("Data/MTb.MapSet.rda")
# my_geneLengths <- mapSet[["geneMap"]] %>% dplyr::select(GENE_ID, NAME, N_EXON_BASES)
# my_geneLengths_f <- my_geneLengths %>%
#   filter(grepl("^Rv[0-9]+[A-Za-z]?$", GENE_ID))
# my_geneLengths_ordered <- my_geneLengths_f[match(Run2_RawReads$X, my_geneLengths_f$GENE_ID), ]
# my_geneLengths_ordered <- my_geneLengths_f %>% mutate(Kilobases = N_EXON_BASES/1000)
# 
# # GoodBiolSamples_wRv_RawReads
# 
# load("Data/MTb.MapSet.rda")
# my_geneLengths <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, N_EXON_BASES)
# my_geneLengths_ordered <- my_geneLengths[match(All_RawReads_f$X, my_geneLengths$GENE_ID), ]
# my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
# 
# Raw_reads_2 <- All_RawReads_f %>% column_to_rownames("X")
# All_RPK <- Raw_reads_2 / my_geneLengths_ordered$Kilobases
# ScalingFactor <- colSums(All_RPK) / 1e6
# All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")


# 9/1/25 Checking with the batch corrected data
# Raw_reads_2 <- counts_corrected
# All_RPK <- Raw_reads_2 / my_geneLengths_ordered$Kilobases
# ScalingFactor <- colSums(All_RPK) / 1e6
# All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")




# Checking how it looks when I use the raw reads with only the coding Rv genes
# my_geneLengths_ordered2 <- my_geneLengths %>%
#   filter(GENE_ID %in% GoodBiolSamples_RawReads_f$X) %>%
#   arrange(match(GENE_ID, GoodBiolSamples_RawReads_f$X)) %>%
#   mutate(Kilobases = N_EXON_BASES/1000)
# Raw_reads_2 <- GoodBiolSamples_RawReads_f %>% column_to_rownames("X")
# All_RPK2 <- Raw_reads_2 / my_geneLengths_ordered2$Kilobases
# ScalingFactor2 <- colSums(All_RPK2) / 1e6
# All_tpm2 <- sweep(All_RPK2, 2, ScalingFactor2, FUN = "/")


