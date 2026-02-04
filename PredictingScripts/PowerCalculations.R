# Use RNASeqPower package to determine power
# E. Lamont 
# 2/3/26



###########################################################
################ COEFFICIENT OF VARIATION #################
# Need this for estimating power

# Coefficient of Variation
# = Standard Deviation / mean

# https://support.bioconductor.org/p/9135351/
# edgeR::estimateCommonDisp

# I guess I do have technical sequencing replicates... maybe I could use them? Not done...

# Try and estimate CV from W0 RawReads 60%Cov
Subset_RawReadsf <- GoodSamples60_RawReadsf %>% 
  column_to_rownames("X") %>%
  select(any_of(W0SputumSamples60_pipeSummary$SampleID2)) 

# Using the EdgeR CV calcutator
y <- DGEList(counts = Subset_RawReadsf)
dge <- estimateCommonDisp(y, verbose = T)
# Disp = 0.71937 , BCV = 0.8482
sqrt(dge$common.dispersion) # [1] 0.8481564

# What if I do a subset of the W0 samples
View(Subset_RawReadsf[1:5])
y2 <- DGEList(counts = Subset_RawReadsf[1:5])
dge2 <- estimateCommonDisp(y2, verbose = T)
# Disp = 0.82477 , BCV = 0.9082 That's not better

# What if I look at the broth samples
Broth_RawReadsf <- All_RawReads_f %>% 
  column_to_rownames("X") %>%
  select(contains("broth"))
y3 <- DGEList(counts = Broth_RawReadsf)
dge3 <- estimateCommonDisp(y3, verbose = T)
# Disp = 0.00554 , BCV = 0.0744

# What if I look at the THP1 spiked samples
THP1Spike_RawReadsf <- All_RawReads_f %>% 
  column_to_rownames("X") %>%
  select(contains("THP1_1e6"))
y4 <- DGEList(counts = THP1Spike_RawReadsf)
dge4 <- estimateCommonDisp(y4, verbose = T)
# Disp = 0.02081 , BCV = 0.1443 
sqrt(dge4$common.dispersion)

# Just the 30 W0 cure sputum used in the R21
W0Cures_R21_SampleIDs <- W0SputumSamples60_pipeSummary %>% 
  filter(Run != "PredictTB_Run4" & Week == "Week 0" & Outcome == "Cure") %>%
  filter(!SampleID %in% c("W0_15035_S25", "W0_12029_S10"))%>% # Not sure why these were left out but they were
  pull(SampleID2)
Subset_RawReadsf <- GoodSamples60_RawReadsf %>% 
  column_to_rownames("X") %>%
  select(all_of(W0Cures_R21_SampleIDs))
y <- DGEList(counts = Subset_RawReadsf)
dge <- estimateCommonDisp(y, verbose = T)
# Disp = 0.62732 , BCV = 0.792 *** USE THIS ONE!

###########################################################
################## MEAN SEQUENCING DEPTH ##################

mean_depth <- mean(rowMeans(Subset_RawReadsf))
mean_depth
# 1627.035

###########################################################
####################### RNASeqPower #######################

# https://bioconductor.org/packages/release/bioc/manuals/RNASeqPower/man/RNASeqPower.pdf
# https://github.com/royfrancis/shiny-rnaseq-power

# Shiny web app
# https://rnaseq-power.serve.scilifelab.se/app/rnaseq-power

# The paper
# https://www.liebertpub.com/doi/10.1089/cmb.2012.0283

BiocManager::install("RNASeqPower")
library(RNASeqPower)

# Calculating effect size, which is the detectable fold change in expression

rnapower(depth = 10, n = 12, n2 = 200, cv = 0.72, alpha = 0.05, power = 0.8)
# [1] 1.924735 should be the effect size (Relative expression effect)

# Best case:
rnapower(depth = 10, n = 21, n2 = 299, cv = 0.72, alpha = 0.05, power = 0.8)
# [1] 1.64437 effect size

rnapower(depth = 1000, n = 12, n2 = 200, cv = 0.792, alpha = 0.05, power = 0.8)
# [1] 1.934773

rnapower(depth = 1000, n = 12, n2 = 299, cv = 0.792, alpha = 0.05, power = 0.8)
# [1] 1.922787

rnapower(depth = 1000, n = 12, cv = 0.792, alpha = 0.05, power = 0.8)
# [1] 2.475805

rnapower(depth = 1000, n = 21, n2 = 299, cv = 0.792, alpha = 0.05, power = 0.8)
# [1] 1.650879


###########################################################
########################### pwr ###########################

install.packages("pwr")
library(pwr)

# This is less appropriate for RNA-seq work! 

pwr.t2n.test(n1 = 12, n2 = 200, d = NULL, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# t test power calculation 
# 
# n1 = 12
# n2 = 200
# d = 0.8364906
# sig.level = 0.05
# power = 0.8
# alternative = two.sided

pwr.t.test(n = NULL, d = 2, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# Two-sample t test power calculation 
# 
# n = 5.089995
# d = 2
# sig.level = 0.05
# power = 0.8
# alternative = two.sided
# 
# NOTE: n is number in *each* group

pwr.t.test(n = NULL, d = 1, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# Two-sample t test power calculation 
# 
# n = 16.71472
# d = 1
# sig.level = 0.05
# power = 0.8
# alternative = two.sided
# 
# NOTE: n is number in *each* group
