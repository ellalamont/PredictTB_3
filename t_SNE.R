# T-SNE
# E. Lamont 
# 4/6/26

# Doesn't look that good...

# https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/
# https://bioconductor.org/packages/devel/bioc/vignettes/M3C/inst/doc/M3Cvignette.pdf


source("Import_data.R")

BiocManager::install("M3C")
library(M3C)


###########################################################
####################### ORGANIZE DATA #####################

# M3C requires
# A matrix or data frame of normalised continuous expression data (e.g. microarray, RNA-seq, methylation arrays, protein arrays) where columns equal samples and rows equal features. M3C’s reference will work better if the feature data is approximately normally distributed across biological replicates.
# For example, for RNA-seq data, VST or rlog transformed count data, log2(CPM), log2(TPM), and
# log2(RPKM), are all acceptable forms of normalisation.
# The data should be filtered to remove features with no or very low signal, and filtered using variance to reduce dimensionality (unsupervised), or p value from a statistical test (supervised). We include a feature filter function that uses variance which is demonstrated later in the vignette. This variance function should be applied after removing the tendency of the variance to increase with the mean (e.g. after VST or log2 transformation).


View(GoodSputum60_tpmf_log2)
GoodSputum60_Metadata


###########################################################
############################ M3C ##########################

pca(GoodSputum60_tpmf_log2, labels = GoodSputum60_Metadata$Type2, legendtextsize = 10,axistextsize = 10,dotsize=2)

tsne(GoodSputum60_tpmf_log2, labels = GoodSputum60_Metadata$Type2)

res <- M3C(GoodSputum60_tpmf_log2, method = 2, objective='PAC',fsize=8, lthick=1, dotsize=1.25)

