# WCGNA
# E. Lamont
# 4/28/26


# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
# https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/


source("Import_data.R")

#install.packages('BiocManager')
library(BiocManager)
#BiocManager::install('WGCNA')
#BiocManager::install('flashClust')

library(WGCNA)
library(flashClust)
library(curl)


###########################################################
########################### DATA ##########################

# Datasheets that still contain W0 and W2, need to keep just the W0 as we go along
# GoodSputum60_tpmf
# GoodSputum60_Metadata

# The tutorial uses gene expression data relative to a control... will try with TPM
# Want samples as rows and genes as columns
expression.data <- as.data.frame(t(GoodSputum60_tpmf %>% select(contains("W0"))))



###########################################################
##################### REMOVE OUTLIERS #####################

# Identifying Outlier Genes
gsg <-goodSamplesGenes(expression.data)
summary(gsg)
# Length Class  Mode   
# goodGenes   4030   -none- logical
# goodSamples   43   -none- logical
# allOK          1   -none- logical
gsg$allOK
# [1] FALSE

# Need to filter out the outliers
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}
# Removing genes: Rv2395A, Rv2561, Rv2562, Rv3098A, Rv3098c


# Identify Outlier Samples
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


###########################################################
################### NETWORK CONSTRUCTION ##################

## Pairwise Gene Co-expression similarity (Pearson)
spt <- pickSoftThreshold(expression.data) 
spt

# Plot the 𝑅2 values as a function of the soft thresholds
# We should be maximizing the 𝑅2 value and minimizing mean connectivity.
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

# Plot mean connectivity as a function of soft thresholds
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

# You can determine the soft power threshold should be set to __ as it is the spt that retains the highest mean connectivity while reaching an 𝑅2 value above 0.80.
# Set it to 4??
softPower <- 4
adjacency <- adjacency(expression.data, power = softPower)

###########################################################
################### MODULE CONSTRUCTION ###################

# Convert the adjacency matrix into a TOM similarity matrix
TOM <- TOMsimilarity(adjacency)

# Create dissimilarity/distance matrix
TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

# Identify modules
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.
# Modules
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
# 6 716 603 533 347 326 257 228 196 188 169 162 119  74  68  33 
# 0 is the number of genes that didn't fit in any module

ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


###########################################################
############# MODULE EIGENGENE IDENTIFICATION #############
# A ME (Module Eigengene) is the standardized gene expression profile for a given module.
# An eigengene is the gene whose expression is representative of the the majority of genes expressed within a module.

MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")


###########################################################
################# EXTERNAL TRAIT MATCHING #################

allTraits <- GoodSputum60_Metadata %>% 
  filter(Week == "Week 0") %>%
  select(SampleID, Week, Outcome)

Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$SampleID)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# Traits need to be continuous
str(datTraits)
# datTraits$Type2 <- ifelse(datTraits$Type2 == "W0 sputum (cure)", 1, 
#                           ifelse(datTraits$Type2 == "W0 sputum (relapse)", 2, 
#                                  ifelse(datTraits$Type2 == "W2 sputum (cure)", 3, 4)))
# datTraits$Week <- ifelse(datTraits$Week == "Week 0", 0, 2)
datTraits$Outcome <- ifelse(datTraits$Outcome == "Relapse", 1, 0)
datTraits$Week <- NULL

###########################################################
################ MODULE TRAIT ASSOCIATIONS ################

# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

###########################################################
############### TARGET GENE IDENTIFICATION ################

# You can use the gene significance along with the genes intramodular connectivity to identify potential target genes associated with a particular trait of interest. For this analysis weight will be the clinical trait.

# Define variable weight containing the weight column of datTrait
Outcome = as.data.frame(datTraits$Outcome)
names(Outcome) = "Outcome"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, Outcome, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Outcome), sep="")
names(GSPvalue) = paste("p.GS.", names(Outcome), sep="")
head(GSPvalue)

# Add the Outcome to existing module eigengenes
MET = orderMEs(cbind(MEs, Outcome))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

par(mar=c(1,1,1,1))
module = "green"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# See the genes in the modules
green_genes <- colnames(expression.data)[mergedColors == "green"]



###########################################################
#################### TRYING NEW STUFF #####################

# Correlate genes with trait
geneTraitCor <- cor(expression.data, datTraits$Outcome, use = "p")

# Module membership (kME)
MM <- cor(expression.data, mergedMEs, use = "p")

# Example: green module
green_genes <- mergedColors == "green"
hub_candidates <- colnames(expression.data)[
  green_genes & abs(MM[, "MEgreen"]) > 0.8 & abs(geneTraitCor) > 0.2
]
hub_candidates

# Example: pink module
pink_genes <- mergedColors == "pink"
hub_candidates <- colnames(expression.data)[
  pink_genes & abs(MM[, "MEpink"]) > 0.8 & abs(geneTraitCor) > 0.2
]
hub_candidates


###########################################################
############# CORRELATE MODULES WITH RELAPSE ##############

moduleTraitCor <- cor(mergedMEs, datTraits$Outcome, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expression.data))

moduleTraitDf <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = moduleTraitCor[,1],
  Pvalue = moduleTraitPvalue[,1])
moduleTraitDf[order(-abs(moduleTraitDf$Correlation)), ]


geneTraitCor <- cor(expression.data, datTraits$Outcome, use = "pairwise.complete.obs")


green_module <- mergedColors == "green"
hub_genes <- colnames(expression.data)[
  green_module &
    abs(MM[, "MEgreen"]) > 0.8 &
    abs(geneTraitCor) > 0.3
]
# [1] "Rv0096"  "Rv0152c" "Rv0341"  "Rv0342"  "Rv0661c" "Rv0885"  "Rv0887c" "Rv1075c" "Rv1187" 
# [10] "Rv1456c" "Rv1647"  "Rv1727"  "Rv1774"  "Rv1812c" "Rv1852"  "Rv1940"  "Rv1966"  "Rv2371"
# [19] "Rv2414c" "Rv2529"  "Rv2688c" "Rv2811"  "Rv3428c" "Rv3432c" "Rv3433c" "Rv3439c" "Rv3449"
# [28] "Rv3637"  "Rv3657c" "Rv3912" 

pink_module <- mergedColors == "pink"
hub_genes <- colnames(expression.data)[
  pink_module &
    abs(MM[, "MEpink"]) > 0.8 &
    abs(geneTraitCor) > 0.3
]
# [1] "Rv0266c" "Rv0371c" "Rv0728c" "Rv1471"  "Rv1472"  "Rv1930c" "Rv2052c" "Rv2236c" "Rv2600" 
# [10] "Rv2601"  "Rv2695"  "Rv2697c" "Rv2698"  "Rv2699c" "Rv3227"  "Rv3801c" "Rv3838c"


