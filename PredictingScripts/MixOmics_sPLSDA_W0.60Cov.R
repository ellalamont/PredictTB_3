# MixOmics sparse PLS-DA W0 samples only
# E. Lamont
# 2/25/26

# Lets do this properly this time with the testing and training sets

# https://mixomics.org/case-studies/splsda-srbct-case-study-2/
# https://mixomics.org/methods/spls-da/



###########################################################
####################### PROCESS DATA ######################

myX <- GoodSamples60_tpmf %>% # Using TPM to start...
  # dplyr::select(-contains("THP1")) %>%
  # dplyr::select(-contains("Broth")) %>% # Lets take out broth as well
  dplyr::select(contains("W0")) %>% 
  t()

# Remove genes (columns) that are all zero
myX <- myX[, colSums(myX) != 0] # Remove columns that are all zero so the scale works for prcomp

# Want this to be Type2 I guess 
my_pipeSummary2 <- my_pipeSummary %>%
  slice(match(rownames(myX), SampleID2)) # Match the rows then just keep the one that are the same
stopifnot(all(rownames(myX) %in% my_pipeSummary2$SampleID2))
myY <- my_pipeSummary2 %>%
  pull(Type2) %>%
  factor()

summary(myY)
# W0 sputum (cure) W0 sputum (relapse) 
# 32                  11 

###########################################################
############ SEPARATE TESTING AND TRAINING SETS ###########

set.seed(23)
train <- sample(1:nrow(myX), 30) # randomly select 30 samples in training
test <- setdiff(1:nrow(myX), train) # rest is part of the test set

# store matrices into training and test set:
myX.train <- myX[train, ]
myX.test <- myX[test,]
myY.train <- myY[train] # 7 relapses
myY.test <- myY[test] # 4 relapses, 9 cures


###########################################################
############# INITIAL sPLS-DA ON TRAINING SET #############

# Initial model
train.splsda <- splsda(myX.train, myY.train, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_1 <- plotIndiv(train.splsda , comp = 1:2, 
                           group = myY.train, ind.names = F,  # colour points by class
                           ellipse = TRUE, # include 95% confidence ellipse for each class
                           legend = TRUE, 
                           title = 'W0 training set PLSDA with confidence ellipses')
train_PLSDA_1$graph
ggsave(train_PLSDA_1$graph,
       file = paste0("PLSDA1_W0.60Cov_TrainSet_v1.pdf"),
       path = "Figures/PredictiveModeling/Mixomics_PLSDA/W0_60Cov",
       width = 8, height = 5, units = "in")

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(train.splsda, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_2 <- plotIndiv(train.splsda, comp = 1:2,
                     group = myY.train, ind.names = FALSE, # colour points by class
                     background = background, # include prediction background for each class
                     legend = TRUE, title = "W0 training set PLSDA with prediction background")
ggsave(train_PLSDA_2$graph,
       file = paste0("PLSDA2_W0.60Cov_TrainSet_v1.pdf"),
       path = "Figures/PredictiveModeling/Mixomics_PLSDA/W0_60Cov",
       width = 8, height = 5, units = "in")


###########################################################
############### TUNE sPLS-DA ON TRAINING SET ##############

## ncomp: number of components ##

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.train <- perf(train.splsda, validation = "Mfold", 
                          folds = 3, nrepeat = 50, # use repeated cross-validation. folds = 3 or 5...
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.train, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")

# BER = Balanced error rate. It's really high! I think this is because it is predicting the cures well but not predicting the relapses so well?

perf.splsda.train$choice.ncomp # what is the optimal value of components according to perf()
# max.dist centroids.dist mahalanobis.dist
# overall        8              4                4
# BER            8              4                4


## keepX: number of variables ##

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.train <- tune.splsda(myX.train, myY.train, ncomp = 8, # calculate for first 8 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX)

plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning

## NEXT: FIGURE OUT WHAT NCOMP CHOICE TO USE





