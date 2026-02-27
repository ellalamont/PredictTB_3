# MixOmics sparse PLS-DA W2 samples only
# E. Lamont
# 2/25/26

# Lets do this properly this time with the testing and training sets

# https://mixomics.org/case-studies/splsda-srbct-case-study-2/
# https://mixomics.org/methods/spls-da/

source("Import_data.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("mixOmics")
library(mixOmics)

###########################################################
####################### PROCESS DATA ######################

myX <- GoodSamples60_tpmf %>% # Using TPM to start...
  # dplyr::select(-contains("THP1")) %>%
  # dplyr::select(-contains("Broth")) %>% # Lets take out broth as well
  dplyr::select(contains("W2")) %>% 
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
# W2 sputum (cure) W2 sputum (relapse) 
# 8                   5 

###########################################################
############ SEPARATE TESTING AND TRAINING SETS ###########

# set.seed(23) # 23
# train <- sample(1:nrow(myX), 30) # randomly select 30 samples in training
# test <- setdiff(1:nrow(myX), train) # rest is part of the test set
# 
# # store matrices into training and test set:
# myX.train <- myX[train, ]
# myX.test <- myX[test,]
# myY.train <- myY[train] # 7 relapses
# myY.test <- myY[test] # 4 relapses, 9 cures
# myY.test

# Make sure they are statified evenly?
set.seed(23)
relapse_idx <- which(myY == "W2 sputum (relapse)")
cure_idx    <- which(myY == "W2 sputum (cure)")
set.seed(23)
train_relapse <- sample(relapse_idx, round(0.7 * length(relapse_idx)))
train_cure    <- sample(cure_idx, round(0.7 * length(cure_idx)))

train <- c(train_relapse, train_cure)
test  <- setdiff(1:length(myY), train)

myX.train <- myX[train, ]
myX.test  <- myX[test,]

myY.train <- myY[train]
myY.test  <- myY[test]
myY.test


###########################################################
############# INITIAL sPLS-DA ON TRAINING SET #############

# Initial model
train.splsda <- splsda(myX.train, myY.train, ncomp = 2)  # Had to change this to 2 for the perf to work

###########################################################
########## PLOT INITIAL sPLS-DA ON TRAINING SET ###########

# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_1 <- plotIndiv(train.splsda , comp = 1:2, 
                           group = myY.train, ind.names = F,  # colour points by class
                           ellipse = TRUE, # include 95% confidence ellipse for each class
                           legend = TRUE, 
                           title = 'W2 training set PLSDA with confidence ellipses')
# ggsave(train_PLSDA_1$graph,
#        file = paste0("PLSDA1_W2.60Cov_TrainSet_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(train.splsda, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_2 <- plotIndiv(train.splsda, comp = 1:2,
                           group = myY.train, ind.names = FALSE, # colour points by class
                           background = background, # include prediction background for each class
                           legend = TRUE, title = "W2 training set PLSDA with prediction background")
# ggsave(train_PLSDA_2$graph,
#        file = paste0("PLSDA2_W2.60Cov_TrainSet_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")


###########################################################
############### TUNE sPLS-DA ON TRAINING SET ##############

## ncomp: number of components ##

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.train <- perf(train.splsda, validation = "loo", 
                          # folds = 3, nrepeat = 10, # use repeated cross-validation. folds = 3 or 5...
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.train, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")

# BER = Balanced error rate. It's really high! I think this is because it is predicting the cures well but not predicting the relapses so well?

perf.splsda.train$choice.ncomp # what is the optimal value of components according to perf()
# max.dist centroids.dist mahalanobis.dist
# overall        1              1                1
# BER            2              2                2


## keepX: number of variables ##

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.train <- tune.splsda(myX.train, myY.train, ncomp = 2, # calculate for first 2 components
                                 validation = 'loo',
                                 # folds = 5, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX)
# Note that the number of components cannot be reliably tuned with nrepeat < 3 or validaion = 'loo'.

plot(tune.splsda.train, col = color.jet(2)) # plot output of variable number tuning
# Note: sd bars cannot be calculated when nrepeat = 1.

tune.splsda.train$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda
# [1] NULL

tune.splsda.train$choice.keepX # what are the optimal values of variables according to tune.splsda()
# comp1 comp2 
# 4     1 
# I think these are genes being used?

# I guess I'll just use what it tells me...
optimal.ncomp.train <- tune.splsda.train$choice.ncomp$ncomp 
# NULL
optimal.keepX.train <- tune.splsda.train$choice.keepX[1:optimal.ncomp.train]
# Error in 1:optimal.ncomp.train : argument of length 0

# Doesn't work when It's 1, has to at least be two
optimal.ncomp.train <- 2
optimal.keepX.train <- tune.splsda.train$choice.keepX[1:optimal.ncomp.train]
# comp1 comp2 
# 4     1


###########################################################
############## FINAL sPLS-DA ON TRAINING SET ##############

# form final model with optimised values for component and variable count
final.splsda.train <- splsda(myX.train, myY.train, 
                             ncomp = optimal.ncomp.train, 
                             keepX = optimal.keepX.train)

###########################################################
########### PLOT FINAL sPLS-DA ON TRAINING SET ############

sPLSDA_FinalModel_Comp12 <- plotIndiv(final.splsda.train , comp = c(1,2),
                                      group = myY.train, ind.names = F, 
                                      ellipse = TRUE, # include 95% confidence ellipse
                                      legend = TRUE,
                                      title = 'W2 Final sPLSDA comp 1&2')
# ggsave(sPLSDA_FinalModel_Comp12$graph,
#        file = paste0("sPLSDA_FinalModel_Comp12_W2.60Cov_TrainSet_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(final.splsda.train, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
sPLSDA_FinalModel_Comp12_v2 <- plotIndiv(final.splsda.train , comp = c(1,2),
                                         group = myY.train, ind.names = F, 
                                         background = background,
                                         legend = TRUE,
                                         title = 'W2 Final sPLSDA comp 1&2')
# ggsave(sPLSDA_FinalModel_Comp12_v2$graph,
#        file = paste0("sPLSDA_FinalModel_Comp12_W2.60Cov_TrainSet_v2.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")

###############################################################
###### VARIABLE PLOTS FROM FINAL sPLS-DA ON TRAINING SET ######

# The stability of a given feature is defined as the proportion of cross validation folds (across repeats) where it was selected for to be used for a given component. Those with the highest stability are likely to be much more “important” for a given component.
# View(perf.splsda.train$features$stable) # Extract stability values
# hmmm everything is 1....

# form new perf() object which utilises the final model
perf.final.splsda.train <- perf(final.splsda.train, 
                                # folds = 3, nrepeat = 50, # use repeated cross-validation
                                validation = "loo", dist = "max.dist",  # use max.dist measure
                                progressBar = FALSE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
# par(mfrow=c(1,3))
plot(perf.final.splsda.train$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.final.splsda.train$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.final.splsda.train$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
# Saved with Export

# Save the variable features table
# my_list <- perf.final.splsda.train$features$stable
# my_list_df <- lapply(my_list, function(x) {
#   data.frame(Feature = x)
# })
# write_xlsx(my_list_df, "Figures/PredictiveModeling/Mixomics_PLSDA/Variable_Stability_FinalsPLSDA_W2.60Cov.xlsx")

## Correlation Circle Plot ##
plotVar(final.splsda.train, comp = c(1,2), cex = 3) # generate correlation circle plot
# Saved with export


###########################################################
################ PREDICTING ON THE TEST SET ###############

# The model is applied on the test set using a specific distance metric. In this case, the Mahalanobis distance was used (arbitrarily).
# Could try different distance matrixes?
# use the model on the Xtest set
predict.splsda.Mahalanobis <- predict(object = final.splsda.train, 
                                      newdata = myX.test,
                                      dist = "mahalanobis.dist",
                                      type = "prob")
# str(predict.splsda.Mahalanobis)
predict.splsda.Mahalanobis$predict
scores_Comp1 <- predict.splsda.Mahalanobis$predict[, "W2 sputum (relapse)", 1]
tapply(scores_Comp1, myY.test, mean)

# predict.splsda.max <- predict(final.splsda.train, myX.test, 
#                                                 dist = "max.dist")
# scores_Comp1 <- predict.splsda.max$predict[, "W2 sputum (relapse)", 1]
# tapply(scores_Comp1, myY.test, mean)

# Training AUC
auc.splsda = auroc(final.splsda.train, roc.comp = 1, print = FALSE) # AUROC for the first component
# ggsave(auc.splsda$graph,
#        file = paste0("ROC_sPLSDA_FinalModel_Comp1_W2.60Cov_TrainingSet_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")
auc.splsda = auroc(final.splsda.train, roc.comp = 2, print = FALSE) # AUROC for all three components

# Test Set AUC
auc.splsda <- auroc(object = final.splsda.train, 
                    newdata = myX.test, outcome.test = myY.test,
                    roc.comp = 1, print = FALSE) # AUROC for all three components
# Looked at using the other components, and just using 1 is best here... Can't figure out how to plot them all on one graph
# ggsave(auc.splsda$graph,
#        file = paste0("ROC_sPLSDA_FinalModel_Comp1_W2.60Cov_TestValidation_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")


###########################################################
############### TEST MODEL WITH ROC FUNCTION ##############
# Like what I was doing with the kFoldCV. This works but can only plot one comp at a time, doesn't combine them.

# # Now test on the validation set
# valiation_probabilities <- predict(
#   final.splsda.train,
#   newdata = myX.test,
#   dist = "mahalanobis.dist",
#   type = "prob") # To get probabilities
# 
# scores_Comp1 <- valiation_probabilities$predict[, "W2 sputum (relapse)", 1]
# tapply(scores_Comp1, myY.test, mean)
# # W2 sputum (cure) W2 sputum (relapse) 
# # 0.2046237           0.2930939 
# 
# # ROC for validation set
# roc_val <- roc(response = myY.test, predictor = scores_Comp1)
# # Setting levels: control = W2 sputum (cure), case = W2 sputum (relapse)
# # Setting direction: controls < cases # So higher scores are relapses 
# plot(roc_val, print.auc = TRUE, main = "MixOmics sPLSDA Comp1")
# 
# # Make the ROC a ggplot object:
# roc_val_df <- data.frame(Sensitivity = roc_val$sensitivities, 
#                          Specificity = roc_val$specificities)
# auc_value <- auc(roc_val)
# roc_Val_plot <- roc_val_df %>% 
#   ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
#   geom_path(linewidth = 1) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   # scale_x_reverse() +
#   coord_equal() +
#   labs(x = "1-Specificity", y = "Sensitivity",
#        title = "MixOmics sPLSDA Comp1",
#        subtitle = paste0("AUC=", round(auc_value,3))) + 
#   theme_bw()
# roc_Val_plot

###########################################################
############# ADD TEST POINTS TO sPLSDA GRAPH #############
# Have to extract the data by hand and make a new graph from it

# Training sample coordinates
train_coords <- final.splsda.train$variates$X[, 1:2]
train_df <- data.frame(Dim1 = train_coords[, 1],
                       Dim2 = train_coords[, 2],
                       Group = myY.train,
                       Set = "Train")

# Use predict() to get variates for test set
pred_test <- predict(final.splsda.train, newdata = myX.test, type = "variates")

predict.splsda.Mahalanobis_variates <- predict(object = final.splsda.train, 
                                               newdata = myX.test,
                                               dist = "mahalanobis.dist",
                                               type = "variates")

# Extract components 1 & 2
test_coords <- predict.splsda.Mahalanobis_variates$variates[, 1:2]
test_df <- data.frame(Dim1 = test_coords[, 1],
                      Dim2 = test_coords[, 2],
                      Group = myY.test, # actual class
                      Set = "Test")

plot_df <- rbind(train_df, test_df) %>%
  mutate(Shape = case_when(Set == "Train" & Group == "W2 sputum (cure)" ~ 21,
                           Set == "Train" & Group == "W2 sputum (relapse)" ~ 24,
                           Set == "Test" & Group == "W2 sputum (cure)" ~ 21,
                           Set == "Test" & Group == "W2 sputum (relapse)" ~ 24),
         Fill = case_when(Set == "Train" ~ NA_character_,  # hollow
                          Set == "Test" ~ Group # filled by class
         ))
# Make Fill a factor for scale_fill_manual
plot_df$Fill <- factor(plot_df$Fill, levels = c("W2 sputum (cure)", "W2 sputum (relapse)"))


sPLSDA_TrainTest <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Group, shape = Shape, fill = Fill)) +
  geom_point(size = 3, stroke = 1, alpha = 0.8) +
  stat_ellipse(data = subset(plot_df, Set == "Train"),
               aes(x = Dim1, y = Dim2, color = Group),
               type = "norm", level = 0.95, linetype = 2) +
  scale_shape_identity() +
  scale_fill_manual(values = c("W2 sputum (cure)" = "#0072B2", "W2 sputum (relapse)" = "#bc5300"), na.value = NA, guide = "none") +
  scale_color_manual(values = c("W2 sputum (cure)" = "#0072B2", "W2 sputum (relapse)" = "#bc5300")) +
  labs(title = "sPLS-DA Comp 1 & 2: Train (hollow) + Test (filled)") +
  theme_bw()
sPLSDA_TrainTest
# ggsave(sPLSDA_TrainTest,
#        file = paste0("PLSDA_W2.60Cov_TrainTest_v1.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W2_60Cov",
#        width = 8, height = 5, units = "in")

###########################################################
################# BOXPLOT OF PROBABILITIES ################

# Choose component (e.g., comp 1)
comp_to_use <- 1
# Extract relapse probability
relapse_prob <- predict.splsda.Mahalanobis$predict[, "W2 sputum (relapse)", comp_to_use]

validation_plot_df <- data.frame(
  True_Outcome = myY.test,
  Relapse_probibility = relapse_prob)
Probability_Val_plot <- ggplot(validation_plot_df, aes(x = True_Outcome, y = Relapse_probibility)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(x = "True outcome",
       y = "Predicted relapse probability",
       title = "MixOmics") +
  theme_bw()
Probability_Val_plot










