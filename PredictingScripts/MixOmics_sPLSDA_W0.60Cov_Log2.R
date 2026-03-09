# MixOmics sparse PLS-DA W0 samples only
# v2: Trying this with Log2 transformed and scaled TPM values
# E. Lamont
# 3/2/26

# Lets do this properly this time with the testing and training sets

# https://mixomics.org/case-studies/splsda-srbct-case-study-2/
# https://mixomics.org/methods/spls-da/

source("Import_data.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("mixOmics")
library(mixOmics)
library(writexl)


###########################################################
######################### SET PATHS #######################

my_nfolds <- 5
my_nrepeats <- 30


my_folderPath <- "Figures/PredictiveModeling/Mixomics_PLSDA/W0_60Cov_Log2"
my_figureTitle <- "Mixomics sPLSDA Log2(TPM+1)"


###########################################################
##################### ORGANIZE DATA #######################

tpm_t <- GoodSamples60_tpmf %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")

model_df <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  dplyr::select(SampleID2, Outcome) %>%
  mutate(Outcome = factor(Outcome, levels = c("Relapse", "Cure"))) %>% # Can't switch these or the train/test split will be different
  inner_join(tpm_t, by = "SampleID2") %>%
  column_to_rownames("SampleID2") %>%
  na.omit()

# Log2 transform
model_df_log2 <- model_df %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1)))


###########################################################
################# SEPARATE TRAIN AND TEST #################

# TPM Transformed Data
set.seed(23)
trainIndexes <- model_df$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_df  <- model_df[trainIndexes, ]
test_df <- model_df[-trainIndexes, ]
# Separate X and Y
myX.train <- train_df %>% dplyr::select(-Outcome)
myY.train <- train_df$Outcome # %>% factor(levels = c("Cure","Relapse")) # I have no idea if I need to do this or not... Switches the coefficient plot but end result is the same
myX.test  <- test_df %>% dplyr::select(-Outcome)
myY.test  <- test_df$Outcome # %>% factor(levels = c("Cure","Relapse"))

# LOG2 Transformed Data
set.seed(23)
trainIndexes <- model_df_log2$Outcome %>%
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_df_log2  <- model_df_log2[trainIndexes, ]
test_df_log2 <- model_df_log2[-trainIndexes, ]
# Separate X and Y
myX.train_log2 <- train_df_log2 %>% dplyr::select(-Outcome)
myY.train_log2 <- train_df_log2$Outcome # %>% factor(levels = c("Cure","Relapse"))
myX.test_log2  <- test_df_log2 %>% dplyr::select(-Outcome)
myY.test_log2  <- test_df_log2$Outcome # %>% factor(levels = c("Cure","Relapse"))

###########################################################
############ REMOVE NEAR ZERO VARIANCE GENES ##############

# TPM data
nzv_Indexes <- caret::nearZeroVar(myX.train) # Based only on the train set! 
if (length(nzv_Indexes) > 0) {
  myX.train_nzv <- myX.train[, -nzv_Indexes]
  myX.test_nzv  <- myX.test[, -nzv_Indexes]
} else {
  myX.train_nzv <- myX.train
  myX.test_nzv  <- myX.test
}

# TPM Log2 data
nzv_Indexes_log2 <- caret::nearZeroVar(myX.train_log2) # Based only on the train set!
# They are the same
if (length(nzv_Indexes_log2) > 0) {
  myX.train_log2_nzv <- myX.train_log2[, -nzv_Indexes_log2]
  myX.test_log2_nzv  <- myX.test_log2[, -nzv_Indexes_log2]
} else {
  myX.train_log2_nzv <- myX.train_log2
  myX.test_log2_nzv  <- myX.test_log2
}

###########################################################
############# INITIAL sPLS-DA ON TRAINING SET #############

# Initial model
set.seed(23)
train.splsda <- mixOmics::splsda(myX.train_log2_nzv, myY.train, ncomp = 10) # set ncomp to 10 for performance assessment later


###########################################################
########## PLOT INITIAL sPLS-DA ON TRAINING SET ###########

# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_1 <- plotIndiv(train.splsda , comp = 1:2, 
                           group = myY.train, ind.names = F,  # colour points by class
                           ellipse = TRUE, # include 95% confidence ellipse for each class
                           legend = TRUE, 
                           title = 'W0 training set sPLSDA with confidence ellipses')
train_PLSDA_1$graph
# ggsave(train_PLSDA_1$graph,
#        file = paste0("sPLSDA_TrainSet.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
set.seed(23)
background = background.predict(train.splsda, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
train_PLSDA_2 <- plotIndiv(train.splsda, comp = 1:2,
                           group = myY.train, ind.names = FALSE, # colour points by class
                           background = background, # include prediction background for each class
                           legend = TRUE, title = "W0 training set sPLSDA with prediction background")
# ggsave(train_PLSDA_2$graph,
#        file = paste0("sPLSDA2_TrainSet.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")


###########################################################
############### TUNE sPLS-DA ON TRAINING SET ##############

## ncomp: number of components ##

# undergo performance evaluation in order to tune the number of components to use
set.seed(23)
perf.splsda.train <- perf(train.splsda, validation = "Mfold", 
                          folds = my_nfolds, nrepeat = my_nrepeats, # use repeated cross-validation. folds = 3 or 5...
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.train, col = color.mixo(5:7), sd = TRUE,
     legend.position = "vertical")

# BER = Balanced error rate. It's really high! I think this is because it is predicting the cures well but not predicting the relapses so well?

perf.splsda.train$choice.ncomp # what is the optimal value of components according to perf()
# max.dist centroids.dist mahalanobis.dist
# overall        1              1                1
# BER            1              1                1


## keepX: number of variables ##

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
set.seed(23)
tune.splsda.train <- tune.splsda(myX.train_log2_nzv, myY.train, ncomp = 2, # calculate for first 2 components
                                 validation = 'Mfold',
                                 folds = my_nfolds, nrepeat = my_nrepeats, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX)

plot(tune.splsda.train, col = color.jet(2)) # plot output of variable number tuning
# Each colored line is a BER. SD is based on repeated cross validation folds.
# As sPLS-DA is an iterative algorithm, values represented for a given component (e.g. comp 1 to 2) include the optimal keepX value chosen for the previous component (comp 1).
# 1-6 looks the best...

tune.splsda.train$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda
# [1] 1


tune.splsda.train$choice.keepX # what are the optimal values of variables according to tune.splsda()
# comp1 comp2 
# 270    10 
# I think these are genes being used?

# I guess I'll just use what it tells me...
optimal.ncomp.train <- tune.splsda.train$choice.ncomp$ncomp 
# 1
optimal.keepX.train <- tune.splsda.train$choice.keepX[1:optimal.ncomp.train]
# comp1 
# 270 

# Doesn't work when It's 1, has to at least be two
optimal.ncomp.train <- 2
optimal.keepX.train <- tune.splsda.train$choice.keepX[1:optimal.ncomp.train]
# comp1 comp2 
# 270    10 


###########################################################
############## FINAL sPLS-DA ON TRAINING SET ##############

# form final model with optimised values for component and variable count
set.seed(23)
final.splsda.train <- mixOmics::splsda(myX.train_log2_nzv, myY.train, 
                                       ncomp = optimal.ncomp.train, 
                                       keepX = optimal.keepX.train)

###########################################################
########### PLOT FINAL sPLS-DA ON TRAINING SET ############

sPLSDA_FinalModel_Comp12 <- plotIndiv(final.splsda.train , comp = c(1,2),
                                      group = myY.train, ind.names = F, 
                                      ellipse = TRUE, # include 95% confidence ellipse
                                      legend = TRUE,
                                      title = 'W0 Final sPLSDA comp 1&2')
# ggsave(sPLSDA_FinalModel_Comp12$graph,
#        file = paste0("sPLSDA_FinalModel_TrainSet.pdf"),
#        path = "Figures/PredictiveModeling/Mixomics_PLSDA/W0_60Cov",
#        width = 8, height = 5, units = "in")

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(final.splsda.train, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
sPLSDA_FinalModel_Comp12_v2 <- plotIndiv(final.splsda.train , comp = c(1,2),
                                         group = myY.train, ind.names = F, 
                                         background = background,
                                         legend = TRUE,
                                         title = 'W0 Final sPLSDA comp 1&2')
# ggsave(sPLSDA_FinalModel_Comp12_v2$graph,
#        file = paste0("sPLSDA_FinalModel_Comp12_W0.60Cov_TrainSet_v2.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")

###############################################################
###### VARIABLE PLOTS FROM FINAL sPLS-DA ON TRAINING SET ######

# The stability of a given feature is defined as the proportion of cross validation folds (across repeats) where it was selected for to be used for a given component. Those with the highest stability are likely to be much more “important” for a given component.
# View(perf.splsda.train$features$stable) # Extract stability values
# hmmm everything is 1....

# form new perf() object which utilises the final model
set.seed(23)
perf.final.splsda.train <- perf(final.splsda.train, 
                                folds = my_nfolds, nrepeat = my_nrepeats, # use repeated cross-validation
                                validation = "Mfold", dist = "max.dist",  # use max.dist measure
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
my_list <- perf.final.splsda.train$features$stable
my_list_df <- lapply(my_list, function(x) {
  data.frame(Feature = x)
})
write_xlsx(my_list_df, paste0(my_folderPath, "/Variable_Stability_Final.xlsx"))


## Correlation Circle Plot ##
plotVar(final.splsda.train, comp = c(1,2), cex = 3) # generate correlation circle plot


###########################################################
################ PREDICTING ON THE TEST SET ###############

# The model is applied on the test set using a specific distance metric. In this case, the Mahalanobis distance was used (arbitrarily).
# Could try different distance matrixes?
# use the model on the Xtest set
set.seed(23)
predict.splsda.Mahalanobis <- predict(object = final.splsda.train, 
                                      newdata = myX.test_log2_nzv,
                                      dist = "mahalanobis.dist",
                                      type = "prob")
# str(predict.splsda.Mahalanobis)
predict.splsda.Mahalanobis$predict
scores_Comp1 <- predict.splsda.Mahalanobis$predict[, "Relapse", 1]
tapply(scores_Comp1, myY.test, mean)

# predict.splsda.max <- predict(final.splsda.train, myX.test, 
#                                                 dist = "max.dist")
# scores_Comp1 <- predict.splsda.max$predict[, "W0 sputum (relapse)", 1]
# tapply(scores_Comp1, myY.test, mean)

# Training AUC
auc.splsda = auroc(final.splsda.train, roc.comp = 1, print = FALSE) # AUROC for the first component
# ggsave(auc.splsda$graph,
#        file = paste0("ROC_sPLSDA_FinalModel_TrainingSet.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")
auc.splsda = auroc(final.splsda.train, roc.comp = 2, print = FALSE) # AUROC for all three components

# Test Set AUC
auc.splsda <- auroc(object = final.splsda.train, 
                    newdata = myX.test_log2_nzv, outcome.test = myY.test,
                    roc.comp = 1, print = FALSE) # AUROC for all three components
# Looked at using the other components, and just using 1 is best here... Can't figure out how to plot them all on one graph
# ggsave(auc.splsda$graph,
#        file = paste0("ROC_sPLSDA_FinalModel_TestSet.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")



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
set.seed(23)
pred_test <- predict(final.splsda.train, newdata = myX.test_log2_nzv, type = "variates")

predict.splsda.Mahalanobis_variates <- predict(object = final.splsda.train, 
                                               newdata = myX.test_nzv,
                                               dist = "mahalanobis.dist",
                                               type = "variates")

# Extract components 1 & 2
test_coords <- predict.splsda.Mahalanobis_variates$variates[, 1:2]
test_df <- data.frame(Dim1 = test_coords[, 1],
                      Dim2 = test_coords[, 2],
                      Group = myY.test, # actual class
                      Set = "Test")

plot_df <- rbind(train_df, test_df) %>%
  mutate(Shape = case_when(Set == "Train" & Group == "Cure" ~ 21,
                           Set == "Train" & Group == "Relapse" ~ 24,
                           Set == "Test" & Group == "Cure" ~ 21,
                           Set == "Test" & Group == "Relapse" ~ 24),
         Fill = case_when(Set == "Train" ~ NA_character_,  # hollow
                          Set == "Test" ~ Group # filled by class
         ))
# Make Fill a factor for scale_fill_manual
plot_df$Fill <- factor(plot_df$Fill, levels = c("Relapse", "Cure"))


sPLSDA_TrainTest <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Group, shape = Shape, fill = Fill)) +
  geom_point(size = 3, stroke = 1, alpha = 0.8) +
  stat_ellipse(data = subset(plot_df, Set == "Train"),
               aes(x = Dim1, y = Dim2, color = Group),
               type = "norm", level = 0.95, linetype = 2) +
  scale_shape_identity() +
  scale_fill_manual(values = c("Cure" = "#0072B2", "Relapse" = "#bc5300"), na.value = NA, guide = "none") +
  scale_color_manual(values = c("Cure" = "#0072B2", "Relapse" = "#bc5300")) +
  labs(title = "sPLS-DA Comp 1 & 2: Train (hollow) + Test (filled)") +
  theme_bw()
sPLSDA_TrainTest
# ggsave(sPLSDA_TrainTest,
#        file = paste0("PLSDA_TrainTest.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")

###########################################################
################# BOXPLOT OF PROBABILITIES ################

# Choose component (e.g., comp 1)
comp_to_use <- 1
# Extract relapse probability
relapse_prob <- predict.splsda.Mahalanobis$predict[, "Relapse", comp_to_use]

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
# ggsave(Probability_Val_plot,
#        file = paste0("Probabilities_TestSet.pdf"),
#        path = my_folderPath,
#        width = 8, height = 5, units = "in")






