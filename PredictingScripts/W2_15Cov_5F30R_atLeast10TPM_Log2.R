# W2
# K-fold cross validation (with logistic regression)
# 5-fold 30-repeats: Do 5-fold cross validation 30 different times, each time with new random fold splits
# log2(TPM+1)
# 15% Txn Cov
# Only Keep top genes that are at TPM >= 10 in >70% of samples!
# 4/7/26

library(caret)
library(glmnet)
library(pROC)
detach("package:CMA", unload = TRUE)
detach("package:DESeq2", unload = T)

###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_nfolds <- 5
my_nrepeats <- 30

my_folderPath <- "Figures/PredictiveModeling/kfold_LR/W2_15Cov/5F30R_TopGenes_log2"
my_figureTitle <- "W2 15Cov 5F30R LR Top Genes Log2"

###########################################################
##################### ORGANIZE DATA #######################

# Keep only W2 sputum (15% Txn Cov)
W2SputumSamples15_pipeSummary <- my_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>%
  filter(N_Genomic > 100000) %>% # Need to filter on something!
  filter(Txn_Coverage_f >= 15) %>%
  filter(Outcome != "Failure") %>%
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome))

W2SputumSamples15_tpmf <- All_tpm_f %>% 
  select(all_of(W2SputumSamples15_pipeSummary$SampleID2))

# Keep only genes that are TPM >= 10 in >70% of samples
keep <- rowMeans(W2SputumSamples15_tpmf >= 10) >= 0.7
W2SputumSamples15_tpmf <- W2SputumSamples15_tpmf[keep, ] # Now only 815 genes

tpm_t <- W2SputumSamples15_tpmf %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")

model_df <- W2SputumSamples15_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>%
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
##################### CREATE SEED LIST ####################
# Because I need to set seed inside the training function cause I keep getting different results

# For tuning alpha and lambda
lassoGrid <- expand.grid(
  alpha = seq(0, 1, by = 0.1), 
  # alpha = 1, # When I want it to be lasso only
  lambda = 10^seq(-4, 1, length = 100))

set.seed(23)
num_resamples <- my_nfolds * my_nrepeats
num_models <- nrow(lassoGrid)
seeds <- vector(mode = "list", length = num_resamples + 1)
for (i in 1:num_resamples) {
  seeds[[i]] <- sample.int(1000000, num_models)
}

# One final seed for final model
seeds[[num_resamples + 1]] <- sample.int(1000000, 1)


###########################################################
############### K-FOLD LOGISTIC REGRESSION ################
# https://rpubs.com/uky994/1328850

# Generate folds with a specific seed so it is repeatable
set.seed(42)
folds <- createMultiFolds(myY.train_log2, k = my_nfolds, times = my_nrepeats)

train_control <- trainControl(
  method = "cv", # "cv" instead of "repeatedcv" because are pre-defining folds for reproducibility
  index = folds, # pre-defined for reproducibility
  verboseIter = F, 
  savePredictions = "final", # caret will store predictions from the best tuned model only
  classProbs = T,  # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
  summaryFunction = twoClassSummary) # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity

# Train the model
set.seed(23)
kfold_model <- train(x = myX.train_log2_nzv, y = myY.train_log2,
                     method = "glmnet",
                     metric = "ROC",
                     tuneGrid = lassoGrid, 
                     standardize = T, # This is in glmnet, telling it to scale (I think)
                     trControl = train_control,
                     seeds = seeds) 

# Check the alpha and best lambda used
kfold_model$bestTune$alpha
# [1]  0
kfold_model$bestTune$lambda
# [1] 3.125716

kfold_model$levels

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# Uses all of them

# write.csv(coef_df, paste0(my_folderPath, "/coef_df.csv"))

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coefficients_plot <- coef_df %>% 
  dplyr::slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene", title = my_figureTitle) + 
  theme_bw()
coefficients_plot
# ggsave(coefficients_plot,
#        file = paste0("Coefficients.pdf"),
#        path = my_folderPath,
#        width = 6, height = 5, units = "in")

###########################################################
########## MODEL ROC/AUC AND PROBABILITY BOXPLOT ##########

# Look at all of the predictions for each sample
# View(kfold_model$pred)
# Doesn't look like it did a very good job...

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse,
                 levels = c("Cure", "Relapse"),  # control first, case second
                 direction = "<")
plot(roc_kfold, print.auc = TRUE, main = paste0(my_figureTitle, " Train set"))

# Print a boxplot of the probabilities
ggplot(kfold_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = paste0(my_figureTitle, " Train set")) +
  theme_bw()

roc_mean <- mean(kfold_model$resample$ROC)
roc_mean
# 0.7057407
roc_sd <- sd(kfold_model$resample$ROC)
roc_sd
# 0.2324968

# These are plotting all the repeats, just plot the average for each sample
kfold_avg <- kfold_model$pred %>%
  dplyr::filter(alpha == kfold_model$bestTune$alpha,
                lambda == kfold_model$bestTune$lambda) %>% 
  group_by(rowIndex) %>%
  summarize(
    obs = dplyr::first(obs),
    mean_Relapse = base::mean(Relapse)
  )

# Print the average ROC
roc_kfold_avg <- roc(response = kfold_avg$obs,
                     predictor = kfold_avg$mean_Relapse,
                     levels = c("Cure", "Relapse"),  # control first, case second
                     direction = "<")
plot(roc_kfold_avg, print.auc = TRUE, main = "Average 3-fold 10-repeat CV logistic regression training set")

# Make the ROC a ggplot object:
roc_kfold_avg_df <- data.frame(Sensitivity = roc_kfold_avg$sensitivities,
                               Specificity = roc_kfold_avg$specificities)
auc_value <- auc(roc_kfold_avg)
alpha_value <- kfold_model$bestTune$alpha
lambda_value <- kfold_model$bestTune$lambda
roc_plot <- roc_kfold_avg_df %>% 
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # scale_x_reverse() +
  coord_equal() +
  labs(x = "1-Specificity", y = "Sensitivity",
       title = paste0(my_figureTitle, " Train set"),
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_plot
# ggsave(roc_plot,
#        file = paste0("ROC_TrainSet.pdf"),
#        path = my_folderPath,
#        width = 6, height = 4, units = "in")

# Print the average boxplot of the probabilities
Probability_plot <- ggplot(kfold_avg, aes(x = obs, y = mean_Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = paste0(my_figureTitle, " Train set")) +
  theme_bw()
Probability_plot
# ggsave(Probability_plot,
#        file = paste0("Probablities_TrainSet.pdf"),
#        path = my_folderPath,
#        width = 6, height = 4, units = "in")

###########################################################
################ BEST PROBABILITY THRESHOLD ###############

# 3/9/26: Need to do this on the training samples!
# From roc_kfold_avg

# Choose the best specificity and sensitivity
bestThreshold <- coords(roc_kfold_avg, x = "best", best.method = "youden", transpose = F)["threshold"]
bestThreshold
# 0.2614351

# Look at all the options
roc_table <- coords(roc_kfold_avg, x = "all", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
roc_table

# # Choose 1 sensitivity but only 0.5 specificity
# bestThreshold <- roc_table %>% 
#   filter(sensitivity == 1) %>%
#   slice_min(abs(specificity - 0.5)) %>% # Choose the values closest to specificity = 0.5 
#   slice_max(threshold) %>%
#   pull(threshold)
# bestThreshold
# # Inf

###########################################################
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
validation_probabilities <- predict(
  kfold_model,
  newdata = myX.test_log2_nzv,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = myY.test_log2, predictor = validation_probabilities$Relapse,
               levels = c("Cure", "Relapse"),  # control first, case second
               direction = "<")
plot(roc_val, print.auc = TRUE, main = paste0(my_figureTitle, " Test set"))

# Make the ROC a ggplot object:
roc_val_df <- data.frame(Sensitivity = roc_val$sensitivities, 
                         Specificity = roc_val$specificities,
                         threshold   = roc_val$thresholds)
auc_value <- auc(roc_val)
roc_Val_plot <- roc_val_df %>% 
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # scale_x_reverse() +
  coord_equal() +
  labs(x = "1-Specificity", y = "Sensitivity",
       title = paste0(my_figureTitle, " Test set"),
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_Val_plot
# ggsave(roc_Val_plot,
#        file = paste0("ROC_TestSet.pdf"),
#        path = my_folderPath,
#        width = 6, height = 4, units = "in")

# Boxplot for validation set
validation_plot_df <- data.frame(
  True_Outcome = myY.test_log2,
  Relapse_probibility = validation_probabilities$Relapse)
Probability_Val_plot <- ggplot(validation_plot_df, aes(x = True_Outcome, y = Relapse_probibility)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(x = "True outcome",
       y = "Predicted relapse probability",
       title = paste0(my_figureTitle, " Test set")) +
  theme_bw()
Probability_Val_plot
# ggsave(Probability_Val_plot,
#        file = paste0("Probabilites_TestSet.pdf"),
#        path = my_folderPath,
#        width = 6, height = 4, units = "in")


# Visualize best threshold
threshold_plot <- roc_val_df %>%
  ggplot(aes(x = threshold)) + 
  geom_line(aes(y = sensitivity, color = "Sensitivity"), size = 1) +
  geom_line(aes(y = specificity, color = "Specificity"), size = 1) +
  geom_vline(xintercept = bestThreshold, linetype = "dashed") +
  labs(x = "Probability Threshold", y = "Performance", color = "", title = paste0(my_figureTitle, " Threshold Performance")) +
  theme_classic()
threshold_plot
# ggsave(threshold_plot,
#        file = paste0("ProbabilityThreshold_TestSet.pdf"),
#        path = my_folderPath,
#        width = 6, height = 4, units = "in")

# Confusion matrix with best threshold
pred_class <- ifelse(validation_probabilities$Relapse >= bestThreshold[[1]], "Relapse", "Cure")
pred_class <- factor(pred_class, levels = c("Relapse", "Cure"))
caret::confusionMatrix(data = pred_class, reference = myY.test_log2, positive = "Relapse")
# Confusion Matrix and Statistics
# 
# Reference
# Prediction Relapse Cure
# Relapse       2    3
# Cure          3    4
# 
# Accuracy : 0.5             
# 95% CI : (0.2109, 0.7891)
# No Information Rate : 0.5833          
# P-Value [Acc > NIR] : 0.8107          
# 
# Kappa : -0.0286         
# 
# Mcnemar's Test P-Value : 1.0000          
#                                           
#             Sensitivity : 0.4000          
#             Specificity : 0.5714          
#          Pos Pred Value : 0.4000          
#          Neg Pred Value : 0.5714          
#              Prevalence : 0.4167          
#          Detection Rate : 0.1667          
#    Detection Prevalence : 0.4167          
#       Balanced Accuracy : 0.4857          
#                                           
#        'Positive' Class : Relapse  