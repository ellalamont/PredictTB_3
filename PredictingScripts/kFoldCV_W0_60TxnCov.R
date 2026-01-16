# K-fold cross validation (with logistic regression)
# 1/14/26

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850



library(caret)
library(glmnet)
library(pROC)

# caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 

###########################################################
##################### ORGANIZE DATA #######################

# Keep only W0 sputum
W0SputumSamples60_pipeSummary <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(Outcome != "Failure")

P_metadata <- W0SputumSamples60_pipeSummary %>% select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)

P_TPM <- GoodSamples60_tpmf %>% select(all_of(W0SputumSamples60_pipeSummary$SampleID2))

# Put everything in one dataframe
P_TPM_t <- P_TPM %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")
my_df2 <- my_df %>% select(-SampleID2)

# Remove na
my_df2 <- na.omit(my_df2) # Didn't change anything

# Outcome must be factor with "Relapse" first
my_df2$Outcome <- factor(my_df2$Outcome, levels = c("Relapse", "Cure"))


###########################################################
################ KEEP OUT A VALIDATION SET ################

# Separate training and testing sets
set.seed(23)
training_samples <- my_df2$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_data  <- my_df2[training_samples, ]
validation_data <- my_df2[-training_samples, ]

# Make sure that Relapse is the positive class for caret
train_data$Outcome <- factor(train_data$Outcome, levels = c("Relapse", "Cure"))
validation_data$Outcome <- factor(validation_data$Outcome, levels = c("Relapse", "Cure"))

###########################################################
############### K-FOLD LOGISTIC REGRESSION ################
# https://rpubs.com/uky994/1328850
set.seed(23)

# 10-fold cross-validation, repeated 3 times
train_control <- trainControl(method = "repeatedcv", 
                              number = 10, # Number of folds, not to high so there is a chance of relapse in each fold (Also Tried 5)
                              repeats = 3, # Just picked something, may be too high (Also Tried 10)
                              verboseIter = F, 
                              savePredictions = "final", # caret will store predictions from the best tuned model only
                              classProbs = T,  # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
                              summaryFunction = twoClassSummary) # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity

# Use when setting alpha = 1 for lasso only. Not during this now, letting it pick it's fav alpha
# lassoGrid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 100)) # A range of lambda values

# Train the model
set.seed(23)
kfold_model <- train(Outcome ~ .,
                     data = train_data,
                     method = "glmnet",
                     metric = "ROC", # tells caret to optimize ROC during hyperparameter tuning
                     # tuneGrid = lassoGrid, Only if doing lasso
                     trControl = train_control,
                     preProcess = c("center", "scale"),
                     tuneLength = 50)

# Warning message:
# In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
# one multinomial or binomial class has fewer than 8  observations; dangerous ground
# In nominalTrainWorkflow(x = x, y = y, wts = weights, info = trainInfo,  :
# There were missing values in resampled performance measures.

# Check the alpha and best lambda used
kfold_model$bestTune
# alpha    lambda
# 2385 0.9632653 0.1425422


###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0625c"     "Rv1031"      "Rv1249c"     "Rv2106"      "Rv2279"      "Rv2355"  
# [8] "Rv2511"      "Rv2986c"     "Rv3472"    

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coefficients_plot <- coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene", title = "10-fold 3-repeat CV logistic regression test set") + 
  theme_bw()
# ggsave(coefficients_plot,
#        file = paste0("10fold3repeatCV_LR_W0_60Cov_Coefficients.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR",
#        width = 6, height = 5, units = "in")

###########################################################
########## MODEL ROC/AUC AND PROBABILITY BOXPLOT ##########

# Look at all of the predictions for each sample
# View(kfold_model$pred)
# Doesn't look like it did a very good job...

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse)
plot(roc_kfold, print.auc = TRUE, main = "10-fold 3-repeat CV logistic regression training set")

# Print a boxplot of the probabilities
ggplot(kfold_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "k-fold logistic regression training set") +
  theme_bw()
  
# These are plotting all the repeats, just plot the average for each sample
kfold_avg <- kfold_model$pred %>%
  group_by(rowIndex) %>%
  summarize(
    obs = first(obs),
    mean_Relapse = mean(Relapse)
  )
  
# Print the average ROC
roc_kfold_avg <- roc(response = kfold_avg$obs,
                     predictor = kfold_avg$mean_Relapse)
plot(roc_kfold_avg, print.auc = TRUE, main = "Average 10-fold 3-repeat CV logistic regression training set")

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
       title = "Average 10-fold 3-repeat CV logistic regression train set",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_plot
# ggsave(roc_plot,
#        file = paste0("10fold3repeatCV_LR_W0_60Cov_ROC.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR",
#        width = 6, height = 4, units = "in")

# Print the average boxplot of the probabilities
Probability_plot <- ggplot(kfold_avg, aes(x = obs, y = mean_Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "Average 10-fold 3-repeat CV logistic regression train set") +
  theme_bw()
Probability_plot
# ggsave(Probability_plot,
#        file = paste0("10fold3repeatCV_LR_W0_60Cov_Probablities.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR",
#        width = 6, height = 4, units = "in")


###########################################################
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
valiation_probabilities <- predict(
  kfold_model,
  newdata = validation_data,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = validation_data$Outcome, predictor = valiation_probabilities$Relapse)
plot(roc_val, print.auc = TRUE, main = "10-fold 3-repeat CV logistic regression validation set")

# Make the ROC a ggplot object:
roc_val_df <- data.frame(Sensitivity = roc_val$sensitivities, 
                         Specificity = roc_val$specificities)
auc_value <- auc(roc_val)
roc_Val_plot <- roc_val_df %>% 
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # scale_x_reverse() +
  coord_equal() +
  labs(x = "1-Specificity", y = "Sensitivity",
       title = "10-fold 3-repeat CV logistic regression validation set",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_Val_plot
# ggsave(roc_Val_plot,
#        file = paste0("10fold3repeatCV_LR_W0_60Cov_ROC_Validation.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR",
#        width = 6, height = 4, units = "in")

# Boxplot for validation set
validation_plot_df <- data.frame(
  True_Outcome = validation_data$Outcome,
  Relapse_probibility = valiation_probabilities$Relapse)
Probability_Val_plot <- ggplot(validation_plot_df, aes(x = True_Outcome, y = Relapse_probibility)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(x = "True outcome",
       y = "Predicted relapse probability",
       title = "10-fold 3-repeat CV logistic regression validation set") +
  theme_bw()
Probability_Val_plot
# ggsave(Probability_Val_plot,
#        file = paste0("10fold3repeatCV_LR_W0_60Cov_Probabilities_Validation.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR",
#        width = 6, height = 4, units = "in")




