# K-fold cross validation (with logistic regression) 3 fold 10 repeats
# Data has been log2 transformed and scaled (I think)
# 3/3/26

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850



library(caret)
library(glmnet)
library(pROC)

# caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 

###########################################################
##################### ORGANIZE DATA #######################

tpm_t <- GoodSamples60_tpmf %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")

model_df <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  dplyr::select(SampleID2, Outcome) %>%
  mutate(Outcome = factor(Outcome, levels = c("Relapse", "Cure"))) %>%
  inner_join(tpm_t, by = "SampleID2") %>%
  column_to_rownames("SampleID2") %>%
  na.omit()

# Log2 transform
model_df_log2 <- model_df %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1)))

###########################################################
################# SEPARATE TRAIN AND TEST #################
  
# Separate training and testing sets
set.seed(23)
trainIndexes <- model_df_log2$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_df_log2  <- model_df_log2[trainIndexes, ]
test_df_log2 <- model_df_log2[-trainIndexes, ]

# Separate X and Y
myX.train_log2 <- train_df_log2 %>% dplyr::select(-Outcome)
myY.train_log2 <- train_df_log2$Outcome

myX.test_log2  <- test_df_log2 %>% dplyr::select(-Outcome)
myY.test_log2  <- test_df_log2$Outcome

###########################################################
######################### SCALE ###########################

# Scaling needs to happen on the train set when separated, then use the same parameters on the test set, to prevent "leakage"

# But glm scales automatically so don't need to scale here?

# Remove near-zero variance genes
nzv_Indexes <- caret::nearZeroVar(myX.train_log2) # Based only on the train set! 
myX.train_log2_nzv <- myX.train_log2[, -nzv_Indexes]
myX.test_log2_nzv  <- myX.test_log2[, -nzv_Indexes]

# Scale the data
myX.train_scaled <- scale(myX.train_log2_nzv) # Scale the train data first!
# Scale the test data based on the attributes of the scaled train data
myX.test_scaled <- scale(
  myX.test_log2_nzv,
  center = attr(myX.train_scaled, "scaled:center"), # column means
  scale = attr(myX.train_scaled, "scaled:scale")) # column SDs




###########################################################
############### K-FOLD LOGISTIC REGRESSION ################
# https://rpubs.com/uky994/1328850

set.seed(23)

# 3-fold cross-validation, repeated 10 times
train_control <- trainControl(
  method = "repeatedcv", 
  number = 3, # Number of folds,
  repeats = 10, # Number of repeats
  verboseIter = F, 
  savePredictions = "final", # caret will store predictions from the best tuned model only
  classProbs = T,  # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
  summaryFunction = twoClassSummary) # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity


# For tuning alpha and lambda
lassoGrid <- expand.grid(
  # alpha = seq(0, 1, by = 0.1),
  alpha = 1, # When I want it to be lasso only
  lambda = 10^seq(-4, 1, length = 100))

# Train the model
set.seed(23)
kfold_model <- train(x = myX.train_log2_nzv, y = myY.train_log2,
                     method = "glmnet",
                     metric = "ROC", # tells caret to optimize ROC during hyperparameter tuning
                     # tuneGrid = lassoGrid, Only if doing lasso
                     trControl = train_control,
                     standardize = T, # This is in glmnet, telling it to scale (I think)
                     # preProcess = c("center", "scale"), # Not really sure if I need this for scaling or not.... If this is selected, the standardize has to be F or it won't work.
                     # tuneLength = 50,
                     tuneGrid = lassoGrid) # , # Only if doing lasso)

# Warning message:
# In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
# one multinomial or binomial class has fewer than 8  observations; dangerous ground

# Check the alpha and best lambda used
kfold_model$bestTune
# alpha      lambda
# 35     1 0.005214008


###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0440"      "Rv0614"      "Rv0733"      "Rv0887c"     "Rv1098c"    
# [7] "Rv1264"      "Rv1430"      "Rv1437"      "Rv1457c"     "Rv1719"      "Rv2030c"    
# [13] "Rv2356c"     "Rv2413c"     "Rv2463"      "Rv2465c"     "Rv2625c"     "Rv2649"     
# [19] "Rv2730"      "Rv2892c"     "Rv2986c"     "Rv3326"      "Rv3569c"     "Rv3781"     
# [25] "Rv3786c"    
# Where is Rv0625c??



# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coefficients_plot <- coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene", title = "3-fold 10-repeat CV logistic regression test set") + 
  theme_bw()
coefficients_plot
# ggsave(coefficients_plot,
#        file = paste0("3fold10repeatCV_LR_W0_60Cov_Coefficients.pdf"),
#        path = "Figures/PredictiveModeling/kfold_LR/W0_60Cov/3F10R",
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




