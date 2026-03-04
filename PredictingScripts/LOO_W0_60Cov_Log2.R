# Leave one out cross validation (then logistic regression?)
# Log transforming the samples
# 3/4/26

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

# TPM Transformed Data
set.seed(23)
trainIndexes <- model_df$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_df  <- model_df[trainIndexes, ]
test_df <- model_df[-trainIndexes, ]
# Separate X and Y
myX.train <- train_df %>% dplyr::select(-Outcome)
myY.train <- train_df$Outcome
myX.test  <- test_df %>% dplyr::select(-Outcome)
myY.test  <- test_df$Outcome

# LOG2 Transformed Data
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
############ REMOVE NEAR ZERO VARIANCE GENES ##############

# TPM data
nzv_Indexes <- caret::nearZeroVar(myX.train) # Based only on the train set! 
myX.train_nzv <- myX.train[, -nzv_Indexes]
myX.test_nzv  <- myX.test[, -nzv_Indexes]

# TPM Log2 data
nzv_Indexes_log2 <- caret::nearZeroVar(myX.train_log2) # Based only on the train set! 
myX.train_log2_nzv <- myX.train_log2[, -nzv_Indexes_log2]
myX.test_log2_nzv  <- myX.test_log2[, -nzv_Indexes_log2]

# They are the same


###########################################################
##################### CREATE SEED LIST ####################
# Because I need to set seed inside the training function cause I keep getting different results

# For tuning alpha and lambda
lassoGrid <- expand.grid(
  alpha = seq(0, 1, by = 0.1),
  # alpha = 1, # When I want it to be lasso only
  lambda = 10^seq(-4, 1, length = 100))

set.seed(23)
num_resamples <- nrow(myX.train_log2_nzv)
num_models <- nrow(lassoGrid)
seeds <- vector(mode = "list", length = num_resamples + 1)
for (i in 1:num_resamples) {
  seeds[[i]] <- sample.int(1000000, num_models)
}

# One final seed for final model
seeds[[num_resamples + 1]] <- sample.int(1000000, 1)



###########################################################
################ LOOCV LOGISTIC REGRESSION ################
# CV will chose alpha and lambda

# Define training control
train_control <- trainControl(
  method = "LOOCV", # Will do leave one out cross validation
  classProbs = T, # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
  summaryFunction = twoClassSummary, # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity
  savePredictions = "final", # caret will store predictions from the best tuned model only
  seeds = seeds
)

# Train the model (LOOCV takes a little bit to run)
set.seed(23)
loocv_model <- train(x = myX.train_log2_nzv, y = myY.train,
                     method = "glmnet", 
                     family = "binomial", # tells glmnet its a logistic regression
                     metric = "ROC", # tells caret to optimize ROC during hyperparameter tuning
                     trControl = train_control, # Provides the strategy to use
                     standardize = T,
                     # preProcess = c("center", "scale"), # Scales the TPM around 0... not sure if this is necessary. When included I get warning about genes having zero variances. And it changes the coefficent values, but not the actual genes used
                     tuneGrid = lassoGrid)

# Hyperparameters:
# Alpha: Controls if it's lasso (1), ridge (0), or elastic net (somewhere between)
# Lambda: Controls how much coefficents are shrunk towards zero, larger lambda means fewer genes are selected
# Hyperparameter tuning picks these two values

# Check the alpha and best lambda used
loocv_model$bestTune$alpha
# [1] 0.1
loocv_model$bestTune$lambda
# [1] 0.1519911

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# What genes did it use?
coef_df <- coef(loocv_model$finalModel, s = lasso_loocv$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame() 
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0582"      "Rv0625c"     "Rv1430"      "Rv2413c"     "Rv2688c"    
# [7] "Rv2896c"   

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene") + 
  theme_bw()
coefficients_plot <- coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene", title = "LOOCV logistic regression train set Log2") + 
  theme_bw()
coefficients_plot
# ggsave(coefficients_plot,
#        file = paste0("LOOCV_LR_W0_60Cov_Coefficients.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W0_60Cov_Log2",
#        width = 6, height = 5, units = "in")

###########################################################
########## MODEL ROC/AUC AND PROBABILITY BOXPLOT ##########

# Look at all of the predictions for each sample
View(loocv_model$pred)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_loocv <- roc(
  response = loocv_model$pred$obs,
  predictor = loocv_model$pred$Relapse,
  levels = c("Cure", "Relapse"),  # control first, case second
  direction = "<")
plot(roc_loocv, print.auc = TRUE, main = "LOOCV LR Train, log2(TPM+1)")

# Make the ROC a ggplot object:
roc_loocv_df <- data.frame(Sensitivity = roc_loocv$sensitivities,
                           Specificity = roc_loocv$specificities)
auc_value <- auc(roc_loocv)
alpha_value <- loocv_model$bestTune$alpha
lambda_value <- loocv_model$bestTune$lambda
roc_plot <- roc_loocv_df %>% 
  ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # scale_x_reverse() +
  coord_equal() +
  labs(x = "1-Specificity", y = "Sensitivity",
       title = "LOOCV LR Train, log2(TPM+1)",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_plot
# ggsave(roc_plot,
#        file = paste0("LOOCV_LR_W0_60Cov_ROC.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W0_60Cov_Log2",
#        width = 6, height = 4, units = "in")


# Print a boxplot of the probabilities
Probability_plot <- ggplot(loocv_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "LOOCV logistic regression training set") +
  theme_bw()
Probability_plot
# ggsave(Probability_plot,
#        file = paste0("LOOCVCV_LR_W0_60Cov_Probablities.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W0_60Cov_Log2",
#        width = 6, height = 4, units = "in")



###########################################################
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
valiation_probabilities <- predict(
  loocv_model,
  newdata = myX.test_log2_nzv,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = myY.test_log2, 
               predictor = valiation_probabilities$Relapse,
               levels = c("Cure", "Relapse"),  # control first, case second
               direction = "<")
plot(roc_val, print.auc = TRUE, main = "LOOCV LR Log2(TPM+1) validation set")

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
       title = "LOOCV LR Log2(TPM+1) validation set",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_Val_plot
# ggsave(roc_Val_plot,
#        file = paste0("LOOCV_LR_W0_60Cov_ROC_Validation.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W0_60Cov_Log2",
#        width = 6, height = 4, units = "in")



# Boxplot for validation set
validation_plot_df <- data.frame(
  True_Outcome = myY.test_log2,
  Relapse_probibility = valiation_probabilities$Relapse)
Probability_Val_plot <- ggplot(validation_plot_df, aes(x = True_Outcome, y = Relapse_probibility)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(x = "True outcome",
       y = "Predicted relapse probability",
       title = "LOOCV LR Log2(TPM+1) validation set") +
  theme_bw()
Probability_Val_plot
# ggsave(Probability_Val_plot,
#        file = paste0("LOOCV_LR_W0_60Cov_Probabilities_Validation.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W0_60Cov_Log2",
#        width = 6, height = 4, units = "in")





