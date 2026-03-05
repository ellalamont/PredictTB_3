# K-fold cross validation (with logistic regression)
# 5-fold 30-repeats: Do 5-fold cross validation 30 different times, each time with new random fold splits
# log2(TPM+1)
# Choosing specific genes that have appeared in other models
# 3/5/26



# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850


library(caret)
library(glmnet)
library(pROC)
detach("package:CMA", unload = TRUE)

# caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 


###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_nfolds <- 5
my_nrepeats <- 30

my_folderPath <- "Figures/PredictiveModeling/kfold_LR/W0_60Cov/5F30R_60Gene_Log2"
my_figureTitle <- "5F30R LR Log2 69-Gene Subset"


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
##################### SUBSET GENES ########################


GeneList <- read_excel("Data/PredictorGenes/ListofPredictorGenes_byHand.xlsx", sheet = "AllGenes") %>%
  pull() %>%
  unique()

model_df <- model_df %>% dplyr::select(Outcome, all_of(GeneList))
model_df_log2 <- model_df_log2 %>% dplyr::select(Outcome, all_of(GeneList))

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
########################### SCALING #######################
# Scaling needs to happen on the train set when separated, then use the same parameters on the test set, to prevent "leakage"

# But glm scales automatically so don't need to scale here?


# Scale the data
# myX.train_scaled <- scale(myX.train_log2_nzv) # Scale the train data first!
# Scale the test data based on the attributes of the scaled train data
# myX.test_scaled <- scale(
#   myX.test_log2_nzv,
#   center = attr(myX.train_scaled, "scaled:center"), # column means
#   scale = attr(myX.train_scaled, "scaled:scale")) # column SDs


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
folds <- createMultiFolds(myY.train, k = my_nfolds, times = my_nrepeats)

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

# Warning message:
# In lognet(x, is.sparse, y, weights, offset, alpha, nobs,  ... :
# one multinomial or binomial class has fewer than 8  observations; dangerous ground

# Check the alpha and best lambda used
kfold_model$bestTune$alpha
# [1]  0.1
kfold_model$bestTune$lambda
# [1] 0.02656088


###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)

# write.csv(coef_df, paste0(my_folderPath, "/coef_df.csv"))

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coefficients_plot <- coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
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
# 0.7013333
roc_sd <- sd(kfold_model$resample$ROC)
roc_sd
# 0.3109547

# These are plotting all the repeats, just plot the average for each sample
kfold_avg <- kfold_model$pred %>%
  group_by(rowIndex) %>%
  summarize(
    obs = first(obs),
    mean_Relapse = mean(Relapse)
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
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
validation_probabilities <- predict(
  kfold_model,
  newdata = myX.test_log2,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = myY.test_log2, predictor = validation_probabilities$Relapse,
               levels = c("Cure", "Relapse"),  # control first, case second
               direction = "<")
plot(roc_val, print.auc = TRUE, main = paste0(my_figureTitle, " Test set"))

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




