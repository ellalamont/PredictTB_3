# Leave one out cross validation (then logistic regression?)
# 1/14/26

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850

source("Import_data.R")

library(caret)
library(glmnet)
library(pROC)

# caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 

###########################################################
##################### ORGANIZE DATA #######################

# Keep only W2 sputum
W2SputumSamples60_pipeSummary <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 2 sputum") 
P_metadata <- W2SputumSamples60_pipeSummary %>% 
  dplyr::select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)

P_TPM <- GoodSamples60_tpmf %>% dplyr::select(all_of(W2SputumSamples60_pipeSummary$SampleID2))

# Remove genes where it's all zero
P_TPM_2 <- P_TPM[rowSums(P_TPM) != 0,] # Remove columns that are all zero so the scale works for prcomp


# Put everything in one dataframe
P_TPM_t <- P_TPM_2 %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")
my_df2 <- my_df %>% 
  dplyr::select(-SampleID2)

# Remove na
my_df2 <- na.omit(my_df2) # Didn't change anything

# Remove genes with near zero variance
NearZeroGenes <- nearZeroVar(my_df2 %>% dplyr::select(-Outcome))
my_df3 <- my_df2 %>% dplyr::select(-NearZeroGenes) # Now 4025 genes

# Outcome must be factor with "Relapse" first
my_df3$Outcome <- factor(my_df3$Outcome, levels = c("Relapse", "Cure"))


###########################################################
################ KEEP OUT A VALIDATION SET ################

# Separate training and testing sets
set.seed(23)
training_samples <- my_df3$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_data  <- my_df3[training_samples, ]
validation_data <- my_df3[-training_samples, ]

# Make sure that Relapse is the positive class for caret
train_data$Outcome <- factor(train_data$Outcome, levels = c("Relapse", "Cure"))
validation_data$Outcome <- factor(validation_data$Outcome, levels = c("Relapse", "Cure"))


###########################################################
################ LOOCV LOGISTIC REGRESSION ################
# CV will chose alpha and lambda

# Define training control
set.seed(23)
train_control <- trainControl(method = "LOOCV", # Will do leave one out cross validation
                              classProbs = T, # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
                              summaryFunction = twoClassSummary, # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity
                              savePredictions = "final" # caret will store predictions from the best tuned model only
)

# Use when setting alpha = 1 for lasso only. Not during this now, letting it pick it's fav alpha
# lassoGrid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 100)) # A range of lambda values

# Train the model (LOOCV takes a little bit to run)
set.seed(23)
loocv_model <- train(Outcome ~., 
                     data = train_data, 
                     method = "glmnet", 
                     family = "binomial", # tells glmnet its a logistic regression
                     metric = "ROC", # tells caret to optimize ROC during hyperparameter tuning
                     trControl = train_control, # Provides the strategy to use
                     preProcess = c("center", "scale"), # Scales the TPM around 0... not sure if this is necessary. When included I get warning about genes having zero variances. And it changes the coefficent values, but not the actual genes used
                     tuneLength = 50) # how many hyperparameter combinations to try (50 values of lambda)

# Hyperparameters:
# Alpha: Controls if it's lasso (1), ridge (0), or elastic net (somewhere between)
# Lambda: Controls how much coefficents are shrunk towards zero, larger lambda means fewer genes are selected
# Hyperparameter tuning picks these two values

# Check the alpha and best lambda used - This is not consistent even when I set the seed...
loocv_model$bestTune
# alpha     lambda
# 2467     1 0.04316874

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# What genes did it use?
coef_df <- coef(loocv_model$finalModel, s = loocv_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame() 
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv1892"      "Rv2627c"     "Rv2665"      "Rv2745c"     "Rv3293"   

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs so positive is higher in relapse
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
  labs(y = "Coefficient", x = "Gene", title = "W2 60Cov LOOCV logistic regression train set") + 
  theme_bw()
coefficients_plot
ggsave(coefficients_plot,
       file = paste0("LOOCV_LR_W2_60Cov_Coefficients.pdf"),
       path = "Figures/PredictiveModeling/LOOCV_LR/W2_60Cov",
       width = 6, height = 5, units = "in")

###########################################################
########## MODEL ROC/AUC AND PROBABILITY BOXPLOT ##########

# Look at all of the predictions for each sample
View(loocv_model$pred)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_loocv <- roc(
  response = loocv_model$pred$obs,
  predictor = loocv_model$pred$Relapse)
plot(roc_loocv, print.auc = TRUE, main = "W2 60Cov LOOCV LR training set")

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
       title = "W2 60Cov LOOCV LR training set",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_plot
# ggsave(roc_plot,
#        file = paste0("LOOCV_LR_W2_60Cov_ROC.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W2_60Cov",
#        width = 6, height = 4, units = "in")


# Print a boxplot of the probabilities
Probability_plot <- ggplot(loocv_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "W2 60Cov LOOCV LR training set") +
  theme_bw()
Probability_plot
# ggsave(Probability_plot,
#        file = paste0("LOOCV_LR_W2_60Cov_Probablities.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W2_60Cov",
#        width = 6, height = 4, units = "in")



###########################################################
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
valiation_probabilities <- predict(
  loocv_model,
  newdata = validation_data,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = validation_data$Outcome, predictor = valiation_probabilities$Relapse)
plot(roc_val, print.auc = TRUE, main = "LOOCV logistic regression validation set")

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
       title = "W2_60cov_LOOCV logistic regression validation set",
       subtitle = paste0("AUC=", round(auc_value,3), "; alpha=", round(alpha_value,3), "; lambda=", round(lambda_value,3))) + 
  theme_bw()
roc_Val_plot
# ggsave(roc_Val_plot,
#        file = paste0("LOOCV_LR_W2_60Cov_ROC_Validation.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W2_60Cov",
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
       title = "W2 60Cov LOOCV LR validation set") +
  theme_bw()
Probability_Val_plot
# ggsave(Probability_Val_plot,
#        file = paste0("LOOCV_LR_W2_60Cov_Probabilities_Validation.pdf"),
#        path = "Figures/PredictiveModeling/LOOCV_LR/W2_60Cov",
#        width = 6, height = 4, units = "in")





