# Leave one out cross validation (then logistic regression?)
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
################ LOOCV LOGISTIC REGRESSION ################
# CV will chose alpha and lambda

# Define training control
train_control <- trainControl(method = "LOOCV", # Will do leave one out cross validation
                              classProbs = T, # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
                              summaryFunction = twoClassSummary, # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity
                              savePredictions = "final" # caret will store predictions from the best tuned model only
                              )

# Use when setting alpha = 1 for lasso only
# lassoGrid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 100)) # A range of lambda values

# Train the model (LOOCV takes a little bit to run)
set.seed(23)
loocv_model <- train(Outcome ~., 
                     data = train_data, 
                     method = "glmnet", 
                     family = "binomial", # tells glmnet its a logistic regression
                     metric = "ROC", # tells caret to optimize ROC during hyperparameter tuning
                     trControl = train_control, # Provides the strategy to use
                     tuneLength = 50) # how many hyperparameter combinations to try (50 values of lambda)

# Hyperparameters:
# Alpha: Controls if it's lasso (1), ridge (0), or elastic net (somewhere between)
# Lambda: Controls how much coefficents are shrunk towards zero, larger lambda means fewer genes are selected
# Hyperparameter tuning picks these two values

# Check the alpha and best lambda used
loocv_model$bestTune
#         alpha     lambda
# 2327 0.944898 0.06921782
# So this isn't quite lasso!

# What genes did it use?
coef_df <- coef(loocv_model$finalModel, s = lasso_loocv$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame() 
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
rownames(coef_df)
# [1] "(Intercept)" "Rv0625c"     "Rv1031"      "Rv1249c"     "Rv1357c"     "Rv1511"      "Rv1765A"    
# [8] "Rv1906c"     "Rv2106"      "Rv2279"      "Rv2355"      "Rv2511"      "Rv2827c"     "Rv2986c"    
# [15] "Rv3019c"     "Rv3187"      "Rv3472"  
colnames(coef_df)[1] <- "V1"

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene") + 
  theme_bw()

# Look at all of the predictions for each sample
View(loocv_model$pred)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_loocv <- roc(
  response = loocv_model$pred$obs,
  predictor = loocv_model$pred$Relapse)
plot(roc_loocv, print.auc = TRUE, main = "LOOCV logistic regression testing set")

# Print a boxplot of the probabilities
ggplot(loocv_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "LOOCV logistic regression testing set") +
  theme_bw()

# Now test on the validation set
valiation_probabilities <- predict(
  loocv_model,
  newdata = validation_data,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = validation_data$Outcome, predictor = valiation_probabilities$Relapse)
plot(roc_val, print.auc = TRUE, main = "LOOCV logistic regression validation set")

# Boxplot for validation set
validation_plot_df <- data.frame(
  True_Outcome = validation_data$Outcome,
  Relapse_probibility = valiation_probabilities$Relapse)
ggplot(validation_plot_df, aes(x = True_Outcome, y = Relapse_probibility)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  labs(x = "True outcome",
    y = "Predicted relapse probability",
    title = "Validation set predicted relapse probabilities") +
  theme_bw()





