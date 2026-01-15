# K-fold cross validation (with logistic regression)
# TRY WITH ALL DATA (NOT JUST 60% TXN COVERAGE)
# 1/54/26

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850



library(caret)
library(glmnet)
library(pROC)

# caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 

###########################################################
##################### ORGANIZE DATA #######################

# Keep only W0 sputum (ALL DATA)
W0SputumSamplesAll_pipeSummary <- my_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(N_Genomic > 100000) %>% # Need to filter on something!
  filter(Outcome != "Failure") %>%
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome))
  

P_metadata <- W0SputumSamplesAll_pipeSummary %>% select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)

P_TPM <- All_tpm_f %>% select(all_of(W0SputumSamplesAll_pipeSummary$SampleID2))

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
  createDataPartition(p = 0.7, list = FALSE)
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
# Warning in preProcess.default(thresh = 0.95, k = 5, freqCut = 19, uniqueCut = 10,  :
#                                 These variables have zero variances: Rv2395A, Rv2561, Rv2562, Rv3098A, Rv3098c

# Check the alpha and best lambda used
kfold_model$bestTune
# alpha    lambda
# 1591 0.6693878 0.1787414


###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0244c"     "Rv0289"      "Rv0576"      "Rv0721"      "Rv0796"     
# [7] "Rv0931c"     "Rv1304"      "Rv1420"      "Rv1555"      "Rv1812c"     "Rv1930c"    
# [13] "Rv2106"      "Rv2279"      "Rv2326c"     "Rv2331A"     "Rv2355"      "Rv2501c"    
# [19] "Rv2503c"     "Rv2649"      "Rv2730"      "Rv2910c"     "Rv3326"      "Rv3333c"    
# [25] "Rv3614c"     "Rv3616c"     "Rv3842c"     
# These genes are pretty different than the 60% coverage sets

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
coef_df %>% 
  slice(-1) %>% # Remove the y-intercept %>%
  mutate(V2 = V1*-1) %>% # Switch the signs
  ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
  geom_col() + 
  coord_flip() + 
  labs(y = "Coefficient", x = "Gene") + 
  theme_bw()

###########################################################
########## MODEL ROC/AUC AND PROBABILITY BOXPLOT ##########

# Look at all of the predictions for each sample
# View(kfold_model$pred)
# Doesn't look like it did a very good job...

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse)
plot(roc_kfold, print.auc = TRUE, main = "k-fold logistic regression testing set")

# Print a boxplot of the probabilities
ggplot(kfold_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "k-fold logistic regression testing set") +
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
plot(roc_kfold_avg, print.auc = TRUE, main = "Average k-fold logistic regression testing set")

# Print the average boxplot of the probabilities
ggplot(kfold_avg, aes(x = obs, y = mean_Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = "Average k-fold logistic regression testing set") +
  theme_bw()

###########################################################
############### TEST MODEL ON VALIDATION SET ##############

# Now test on the validation set
valiation_probabilities <- predict(
  kfold_model,
  newdata = validation_data,
  type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = validation_data$Outcome, predictor = valiation_probabilities$Relapse)
plot(roc_val, print.auc = TRUE, main = "k-fold logistic regression validation set")

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






