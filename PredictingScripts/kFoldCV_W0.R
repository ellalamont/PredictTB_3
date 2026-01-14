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

# Outcome must be factor with "Relapse" first (for caret convention??)
my_df2$Outcome <- factor(my_df2$Outcome, levels = c("Relapse", "Cure"))


###########################################################
###################### USE ALL DATA #######################
# https://rpubs.com/uky994/1328850
set.seed(23)

# 10-fold cross-validation, repeated 3 times
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = F, savePredictions = "final", classProbs = T, summaryFunction = twoClassSummary)

# --- Lasso Regression with caret ---
# For Lasso, we specify a tuneGrid where alpha is fixed at 1.
# caret will then search for the best lambda.

lassoGrid <- expand.grid(alpha = 1,
                         lambda = 10^seq(-3, 1, length = 100)) # A range of lambda values

cat("\n--- Training Lasso Regression with caret ---\n")

lasso_caret <- train(
  Outcome ~ .,
  data = my_df2,
  method = "glmnet",
  tuneGrid = lassoGrid,
  trControl = train_control,
  preProcess = c("center", "scale"),
  metric = "ROC"
)

# Print the model results
print(lasso_caret)
# Get the best lambda found by caret
cat("Optimal lambda for Lasso (caret):", lasso_caret$bestTune$lambda, "\n") # 0.1668101 

# Get the coefficients from the best model
coef_df <- coef(lasso_caret$finalModel, s = lasso_caret$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
rownames(coef_df)
# [1] "(Intercept)" "Rv0523c"     "Rv0625c"     "Rv1031"      "Rv1775"      "Rv2019"
# [7] "Rv2331A"     "Rv2896c" 
  
# Print the ROC: Not sure if this is plotting in the right direction...
best_lambda <- lasso_caret$bestTune$lambda
pred_best <- lasso_caret$pred %>%
  dplyr::filter(lambda == best_lambda)
roc_kfold <- roc(
  response = pred_best$obs,
  predictor = pred_best$Relapse,  # PROBABILITY
  levels = c("Cure", "Relapse"),
  direction = "<")
plot(roc_kfold, print.auc = TRUE)
  
boxplot(
  Relapse ~ obs,
  data = pred_best,
  ylab = "Predicted probability of relapse",
  xlab = "True outcome",
  main = "Out-of-fold relapse probabilities (10-fold CV)"
)
  

