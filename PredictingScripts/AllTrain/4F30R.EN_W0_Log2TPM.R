# K-fold cross validation (with logistic regression)
# 4-fold 30-repeats: Do 4-fold cross validation 30 different times, each time with new random fold splits
# Elastic net
# log2(TPM+1)
# Not keeping out a validation set
# 5/15/26

###****************###
# RESULTS:
# Alpha = 0
# Area under the curve: 0.622
###****************###

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850


library(caret) # caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 
library(glmnet)
library(pROC)
# detach("package:CMA", unload = TRUE)


###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_nfolds <- 4
my_nrepeats <- 30

my_folderPath <- "Figures/PredictiveModeling/AllTrain/4F30R_W0_Log2TPM"
my_figureTitle <- "4F30R LR Log2(TPM+1)"


###########################################################
##################### ORGANIZE DATA #######################

tpm_t <- GoodSamples60_tpmf %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")

model_df <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  dplyr::select(SampleID2, Outcome) %>%
  mutate(Outcome = factor(Outcome, levels = c("Cure", "Relapse"))) %>%
  inner_join(tpm_t, by = "SampleID2") %>%
  column_to_rownames("SampleID2") %>%
  na.omit()

# Log2 transform
model_df_log2 <- model_df %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1)))

myX <- model_df_log2 %>% dplyr::select(-Outcome)
myY <- model_df_log2$Outcome


###########################################################
########### K-FOLD LOGISTIC REGRESSION FUNCTION ###########
# https://rpubs.com/uky994/1328850

run_glmnet_cv <- function(X, Y, nfolds = 5, nrepeats = 30) {
  
  ## Create folds ##
  set.seed(42)
  folds <- createMultiFolds(Y, k = nfolds, times = nrepeats)
  
  ## Train control ##
  train_control <- trainControl(
    method = "cv",
    index = folds,
    classProbs = T, # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
    summaryFunction = twoClassSummary, # # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity
    savePredictions = "final") # # caret will store predictions from the best tuned model only
  
  ## Tuning Grid ##
  lassoGrid <- expand.grid(
    alpha = seq(0, 1, by = 0.25),
    # alpha = 1, # When I want it to be lasso only
    lambda = 10^seq(-3, 1, length = 10))
  
  ## Train model ##
  set.seed(23)
  model <- train(
    x = X, y = Y,
    method = "glmnet",
    metric = "ROC",
    tuneGrid = lassoGrid,
    trControl = train_control,
    preProcess = "nzv", # Remove near zero genes
    standardize = TRUE)
  
  return(model)
}


###########################################################
####################### RUN THE MODEL #####################

kfold_model <- run_glmnet_cv(myX, myY)

# Check the alpha and best lambda used
kfold_model$bestTune$alpha
# [1]  0
kfold_model$bestTune$lambda
# [1] 10

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# It's all of them

# Plot the coefficients (it's backwards so negative is actually higher in relapse)
# coefficients_plot <- coef_df %>% 
#   slice(-1) %>% # Remove the y-intercept %>%
#   mutate(V2 = V1*-1) %>% # Switch the signs
#   ggplot(aes(x = reorder(rownames(.), V2), y = V2)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(y = "Coefficient", x = "Gene", title = my_figureTitle) + 
#   theme_bw()
# coefficients_plot
# ggsave(coefficients_plot,
#        file = paste0("Coefficients.pdf"),
#        path = my_folderPath,
#        width = 6, height = 5, units = "in")


###########################################################
##################### MODEL ROC/AUC #######################

# Look at all of the predictions for each sample
# View(kfold_model$pred)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse,
                 levels = c("Cure", "Relapse"),  # (Control, case)
                 direction = "<")
plot(roc_kfold, print.auc = TRUE, main = paste0(my_figureTitle, " Train set"))

auc(roc_kfold)
# Area under the curve: 0.6223

roc_mean <- mean(kfold_model$resample$ROC)
roc_mean
# 0.502381
roc_sd <- sd(kfold_model$resample$ROC)
roc_sd
# 0.2249632

# These are plotting all the repeats, just plot the average for each sample
kfold_avg <- kfold_model$pred %>%
  group_by(rowIndex) %>%
  summarize(
    obs = first(obs),
    mean_Relapse = mean(Relapse))

# Print the average ROC
roc_kfold_avg <- roc(response = kfold_avg$obs,
                     predictor = kfold_avg$mean_Relapse,
                     levels = c("Cure", "Relapse"),  # control first, case second
                     direction = "<")
plot(roc_kfold_avg, print.auc = TRUE, main = "Average 4-fold 30-repeat CV logistic regression")

