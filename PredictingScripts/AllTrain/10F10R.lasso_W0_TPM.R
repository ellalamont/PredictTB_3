# K-fold cross validation (with logistic regression)
# 10-fold 10-repeats
# TPM
# lasso
# Not keeping out a validation set
# 5/15/26

###****************###
# RESULTS:
# Parameters: 10-FOLD, 10 REPEAT, LASSO, TPM, W0
# Alpha = 0.5
# Area under the curve: 0.8295
# View(kfold_model$pred) looks better at predicting relapse
# Genes used: "Rv0625c"     "Rv2331A" 
###****************###

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850


library(caret) # caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 
library(glmnet)
library(pROC)
# detach("package:CMA", unload = TRUE)


###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_nfolds <- 10
my_nrepeats <- 10

my_folderPath <- "Figures/PredictiveModeling/AllTrain/10F10R.lasso_W0_TPM"
my_figureTitle <- "10F10R LR TPM"


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

myX <- model_df %>% dplyr::select(-Outcome)
myY <- model_df$Outcome


###########################################################
########### K-FOLD LOGISTIC REGRESSION FUNCTION ###########
# https://rpubs.com/uky994/1328850

run_glmnet_cv <- function(X, Y, nfolds, nrepeats) {
  
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
    # alpha = seq(0, 1, by = 0.25),
    alpha = 1, # When I want it to be lasso only
    lambda = 10^seq(-3, 1, length = 10))
  
  ## Train model ##
  set.seed(23)
  model <- train(
    x = X,
    y = Y,
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

kfold_model <- run_glmnet_cv(myX, myY, nfolds = my_nfolds, nrepeats = my_nrepeats)

# Check the alpha and best lambda used
kfold_model$bestTune$alpha
# [1] 1
kfold_model$bestTune$lambda
# [1] 0.1668101

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0625c"     "Rv2331A"  

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
View(kfold_model$pred)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse,
                 levels = c("Cure", "Relapse"),  # control first, case second
                 direction = "<")
plot(roc_kfold, print.auc = TRUE, main = paste0(my_figureTitle))

roc_mean <- mean(kfold_model$resample$ROC)
roc_mean
# 0.8282275
roc_sd <- sd(kfold_model$resample$ROC)
roc_sd
# 0.1686544

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
plot(roc_kfold_avg, print.auc = TRUE, main = "Average 4-fold 30-repeat")

auc(roc_kfold_avg)
# Area under the curve: 0.8295



##########################################################################
# Do it with permutations to see what it would be randomly

###########################################################
################### PERMUTATION OF NULL ###################
# set.seed(42)
# n_perm <- 10
# 
# perm_auc <- numeric(n_perm)
# 
# for (i in 1:n_perm) {
#   cat("Permutation:", i, "\n")
#   # Shuffle labels
#   perm_Y <- sample(myY)
#   # Train model
#   perm_model <- run_glmnet_cv(myX, perm_Y)
#   # Store AUC
#   perm_auc[i] <- mean(perm_model$resample$ROC)
# }
# 
# perm_auc
# # [1] 0.5421693 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
# # That's pretty good at 0.5....