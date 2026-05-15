# Leave one out cross validation (with logistic regression)
# log2(TPM+1)
# Not keeping out a validation set
# 5/15/26

# Results: This is NOT good! Don't do LOOCV!

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850


library(caret) # caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 
library(glmnet)
library(pROC)
# detach("package:CMA", unload = TRUE)


###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_folderPath <- "Figures/PredictiveModeling/AllTrain/LOO_W0_Log2TPM"
my_figureTitle <- "LOO LR Log2(TPM+1)"


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

myX <- model_df_log2 %>% dplyr::select(-Outcome)
myY <- model_df_log2$Outcome


###########################################################
################## CREATE LOOCV FUNCTION ##################

run_glmnet_loocv <- function(X, Y) {
  
  ## Train control ##
  train_control <- trainControl(
    method = "LOOCV",
    classProbs = T, # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
    summaryFunction = twoClassSummary, # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity
    savePredictions = "final") # caret will store predictions from the best tuned model only
  
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
    preProcess = "nzv",
    standardize = TRUE)
  
  return(model)
}

###########################################################
##################### RUN THE FUNCTION ####################

loocv_model <- run_glmnet_loocv(myX, myY)

# Check the alpha and best lambda used
loocv_model$bestTune$alpha
# [1]  1 # It better be 1 because I set it to be 1 
loocv_model$bestTune$lambda
# [1] 0.05994843

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(loocv_model$finalModel, s = loocv_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
# [1] "(Intercept)" "Rv0046c"     "Rv0289"      "Rv0396"      "Rv0614"      "Rv0625c"    
# [7] "Rv1098c"     "Rv1439c"     "Rv1595"      "Rv1721c"     "Rv2031c"     "Rv2106"     
# [13] "Rv2351c"     "Rv2356c"     "Rv2475c"     "Rv2601"      "Rv2770c"     "Rv2986c"    
# [19] "Rv3219"      "Rv3569c"     "Rv3616c" 

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

# Print the ROC: Not sure if this is plotting in the right direction...
my_roc <- roc(response = loocv_model$pred$obs,
                 predictor = loocv_model$pred$Relapse,
                 levels = c("Relapse", "Cure"),  
                 direction = ">")
plot(my_roc, print.auc = TRUE, main = paste0(my_figureTitle, " Train set"))

auc(my_roc)
# Area under the curve: 0.5341

# Print a boxplot of the probabilities
ggplot(loocv_model$pred, aes(x = obs, y = Relapse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(x = "Observed outcome", 
       y = "Predicted relapse probability",
       title = paste0(my_figureTitle, " Train set")) +
  theme_bw()

