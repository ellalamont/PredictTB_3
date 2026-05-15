# K-fold cross validation (with logistic regression)
# 5-fold 30-repeats: Do 5-fold cross validation 30 different times, each time with new random fold splits
# log2(TPM+1)
# Not keeping out a validation set
# 5/14/26

# Results: It's really bad when alpha = 1 (lasso).... Which I don't understand why NOT doing the validation set thing would make it worse....
# Results: If I don't set alpha, it makes it 0, and is also pretty bad...

# https://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
# https://rpubs.com/uky994/1328850


library(caret) # caret package (train function) handles the data splitting, resampling, model training, and performance evaluation 
library(glmnet)
library(pROC)
# detach("package:CMA", unload = TRUE)


###########################################################
############# SET FOLD AND REPEAT PARAMETERS ##############

my_nfolds <- 5
my_nrepeats <- 30

my_folderPath <- "Figures/PredictiveModeling/AllTrain/5F30R_W0_Log2TPM"
my_figureTitle <- "5F30R LR Log2(TPM+1)"


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
##################### CREATE SEED LIST ####################
# Because I need to set seed inside the training function cause I keep getting different results

# For tuning alpha and lambda
lassoGrid <- expand.grid(
  alpha = seq(0, 1, by = 0.1), 
  # alpha = 1, # When I want it to be lasso only
  lambda = 10^seq(-3, 1, length = 20))

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
folds <- createMultiFolds(myY, k = my_nfolds, times = my_nrepeats)

train_control <- trainControl(
  method = "cv", # "cv" instead of "repeatedcv" because are pre-defining folds for reproducibility
  index = folds, # pre-defined for reproducibility
  verboseIter = F, 
  savePredictions = "final", # caret will store predictions from the best tuned model only
  classProbs = T,  # Caret will produce class probabilities not just predicted classes (cure/relapse). Needed to compute ROC/AUC
  seeds = seeds,
  summaryFunction = twoClassSummary) # For binary classification. Means caret will compute ROC/AUC, sensitivity, specificity

# Train the model
set.seed(23)
kfold_model <- train(x = myX, y = myY,
                     preProcess = c("nzv"), # to remove near zero genes
                     method = "glmnet",
                     metric = "ROC",
                     tuneGrid = lassoGrid, 
                     standardize = T, # This is in glmnet, telling it to scale (I think)
                     trControl = train_control) 
# Wow no warning messages! 

# Check the alpha and best lambda used
kfold_model$bestTune$alpha
# [1]  1 # It better be 1 because I set it to be 1 
kfold_model$bestTune$lambda
# [1] 0.1274275

###########################################################
############### LOOK AT GENES USED IN MODEL ###############

# Get the coefficients from the best model
coef_df <- coef(kfold_model$finalModel, s = kfold_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
colnames(coef_df)[1] <- "V1"
rownames(coef_df)
#  [1] "(Intercept)" "Rv0289"      "Rv0614"      "Rv0625c"     "Rv1439c"     "Rv1595"      "Rv1721c"    
# [8] "Rv2031c"     "Rv2601"      "Rv2986c"     "Rv3616c"  

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

# Print the ROC: Not sure if this is plotting in the right direction...
roc_kfold <- roc(response = kfold_model$pred$obs,
                 predictor = kfold_model$pred$Relapse,
                 levels = c("Relapse", "Cure"),  
                 direction = ">")
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
# 0.5112434
roc_sd <- sd(kfold_model$resample$ROC)
roc_sd
# 0.2549382

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
plot(roc_kfold_avg, print.auc = TRUE, main = "Average 5-fold 30-repeat CV logistic regression training set")

