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
###################### USE ALL DATA #######################
# Start with using all the data to train so I understand what is going one:

# Define training control
train_control <- trainControl(method = "LOOCV", classProbs = T, summaryFunction = twoClassSummary, savePredictions = "final")

# Train the model
set.seed(23)
lasso_loocv <- train(Outcome ~., data = my_df2, method = "glmnet", family = "binomial", metric = "ROC", trControl = train_control, tuneLength = 50)

# Check the alpha and best lambda used
lasso_loocv$bestTune # Used alpha = 1 (lasso) even though I didn't specifically tell it to
#     alpha    lambda
# 2490     1 0.2195143

lasso_loocv$results %>% arrange(desc(ROC)) %>% head()
max(lasso_loocv$results$ROC)

# Print the ROC: Not sure if this is plotting in the right direction...
roc_loocv <- roc(
  response = lasso_loocv$pred$obs,
  predictor = lasso_loocv$pred$Relapse,
  levels = c("Cure", "Relapse")
)

plot(roc_loocv, print.auc = TRUE)

# What genes did it use?
coef_df <- coef(lasso_loocv$finalModel, s = lasso_loocv$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame() 
coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
rownames(coef_df)
# It's just choosing Rv0625c again???

print(lasso_loocv)
summary(lasso_loocv)










