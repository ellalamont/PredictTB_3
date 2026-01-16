# K nearest neighbors
# 1/16/26

# https://www.sthda.com/english/articles/35-statistical-machine-learning-essentials/142-knn-k-nearest-neighbors-essentials/
# https://www.geeksforgeeks.org/machine-learning/k-nearest-neighbours/


# KNN assigns category based on nearby points (neighbors)
# K just tells the algorithm how many nearby points to look at when making a decision
# k-fold cross validation to choose k
# Uses distance metrics to identify nearest points. Normally Euclidean distance but there are others
# No training steps!
# May not be so good for high-dimensional data (no feature selection?)

library(caret)

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
################## COMPUTE KNN CLASSIFIER #################

# Need to choose the best value of k that minimized cross-validation error and fits the best KNN model (in caret package)
# Using Cross Validation to find the best k

set.seed(23)
train_control <- trainControl(method = "cv",
                              number = 10, # Number of folds
                              classProbs = T, # Needed for ROC/AUC
                              savePredictions = "final",
                              summaryFunction = twoClassSummary,
                              sampling = "down") # so when picking groups some relapses are always included

set.seed(23)
kNN_model <- train(Outcome ~ .,
                   data = train_data,
                   method = "knn",
                   trControl = train_control,
                   preProcess = c("center", "scale"),
                   tuneLength = 20)

# Warning in preProcess.default(thresh = 0.95, k = 5, freqCut = 19, uniqueCut = 10,  :
# These variables have zero variances: Rv2395A, Rv2561, Rv2562, Rv3098A, Rv3098c

# Plot model accuracy vs different values of k
plot(kNN_model)
# Well this looks weird

# Print the best tuning parameter k that maximizes model accuracy
kNN_model$bestTune$k
# 15

# Look at the ROC for training probabilities (not sure if this is right!)
train_probabilities <- predict(kNN_model,
                               newdata = train_data,
                               type = "prob")
roc_train <- roc(response = train_data$Outcome, 
                 predictor = train_probabilities$Relapse)
plot(roc_train, print.auc = TRUE, main = "kNN train set")

###########################################################
########### MAKE PREDICTIONS ON VALIDATION DATA ###########

valiation_probabilities <- predict(kNN_model,
                                   newdata = validation_data,
                                   type = "prob") # To get probabilities

# ROC for validation set
roc_val <- roc(response = validation_data$Outcome, predictor = valiation_probabilities$Relapse,
               direction = ">") # Not sure what direction this should be!
plot(roc_val, print.auc = TRUE, main = "kNN validation set")

# Okay well kNN is clearly not as good as k-fold logistic regression!








