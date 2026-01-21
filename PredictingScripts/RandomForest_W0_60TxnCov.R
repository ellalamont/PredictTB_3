# Random Forest
# 1/16/26

# https://www.sthda.com/english/articles/35-statistical-machine-learning-essentials/140-bagging-and-random-forest-essentials/
# https://www.usu.edu/math/jrstevens/stat5570/3.4.forests.pdf
# https://genesrf.iib.uam.es/
# https://www.rdocumentation.org/packages/varSelRF/versions/0.7-8/topics/varSelRF


library(caret)
# install.packages("randomForest")
library(randomForest)

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

# Remove genes with near zero variance
NearZeroGenes <- nearZeroVar(my_df2 %>% select(-Outcome))
my_df3 <- my_df2 %>% select(-NearZeroGenes) # Now 4025 genes

# Saving for trying in https://genesrf.iib.uam.es/
# my_df_tosave <- my_df %>% select(-NearZeroGenes) %>% t()
# write.table(my_df_tosave, file = "Data_for_online_RF_Tool.txt", sep = "\t")

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
################## COMPUTE RF CLASSIFIER ##################

# Selecting the optimal number (mtry) of predictor variables using caret

set.seed(23)
RF_model <- train(Outcome ~ .,
                  data = train_data,
                  method = "rf",
                  trControl = trainControl("cv", number = 10),
                  importance = T)
# Warning in preProcess.default(thresh = 0.95, k = 5, freqCut = 19, uniqueCut = 10,  :
# These variables have zero variances: Rv2395A, Rv2561, Rv2562, Rv3098A, Rv3098c
# Warning message:
# In nominalTrainWorkflow(x = x, y = y, wts = weights, info = trainInfo,  :
#                           There were missing values in resampled performance measures.

# Best tuning parameter
RF_model$bestTune
# 2

# Final model
RF_model$finalModel
# Wow it really didn't do very well

###########################################################
########### MAKE PREDICTIONS ON VALIDATION DATA ###########

validation_classes <- predict(RF_model,
                                   newdata = validation_data)


###########################################################
############### FIND IMPORTANCE OF VARIABLES ##############

importance(RF_model$finalModel)
# MeanDecreaseAccuracy, which is the average decrease of model accuracy in predicting the outcome of the out-of-bag samples when a specific variable is excluded from the model.
# MeanDecreaseGini, which is the average decrease in node impurity that results from splits over that variable. The Gini impurity index is only used for classification problem. In the regression the node impurity is measured by training set RSS. These measures, calculated using the training set, are less reliable than a measure calculated on out-of-bag data. See Chapter @ref(decision-tree-models) for node impurity measures (Gini index and RSS).

# Plot MeanDecreaseAccuracy
varImpPlot(RF_model$finalModel, type = 1)
# Plot MeanDecreaseGini
varImpPlot(RF_model$finalModel, type = 2)
# More to the right means more important

# Show importance of variables in percentage
varImp(RF_model)

