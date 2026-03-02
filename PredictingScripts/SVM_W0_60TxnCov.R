# Support Vector Machines!
# E. Lamont 
# 3/2/26

# https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.28567
# https://stackoverflow.com/questions/17529537/example-for-svm-feature-selection-in-r
# https://www.geeksforgeeks.org/machine-learning/svm-feature-selection-in-r-with-example/


# RFE for Feature Selection?

library(caret)
library(e1071)

###########################################################
################ TESTING AND TRAINING SETS ################

source("Separate_Test_Train.R") # Using the Mixomics method
View(myX.train)
str(myX.train)
myY.train
class(myY.train) # Needs to be a factor

# SVM want's log transformed TPM
myX.train_Log2 <- log2(myX.train + 1)
myX.test_Log2  <- log2(myX.test + 1)

# SVM required scaling
myX.train_scaled <- scale(myX.train_Log2)
myX.test_scaled <- scale(myX.test_Log2,
  center = attr(myX.train_scaled, "scaled:center"),
  scale  = attr(myX.train_scaled, "scaled:scale"))



###########################################################
################## RFE FEATURE SELECTION ##################

# Define RFE control setup
control <- rfeControl(functions = caretFuncs,  # Default caret functions for models
                      method = "repeatedcv",  # Cross-validation
                      number = 3,  # Number of folds in cross-validation
                      repeats = 10,
                      verbose = FALSE)  # Control output verbosity

# Running RFE with SVM using a linear kernel
svm_rfe <- caret::rfe(x = myX.train_scaled, 
               y = myY.train,
               sizes = c(2, 5, 10, 20),  # Number of features to select
               rfeControl = control,
               method = "svmLinear")  # Ensures we are using SVM with a linear kernel, could also be svmRadial?

# Error in { : task 1 failed - "argument 1 is not a vector"

# STOPPED! Keeps crashing because there are too many genes!

###########################################################
################## PCA FEATURE SELECTION ##################

pca_res <- prcomp(myX.train_scaled, center = TRUE, scale. = TRUE)

# Choose top PCs explaining ~80% variance (or top N)
explained_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
num_pcs <- min(which(explained_var >= 0.8))  # e.g., 80% variance
cat("Using top", num_pcs, "PCs for SVM\n")

train_pcs <- as.data.frame(pca_res$x[, 1:num_pcs])

# Project test set onto same PCs
test_pcs <- as.data.frame(predict(pca_res, newdata = myX.test_scaled)[, 1:num_pcs])


###########################################################
##################### TRAIN SVM ON PCA ####################

svm_control <- trainControl(
  method = "repeatedcv",
  number = 3,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

svm_model <- train(
  x = train_pcs,
  y = myY.train,
  method = "svmLinear",
  trControl = svm_control,
  metric = "ROC"  # maximize AUC
)

# Still not working... gave up









