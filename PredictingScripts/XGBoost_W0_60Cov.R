# Try XGBoost
# E. Lamont
# 3/3/26

# This method really needs a feature selection

install.packages("xgboost")
library(xgboost)

# https://xgboost.readthedocs.io/en/stable/tutorials/model.html
# https://xgboost.readthedocs.io/en/release_1.5.0/R-package/xgboostPresentation.html
# https://xgboost.readthedocs.io/en/release_1.5.0/tutorials/param_tuning.html

###########################################################
################ TESTING AND TRAINING SETS ################

source("PredictingScripts/Separate_Test_Train.R") # Using the Mixomics method
View(myX.train)
str(myX.train)
class(myX.train) # matrix
myY.train
class(myY.train) # Needs to be a factor

# SVM want's log transformed TPM
myX.train_Log2 <- log2(myX.train + 1)
myX.test_Log2  <- log2(myX.test + 1)

dtrain <- xgb.DMatrix(data = as.matrix(myX.train_Log2), label = as.numeric(myY.train) - 1)
dtest  <- xgb.DMatrix(data = as.matrix(myX.test_Log2), label = as.numeric(myY.test) - 1)

###########################################################
##################### ONLINE TUTORIAL #####################

bstDense <- xgboost(data = myX.train_Log2, label = myY.train, 
                    max.depth = 2, # Depth of trees
                    eta = 1, 
                    nthread = 2,  # Number of CPU threads
                    nrounds = 2, # Number of passes on the data
                    objective = "binary:logistic") # Use a binary classification model
bstDense

pred <- predict(bstDense, myX.test_Log2)

# size of the prediction vector
print(length(pred))

###########################################################
######################## OTHER WAY ########################

params <- list(
  objective = "binary:logistic", # Use a binary classification model
  eval_metric = "auc",
  max_depth = 6, # Depth of trees
  eta = 0.1,
  gamma = 0,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1
)

# Cross validation
set.seed(23)
cv <- xgb.cv(
  params = params,
  data = dtrain, 
  nrounds = 200,
  nfold = 3, # 3-fold CV
  stratified = T, # To try and preserve class balances
  early_stopping_rounds = 15,
  print_every_n = 10,
  verbose = 1
)
# Multiple eval metrics are present. Will use test_auc for early stopping.
# Will train until test_auc hasn't improved in 15 rounds.
# 
# [1]	train-auc:0.919498±0.051136	test-auc:0.566799±0.112805 
# [11]	train-auc:1.000000±0.000000	test-auc:0.695767±0.072014 
# [21]	train-auc:1.000000±0.000000	test-auc:0.652116±0.086425 
# Stopping. Best iteration:
# [22]	train-auc:1.000000±0.000000	test-auc:0.652116±0.086425
# 
# [22]	train-auc:1.000000±0.000000	test-auc:0.652116±0.086425 

eval_log <- cv$evaluation_log
best_nrounds <- eval_log$iter[which.max(eval_log$test_auc_mean)]
cat("Best # of Rounds:", best_nrounds, "\n")
# Best # of Rounds: 7 

# Train the final model
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)

pred_probs <- predict(xgb_model, dtest)

# Feature importance
importance_matrix <- xgb.importance(model = xgb_model)
print(importance_matrix)
# Feature         Gain      Cover Frequency
# <char>        <num>      <num>     <num>
# 1: Rv2031c 0.2101062316 0.12074444       0.1
# 2:  Rv1430 0.1733614839 0.12144948       0.1
# 3:  Rv0733 0.1471140900 0.11149379       0.1
# 4: Rv2796c 0.1421258810 0.10852033       0.1
# 5: Rv0625c 0.1139301560 0.10436591       0.1
# 6: Rv1098c 0.1090888437 0.10182473       0.1
# 7:  Rv2601 0.1009370829 0.09968994       0.1
# 8:  Rv0041 0.0018585295 0.07498359       0.1
# 9: Rv2381c 0.0013123409 0.07482158       0.1
# 10:  Rv0002 0.0001653605 0.08210622       0.1
xgb.plot.importance(importance_matrix, top_n = 20)




