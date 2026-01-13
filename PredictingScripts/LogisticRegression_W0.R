# Logistic Regression W0 data
# 1/12/26

# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# https://www.sthda.com/english/articles/36-classification-methods-essentials/150-stepwise-logistic-regression-essentials-in-r/
# https://www.sthda.com/english/articles/36-classification-methods-essentials/143-evaluation-of-classification-model-accuracy-essentials/

# https://glmnet.stanford.edu/articles/glmnet.html


# Could I train on all these data, then test on W2 or the partial transcriptomes?

# 38 cure, 12 relapse

source("Import_data.R") 

# install.packages("caret")
library(caret)
# install.packages("glmnet")
library(glmnet)
# install.packages("pROC")
library(pROC)

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
my_df2 <- my_df %>% select(-SampleID2) # I think I might need the dataframe without the SampleID2...? Because the tutorial data doesn't have this

# Remove na
my_df2 <- na.omit(my_df2) # Didn't change anything


###########################################################
############### SPLIT INTO TEST/TRAIN TEST ################

set.seed(42)
training.samples <- my_df2$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train.data  <- my_df2[training.samples, ]
test.data <- my_df2[-training.samples, ]


###########################################################
############### CREATE MATRIX OF PREDICTORS ###############

# Dumy code categorical predictor variables
x <- model.matrix(Outcome ~ ., train.data)[,-1]
# Convert the outcome (class) to a numerical variable
y <- ifelse(train.data$Outcome == "Relapse", 1, 0)

###########################################################
####################### FIND LAMBDA #######################

# The larger the Lamba the more coefficients are shrunk towards zero
# Use cross validation to find the best lamba

set.seed(42) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
# Warning messages:
#   1: In lognet(x, is.sparse, y, weights, offset, alpha, nobs, nvars,  :
#                  one multinomial or binomial class has fewer than 8  observations; dangerous ground

plot(cv.lasso)
cv.lasso$lambda.min # 0.05202883. Gives the most accurate lamba
cv.lasso$lambda.1se # 0.1826832. More conservative/simple, fewer predictors, may be less accurate

# Look at the regression coefficients
coef(cv.lasso, cv.lasso$lambda.1se)

coef.1se_df <- as.data.frame(as.matrix(coef(cv.lasso, s = "lambda.1se"))) %>% 
  filter(lambda.1se != 0)
# View(coef.1se_df)
# Only keeps 4 genes: Rv0046c, Rv0345, Rv0625c, Rv1249c

coef.min_df <- as.data.frame(as.matrix(coef(cv.lasso, s = "lambda.min"))) %>% 
  filter(lambda.min != 0)
# View(coef.min_df)
# Keeps 12 genes... so proceed with this because 12 isn't that many....
# [1] "(Intercept)" "Rv0046c"     "Rv0345"      "Rv0361"      "Rv0625c"     "Rv1031"      "Rv1032c"    
# [8] "Rv1083"      "Rv1249c"     "Rv2104c"     "Rv2986c"     "Rv3219"      "Rv3786c"  

###########################################################
################### COMPUTE LASSO MODEL ###################

# Final model with lambda.min
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)

# Make prediction on test data
x.test <- model.matrix(Outcome ~., test.data)[,-1]
probabilities <- lasso.model %>% predict(newx = x.test) # This gives log odds?
predicted.classes <- ifelse(probabilities > 0.5, "Relapse", "Cure")
# Model accuracy
observed.classes <- test.data$Outcome
mean(predicted.classes == observed.classes)
# It guessed they would all be cure...

# Adjust some arguments:
probabilities2 <- lasso.model %>% predict(newx = x.test, type = "response") # By setting type = "response" this now gives probabilities instead of log odds?
predicted.classes2 <- ifelse(probabilities2 > 0.2, "Relapse", "Cure") # Change this value a bit because relapses are so rare
observed.classes <- test.data$Outcome
mean(predicted.classes2 == observed.classes) # 0.8333333

###########################################################
############# HOW WELL THE LASSO MODEL WORKED #############

# Classification Accuracy: Proportion of observations that have been correctly classified
# Classification Error Rate: Proportion of observations that have been misclassified (Error Rate = 1 - accuracy)
Accuracy <- mean(observed.classes == predicted.classes2) # 0.8333333
Error <- mean(observed.classes != predicted.classes2) # 0.1666667

# Confusion matrix: Numbers
table(observed.classes, predicted.classes2)
# Confusion matrix: proportion
table(observed.classes, predicted.classes2) %>% prop.table() %>% round(digits = 3)

# Other metrics
predicted.classes3 <- as.vector(predicted.classes2) %>% factor()
confusionMatrix(predicted.classes3, observed.classes, positive = "Relapse")

# ROC and AUC
res.roc <- roc(observed.classes, as.vector(probabilities2))
plot.roc(res.roc, print.auc = TRUE)
res.roc <- roc(observed.classes, as.vector(probabilities))
plot.roc(res.roc, print.auc = TRUE)

###########################################################
################### COMPUTE RIDGE MODEL ###################

set.seed(42) 
cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial")
# Warning messages:
# 1: In lognet(x, is.sparse, y, weights, offset, alpha, nobs, nvars,  :
#                one multinomial or binomial class has fewer than 8  observations; dangerous ground
plot(cv.ridge)
cv.ridge$lambda.min # 191.3818. Gives the most accurate lamba
cv.ridge$lambda.1se # 334.4451. More conservative/simple, fewer predictors, may be less accurate

# Look at the regression coefficients
# Actually don't bother because ridge regression keeps everything

# With lambda.min and alpha = 0 which is ridge regression (uses all variables)
ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min)
ridge_probabilities <- ridge.model %>% predict(newx = x.test, type = "response") # By setting type = "response" this now gives probabilities instead of log odds?
ridge_predicted.classes <- ifelse(ridge_probabilities > 0.2, "Relapse", "Cure") # Change this value a bit because relapses are so rare
observed.classes <- test.data$Outcome
mean(ridge_predicted.classes == observed.classes) # 0.25

# So ridge regression clearly didn't work... model is overfit with using all the genes...

# Try ridge model with weights, give more weight to the relapse cases(?)
w <- ifelse(y == 1, 299/21, 1)
ridge.model_weighted <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min, weights = w)
ridge_probabilities_weighted <- ridge.model_weighted %>% predict(newx = x.test, type = "response")
# This clearly isn't going to be any better. Ridge regression won't work! 

###########################################################
############ STEPWISE LOGISTIC REGRESSION - BAD ###########

# I don't think this is so good when p>>n but will try anyway for practice

# Have to do the full model first??
full.model <- glm(Outcome ~., data = train.data, family = binomial)
coef(full.model)

library(MASS)
# step.model <- full.model %>% stepAIC(trace = FALSE) # Had to stop this, was taking too long
coef(step.model)


###########################################################
############# LASSO LOGISTIC W/PARTIAL AS TEST ############
#1/13/26

# Use all the complete transcriptomes for the train set, then test on the partial transcriptomes
# There should be 13 cure and 7 relapse in the partial transcriptomes

# Grab the partial transcriptomes
W0SputumSamplesPartial_metadata <- my_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(Outcome != "Failure") %>%
  filter(!SampleID2 %in% W0SputumSamples60_pipeSummary$SampleID2) %>% # Remove everything that is a full txn
  filter(N_Genomic > 100000) %>% # Just cause I really don't want to use the lower ones
  mutate(Outcome = ifelse(Outcome == "Prob Relapse", "Relapse", Outcome)) %>%
  select(SampleID2, Outcome)
# Okay now there are 8 cures and 6 relapses. Is this ratio bad?

W0SputumSamplesPartial_metadata$Outcome <- factor(W0SputumSamplesPartial_metadata$Outcome)

Partial_TPM <- All_tpm_f %>% select(all_of(W0SputumSamplesPartial_metadata$SampleID2))

# Put everything in one dataframe
Partial_TPM_t <- Partial_TPM %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
Partial_df <- inner_join(W0SputumSamplesPartial_metadata, Partial_TPM_t, by = "SampleID2")
Partial_df2 <- Partial_df %>% select(-SampleID2) # I think I might need the dataframe without the SampleID2...? Because the tutorial data doesn't have this

# Remove na
Partial_df2 <- na.omit(Partial_df2) # Didn't change anything

# Create matrix of predictors (Still not clear on this step...)
x_Full <- model.matrix(Outcome ~ ., my_df2)[,-1]
# Convert the outcome (class) to a numerical variable
y_Full <- ifelse(my_df2$Outcome == "Relapse", 1, 0)

# Determine the best Lambda
set.seed(42) 
cv.lasso_Full <- cv.glmnet(x_Full, y_Full, alpha = 1, family = "binomial")
plot(cv.lasso_Full)
cv.lasso_Full$lambda.min # 0.1695825. Gives the most accurate lamba
cv.lasso_Full$lambda.1se # 0.2700233. More conservative/simple, fewer predictors, may be less accurate
coef.min_df_Full <- as.data.frame(as.matrix(coef(cv.lasso_Full, s = "lambda.min"))) %>% filter(lambda.min != 0)
rownames(coef.min_df_Full) # "Rv0523c"     "Rv0625c"     "Rv1775"      "Rv2019"      "Rv2331A"     "Rv2896c"  
coef.1se_df_Full <- as.data.frame(as.matrix(coef(cv.lasso_Full, s = "lambda.1se"))) %>% filter(lambda.1se != 0)
rownames(coef.1se_df_Full) # "Rv0625c" Wow only one don't use this one! 

# Final model with lambda.min
lasso.model_Full <- glmnet(x_Full, y_Full, alpha = 1, family = "binomial", lambda = cv.lasso_Full$lambda.min)

# Make prediction on test data
x.test_Partial <- model.matrix(Outcome ~., Partial_df2)[,-1]
probabilities_Partial <- lasso.model_Full %>% predict(newx = x.test_Partial, type = "response") # By setting type = "response" this now gives probabilities instead of log odds?
predicted.classes_Partial <- ifelse(probabilities_Partial > 0.2, "Relapse", "Cure") # Change this value a bit because relapses are so rare
observed.classes_Partial <- Partial_df2$Outcome
mean(predicted.classes_Partial == observed.classes_Partial) # 0.3571429
# Wow this was really bad!!





