# Logistic Regression W0 data
# 1/12/26

# https://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# https://www.sthda.com/english/articles/36-classification-methods-essentials/150-stepwise-logistic-regression-essentials-in-r/
# https://www.sthda.com/english/articles/36-classification-methods-essentials/143-evaluation-of-classification-model-accuracy-essentials/


# 38 cure, 12 relapse

source("Import_data.R") 

# install.packages("caret")
library(caret)
# install.packages("glmnet")
library(glmnet)
# install.packages("pROC")
library(pROC)


###########################################################
#################### STHDA TUTORIAL #######################

install.packages("mlbench")
library(mlbench)
data("PimaIndiansDiabetes2", package = "mlbench")
head(PimaIndiansDiabetes2)
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)
sample_n(PimaIndiansDiabetes2, 3)
set.seed(123)
Xtraining.samples <- PimaIndiansDiabetes2$diabetes %>% 
  createDataPartition(p = 0.8, list = FALSE)
Xtrain.data  <- PimaIndiansDiabetes2[Xtraining.samples, ]
Xtest.data <- PimaIndiansDiabetes2[-Xtraining.samples, ]
# Dumy code categorical predictor variables
Xx <- model.matrix(diabetes~., Xtrain.data)[,-1]
# Convert the outcome (class) to a numerical variable
Xy <- ifelse(Xtrain.data$diabetes == "pos", 1, 0)
set.seed(123) 
Xcv.lasso <- cv.glmnet(Xx, Xy, alpha = 1, family = "binomial")
plot(Xcv.lasso)
Xcv.lasso$lambda.min
coef(Xcv.lasso, Xcv.lasso$lambda.min)

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
cv.lasso$lambda.1se # 0.1826832. More conservative, fewer predictors, maybe better for biology

# Look at the regression coefficients
coef(cv.lasso, cv.lasso$lambda.1se)

coef.1se_df <- as.data.frame(as.matrix(coef(cv.lasso, s = "lambda.1se"))) %>% 
  filter(lambda.1se != 0)
View(coef.1se_df)
# Only keeps 4 genes: Rv0046c, Rv0345, Rv0625c, Rv1249c

coef.min_df <- as.data.frame(as.matrix(coef(cv.lasso, s = "lambda.min"))) %>% 
  filter(lambda.min != 0)
View(coef.min_df)
# Keeps 12 genes... so proceed with this because 12 isn't that many....
# [1] "(Intercept)" "Rv0046c"     "Rv0345"      "Rv0361"      "Rv0625c"     "Rv1031"      "Rv1032c"    
# [8] "Rv1083"      "Rv1249c"     "Rv2104c"     "Rv2986c"     "Rv3219"      "Rv3786c"  