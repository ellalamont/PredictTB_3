# CMA W0 60% Txn Cov
# E. Lamont
# 2/27/26

# https://link.springer.com/article/10.1186/1471-2105-9-439
# https://bioconductor.org/packages/2.3/bioc/html/CMA.html
# https://bioconductor.org/packages/2.3/bioc/vignettes/CMA/inst/doc/CMA_vignette.pdf

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
BiocManager::install("CMA")
library(CMA)

source("Import_data.R")

###########################################################
################ TESTING AND TRAINING SETS ################

source("Separate_Test_Train.R") # Using the Mixomics method
View(myX.train)
str(myX.train)
myY.train
class(myY.train) # Needs to be a factor


###########################################################
################ NOT SURE WHAT"S HAPPENING ################


fiveCV5iter <- GenerateLearningsets(y = myY.train, method = "CV",
                                    fold = 5, niter = 5, strat = TRUE)

class_dlda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                             classifier = dldaCMA)

genesel_da <- GeneSelection(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                            method = "t.test", scheme = "one-vs-all")

class_lda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                            classifier = ldaCMA, genesel = genesel_da, nbgene = 10)

class_fda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                            classifier = fdaCMA, genesel = genesel_da, nbgene = 10, comp = 2)

class_qda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                            classifier = qdaCMA, genesel = genesel_da, nbgene = 1)

class_scda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                             classifier = scdaCMA, tuninglist = list(grids = list()))

class_plsda <- classification(X = myX.train, y = myY.train, learningsets = fiveCV5iter,
                              classifier = pls_ldaCMA)

dalike <- list(class_dlda, class_lda, class_fda, class_qda)
par(mfrow = c(3, 1))
comparison <- compare(dalike, plot = TRUE, measure = c("misclassification",
                                                         "brier score", "average probability"))
print(comparison)


###########################################################
######################## TRYING AGAIN #####################

# Log2-transform RNA-seq TPMs
myX_log <- log2(myX + 1)

# Standardize (z-score)
myX_scaled <- scale(myX_log)

# Make sure Y is a factor with relapse/cure labels
myY <- factor(myY, levels = c("W0 sputum (cure)", "W0 sputum (relapse)"))

# 2️⃣ Stratified train/test split (~70/30)
# -------------------------------
set.seed(23)

train_idx <- createDataPartition(myY, p = 0.7, list = FALSE)
myX.train <- myX_scaled[train_idx, ]
myX.test  <- myX_scaled[-train_idx, ]
myY.train <- myY[train_idx]
myY.test  <- myY[-train_idx]

table(myY.train)
table(myY.test)

# -------------------------------
# 3️⃣ Generate learning sets (nested CV)
# -------------------------------
set.seed(23)
fiveCV5iter <- GenerateLearningsets(
  y = myY.train,
  method = "CV",
  fold = 5,
  niter = 5,
  strat = TRUE   # preserves class balance in each fold
)

# -------------------------------
# 4️⃣ Gene selection (optional but recommended)
# -------------------------------
# Example: t-test, one-vs-all
genesel_da <- GeneSelection(
  X = myX.train,
  y = myY.train,
  learningsets = fiveCV5iter,
  method = "t.test",
  scheme = "one-vs-all"
)

# -------------------------------
# 5️⃣ Train classifiers
# -------------------------------
# Discriminant classifiers
class_dlda <- classification(X = myX.train, y = myY.train,
                             learningsets = fiveCV5iter,
                             classifier = dldaCMA)

class_lda <- classification(X = myX.train, y = myY.train,
                            learningsets = fiveCV5iter,
                            classifier = ldaCMA,
                            genesel = genesel_da,
                            nbgene = 10)

class_fda <- classification(X = myX.train, y = myY.train,
                            learningsets = fiveCV5iter,
                            classifier = fdaCMA,
                            genesel = genesel_da,
                            nbgene = 10,
                            comp = 2)

class_qda <- classification(X = myX.train, y = myY.train,
                            learningsets = fiveCV5iter,
                            classifier = qdaCMA,
                            genesel = genesel_da,
                            nbgene = 1)

class_scda <- classification(X = myX.train, y = myY.train,
                             learningsets = fiveCV5iter,
                             classifier = scdaCMA,
                             tuninglist = list(grids = list()))

class_plsda <- classification(X = myX.train, y = myY.train,
                              learningsets = fiveCV5iter,
                              classifier = pls_ldaCMA)
# Error in pls.regression(Xlearn, transformy(Ylearn + 1), ncomp = comp) : 
# could not find function "pls.regression"


# -------------------------------
# 6️⃣ Compare performance
# -------------------------------
dalike <- list(class_dlda, class_lda, class_fda, class_qda, class_scda)

par(mfrow = c(3, 1))
comparison <- compare(dalike, plot = TRUE,
                      measure = c("misclassification", "brier score", "average probability"))
print(comparison)

# -------------------------------
# 7️⃣ Optional: Evaluate on test set
# -------------------------------
# Example for one classifier
pred_test <- predict(class_lda, newdata = myX.test)

conf_mat <- table(Predicted = pred_test$class, True = myY.test)
print(conf_mat)

roc_obj <- pROC::roc(myY.test, as.numeric(pred_test$prob[,2]))  # probability for relapse
auc_val <- pROC::auc(roc_obj)
cat("Test set AUC:", auc_val, "\n")

