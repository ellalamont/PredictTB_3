# Try this Random Forest package varSelRF
# E. Lamont
# 2/27/26

# https://www.usu.edu/math/jrstevens/stat5570/3.4.forests.pdf
# https://genesrf.iib.uam.es/
# https://www.rdocumentation.org/packages/varSelRF/versions/0.7-8/topics/varSelRF

BiocManager::install("varSelRF")
library(varSelRF)

# varSelRF: Variable selection from random forests using OOB error
# Using the OOB error as minimization criterion, carry out variable elimination from random forest, by successively eliminating the least important variables (with importance as returned from random forest).


###########################################################
################ TESTING AND TRAINING SETS ################

source("Separate_Test_Train.R") # Using the Mixomics method
View(myX.train)
myY.train


###########################################################
####################### TRY varSelRF ######################

set.seed(23)
rfsel <- varSelRF(xdata = myX.train, Class = myY.train, 
                  c.sd = 1, mtryFactor = 1, ntree = 5000,
                  ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
                  whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE,
                  returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)

rf.sig.gn <- rfsel$selected.vars
# [1] "Rv0253"  "Rv0625c" "Rv2363"  "Rv2688c"

rfsel$selected.model
# [1] "Rv0253 + Rv0625c + Rv2363 + Rv2688c"

rfsel$firstForest






