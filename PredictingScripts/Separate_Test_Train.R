# Methods to separate the testing and training sets that will work with different models 



###########################################################
##################### MIXOMICS METHOD #####################

myX <- GoodSamples60_tpmf %>% # Using TPM to start...
  # dplyr::select(-contains("THP1")) %>%
  # dplyr::select(-contains("Broth")) %>% # Lets take out broth as well
  dplyr::select(contains("W0")) %>% 
  t()

# Remove genes (columns) that are all zero
myX <- myX[, colSums(myX) != 0] # Remove columns that are all zero so the scale works for prcomp

# Want this to be Type2 I guess 
my_pipeSummary2 <- my_pipeSummary %>%
  slice(match(rownames(myX), SampleID2)) # Match the rows then just keep the one that are the same
stopifnot(all(rownames(myX) %in% my_pipeSummary2$SampleID2))
myY <- my_pipeSummary2 %>%
  pull(Type2) %>%
  factor()

summary(myY)
# W0 sputum (cure) W0 sputum (relapse) 
# 32                  11 


# set.seed(23) # 23
# train <- sample(1:nrow(myX), 30) # randomly select 30 samples in training
# test <- setdiff(1:nrow(myX), train) # rest is part of the test set
# 
# # store matrices into training and test set:
# myX.train <- myX[train, ]
# myX.test <- myX[test,]
# myY.train <- myY[train] # 7 relapses
# myY.test <- myY[test] # 4 relapses, 9 cures
# myY.test

# Make sure they are statified evenly?
set.seed(23)
relapse_idx <- which(myY == "W0 sputum (relapse)")
cure_idx    <- which(myY == "W0 sputum (cure)")
set.seed(23)
train_relapse <- sample(relapse_idx, round(0.7 * length(relapse_idx)))
train_cure    <- sample(cure_idx, round(0.7 * length(cure_idx)))

train <- c(train_relapse, train_cure)
test  <- setdiff(1:length(myY), train)

myX.train <- myX[train, ]
myX.test  <- myX[test,]

myY.train <- myY[train]
myY.test  <- myY[test]
myY.test

###########################################################
###################### CARET METHOD #######################


# Keep only W0 sputum
W0SputumSamples60_pipeSummary <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  filter(Outcome != "Failure")

P_metadata <- W0SputumSamples60_pipeSummary %>% dplyr::select(SampleID2, Outcome)
P_metadata$Outcome <- factor(P_metadata$Outcome)

P_TPM <- GoodSamples60_tpmf %>% dplyr::select(all_of(W0SputumSamples60_pipeSummary$SampleID2))

# Put everything in one dataframe
P_TPM_t <- P_TPM %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")
my_df <- inner_join(P_metadata, P_TPM_t, by = "SampleID2")
my_df2 <- my_df %>% dplyr::select(-SampleID2)

# Remove na
my_df2 <- na.omit(my_df2) # Didn't change anything

# Outcome must be factor with "Relapse" first
my_df2$Outcome <- factor(my_df2$Outcome, levels = c("Relapse", "Cure"))

# Separate training and testing sets
set.seed(23)
training_samples <- my_df2$Outcome %>% 
  createDataPartition(p = 0.7, list = FALSE) # 0.7 gives the test data 3 relapses and the train data 9 relapses
train_data  <- my_df2[training_samples, ]
validation_data <- my_df2[-training_samples, ]

# Make sure that Relapse is the positive class for caret
train_data$Outcome <- factor(train_data$Outcome, levels = c("Relapse", "Cure"))
validation_data$Outcome <- factor(validation_data$Outcome, levels = c("Relapse", "Cure"))






