# PERMDISP Variation 
# 4/13/26
# E. Lamont

# "All happy families are alike; each unhappy family is unhappy in its own way."


# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/
# "PERMDISP is a multivariate extension of Levene’s test (Anderson 2006) to examine whether groups differ in plot-to-plot variability"
# "In essence, PERMDISP involves calculating the distance from each data point to its group centroid and then testing whether those distances differ among the groups."

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9))

source("Import_data.R")
# GoodSputum60_tpmf_log2
# GoodSputum60_Metadata

source("TimeCourseScripts/Delta_Weeks0to2.R")

library(vegan)


###########################################################
################## W0 Log2(TPM) euclidean #################

meta_W0 <- GoodSputum60_Metadata %>% filter(Week == "Week 0")

samples_W0 <- meta_W0$SampleID
log2TPM_W0 <- GoodSputum60_tpmf_log2[, samples_W0]

# remove genes with zero variance
log2TPM_W0_f <- log2TPM_W0[apply(log2TPM_W0, 1, var) != 0, ]
# Now only 4025 genes

log2TPM_W0_t <- t(log2TPM_W0_f)

# Make a distance matrix
dist_W0 <- vegdist(log2TPM_W0_t, method = "euclidean")

# Run PERMDISP
bd_W0 <- betadisper(dist_W0, meta_W0$Outcome)
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W0, group = meta_W0$Outcome)
# 
# No. of Positive Eigenvalues: 42
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Cure Relapse 
# 81.68   86.04 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 42 eigenvalues)
# PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
# 53840 31672 25172 18976 16995 14535 13005 12807 

permutest(bd_W0)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq     F N.Perm Pr(>F)
# Groups     1    155  155.31 0.128    999  0.728
# Residuals 41  49752 1213.46        

boxplot(bd_W0, main = "W0 Distance to Centroid")

plot_df_W0 <- data.frame(
  distance = bd_W0$distances,
  Outcome = meta_W0$Outcome,
  Timepoint = "W0")


###########################################################
################# W0 Log2(TPM) CORRELATION ################

meta_W0 <- GoodSputum60_Metadata %>% filter(Week == "Week 0")

# If only want to look at the ones that have matching W2, do this:
meta_W0 <- meta_W0 %>% dplyr::filter(SampleID %in% W0_subset_names) 

samples_W0 <- meta_W0$SampleID
log2TPM_W0 <- GoodSputum60_tpmf_log2[, samples_W0]

# remove genes with zero variance
log2TPM_W0_f <- log2TPM_W0[apply(log2TPM_W0, 1, var) != 0, ]
# Now only 4025 genes

# Make a correlation distance matrix
cor_mat_W0 <- cor(log2TPM_W0_f, method = "pearson")
dist_W0 <- as.dist(1 - cor_mat_W0)

# Run PERMDISP
bd_W0 <- betadisper(dist_W0, meta_W0$Outcome, type = "centroid", add = TRUE)
bd_W0
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W0, group = meta_W0$Outcome, type = "centroid", add
#                  = TRUE)
# 
# No. of Positive Eigenvalues: 41
# No. of Negative Eigenvalues: 0
# 
# Average distance to centroid:
#   Cure Relapse 
# 0.2882  0.3266 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 41 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 0.4859 0.3470 0.2508 0.1985 0.1786 0.1659 0.1489 0.1434 

permutest(bd_W0)
permutest_bd_W0 <- permutest(bd_W0)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.010736 0.0107357 1.8291    999  0.195
# Residuals 41 0.240640 0.0058693           

boxplot(bd_W0, main = "W0 Distance to Centroid")
plot(bd_W0)

anova(bd_W0)
# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     1 0.012025 0.0120247  2.4681 0.1239
# Residuals 41 0.199749 0.0048719 

# PERMANOVA
adonis2(dist_W0 ~ meta_W0$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_W0 ~ meta_W0$Outcome)
# Df SumOfSqs      R2     F Pr(>F)
# Model     1   0.0684 0.02641 1.112  0.282
# Residual 41   2.5219 0.97359             
# Total    42   2.5903 1.00000 

# Of distances?
adonis2(dist(bd_W0$distances) ~ meta_W0$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist(bd_W0$distances) ~ meta_W0$Outcome)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1 0.012025 0.05678 2.4681   0.14
# Residual 41 0.199749 0.94322              
# Total    42 0.211774 1.00000  

###########################################################
########## W0 Log2(TPM) CORRELATION - PCoA FIG ############

# Pull out the distances of each sample
scores_W0 <- as.data.frame(bd_W0$vectors)
scores_W0$SampleID <- rownames(scores_W0)
scores_W0$Outcome <- meta_W0$Outcome

# % variation explained = proportion of distance variation captured
eig <- bd_W0$eig
var_exp <- eig / sum(eig)
percent_var <- round(var_exp * 100, 1)
percent_var[1:5]

# Determine where the centroids are 
centroids_W0 <- scores_W0 %>%
  dplyr::group_by(Outcome) %>%
  dplyr::summarise(PCoA1 = mean(PCoA1), 
                   PCoA2 = mean(PCoA2))

scores_W0 <- merge(scores_W0, centroids_W0, by = "Outcome", suffixes = c("", "_centroid"))

# Make hull df for polygon on plot
hulls <- scores_W0 %>%
  group_by(Outcome) %>%
  slice(chull(PCoA1, PCoA2))

# Make the figure
W0_PermDisp_PCoA_fig1 <- scores_W0 %>% 
  ggplot(aes(x = PCoA1, y = PCoA2, fill = Outcome, shape = Outcome)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8) +
  geom_polygon(data = hulls, aes(color = Outcome, group = Outcome), 
               fill = NA, linewidth = 0.6, show.legend = F) +
  geom_segment(data = scores_W0, aes(xend = PCoA1_centroid, yend = PCoA2_centroid, color = "black"), alpha = 0.3) +
  geom_point(data = centroids_W0, aes(x = PCoA1, y = PCoA2, color = Outcome), 
             size = 4, shape = 4, stroke = 1.5, show.legend = F) +
  scale_fill_manual(values = c(`Cure` = "#0072B2", `Relapse` = "#bc5300")) +  
  scale_color_manual(values = c(`Cure` = "#0072B2", `Relapse` = "#bc5300")) +  
  scale_shape_manual(values = c(`Cure` = 21, `Relapse` = 21)) + 
  labs(title = "W0 Subset: Log2(TPM+1) → Pearson Corr → PERMDISP",
       subtitle = paste0("Avg dist to centroid: Cure=", round(bd_W0$group.distances[[1]], 3), "; Relapse=", round(bd_W0$group.distances[[2]], 3), "; p=", permutest_bd_W0$tab$`Pr(>F)`[1]),
       x = "PCoA1",
       y = "PCoA2") +
  my_plot_themes
W0_PermDisp_PCoA_fig1
# ggsave(W0_PermDisp_PCoA_fig1,
#        file = paste0("W0Subset_Corr_PCoA_v1.png"),
#        path = "Figures/PERMDISP",
#        dpi = 600,
#        width = 8, height = 6, units = "in")


###########################################################
##################### W0 TPM euclidean ####################

meta_W0 <- GoodSputum60_Metadata %>% filter(Week == "Week 0")

samples_W0 <- meta_W0$SampleID
TPM_W0 <- GoodSputum60_tpmf[, samples_W0]

# remove genes with zero variance
TPM_W0_f <- TPM_W0[apply(TPM_W0, 1, var) != 0, ]
# Now only 4025 genes

TPM_W0_t <- t(TPM_W0_f)

# Make a distance matrix
dist_W0 <- vegdist(TPM_W0_t, method = "euclidean")

# Run PERMDISP
bd_W0 <- betadisper(dist_W0, meta_W0$Outcome)
bd_W0
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W0, group = meta_W0$Outcome)
# 
# No. of Positive Eigenvalues: 42
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Cure Relapse 
# 16591   17243 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 42 eigenvalues)
# PCoA1     PCoA2     PCoA3     PCoA4     PCoA5     PCoA6     PCoA7     PCoA8 
# 3.135e+09 1.884e+09 1.387e+09 9.500e+08 7.998e+08 7.314e+08 6.768e+08 6.299e+08 

permutest(bd_W0)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df     Sum Sq  Mean Sq     F N.Perm Pr(>F)
# Groups     1    3475214  3475214 0.083    999  0.773
# Residuals 41 1715882047 41850782   

boxplot(bd_W0, main = "W0 Distance to Centroid")

plot_df_W0 <- data.frame(
  distance = bd_W0$distances,
  Outcome = meta_W0$Outcome,
  Timepoint = "W0")



###########################################################
##################### W0 TPM CORRELATION ##################

meta_W0 <- GoodSputum60_Metadata %>% filter(Week == "Week 0")

samples_W0 <- meta_W0$SampleID
TPM_W0 <- GoodSputum60_tpmf[, samples_W0]

# remove genes with zero variance
TPM_W0_f <- TPM_W0[apply(TPM_W0, 1, var) != 0, ]
# Now only 4025 genes

# Make a correlation distance matrix
cor_mat_W0 <- cor(TPM_W0_f, method = "pearson")
dist_W0 <- as.dist(1 - cor_mat_W0)

# Run PERMDISP
bd_W0 <- betadisper(dist_W0, meta_W0$Outcome, type = "centroid", add = TRUE)
bd_W0
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W0, group = meta_W0$Outcome, type = "centroid", add
#                  = TRUE)
# 
# No. of Positive Eigenvalues: 41
# No. of Negative Eigenvalues: 0
# 
# Average distance to centroid:
#   Cure Relapse 
# 0.2836  0.3167 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 41 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 0.6266 0.4537 0.2997 0.2473 0.1723 0.1511 0.1325 0.1264 

permutest(bd_W0)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00895 0.0089500 1.7387    999  0.206
# Residuals 41 0.21105 0.0051475           

boxplot(bd_W0, main = "W0 Distance to Centroid")
plot(bd_W0)

anova(bd_W0)
# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups     1 0.00895 0.0089500  1.7387 0.1946
# Residuals 41 0.21105 0.0051475  

# PERMANOVA
adonis2(dist_W0 ~ meta_W0$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_W0 ~ meta_W0$Outcome)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.07069 0.04149 1.7745  0.103
# Residual 41  1.63332 0.95851              
# Total    42  1.70401 1.00000    

# Of distances?
adonis2(dist(bd_W0$distances) ~ meta_W0$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist(bd_W0$distances) ~ meta_W0$Outcome)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.00895 0.04068 1.7387  0.201
# Residual 41  0.21105 0.95932              
# Total    42  0.22000 1.00000    

plot_df_W0 <- data.frame(
  distance = bd_W0$distances,
  Outcome = meta_W0$Outcome,
  Timepoint = "W0")


###########################################################
################# W2 Log2(TPM) euclidean ##################

meta_W2 <- GoodSputum60_Metadata %>% filter(Week == "Week 2")

samples_W2 <- meta_W2$SampleID
log2TPM_W2 <- GoodSputum60_tpmf_log2[, samples_W2]

# remove genes with zero variance
log2TPM_W2_f <- log2TPM_W2[apply(log2TPM_W2, 1, var) != 0, ]
# Now only 4024 genes

log2TPM_W2_t <- t(log2TPM_W2_f)

# Make a distance matrix
dist_W2 <- vegdist(log2TPM_W2_t, method = "euclidean")
# dist_W2 <- vegdist(log2TPM_W2_t, method = "bray")

# Run PERMDISP
bd_W2 <- betadisper(dist_W2, meta_W2$Outcome)
bd_W2
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W2, group = meta_W2$Outcome)
# 
# No. of Positive Eigenvalues: 12
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Cure Relapse 
# 142.8   176.3 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 12 eigenvalues)
# PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
# 74491 56239 44510 34829 33034 29678 27943 19766  

permutest(bd_W2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1   3447    3447 2.3822    999  0.148
# Residuals 11  15917    1447  

boxplot(bd_W2, main = "W2 Distance to Centroid")
plot(bd_W2)

plot_df_W2 <- data.frame(
  distance = bd_W2$distances,
  Outcome = meta_W2$Outcome,
  Timepoint = "W2")


###########################################################
################# W2 Log2(TPM) CORRELATION ################
# ** THIS ONE **

meta_W2 <- GoodSputum60_Metadata %>% filter(Week == "Week 2")

samples_W2 <- meta_W2$SampleID
log2TPM_W2 <- GoodSputum60_tpmf_log2[, samples_W2]

# remove genes with zero variance
log2TPM_W2_f <- log2TPM_W2[apply(log2TPM_W2, 1, var) != 0, ]
# Now only 4024 genes

# Make a correlation distance matrix
cor_mat_W2 <- cor(log2TPM_W2_f, method = "pearson")
dist_W2 <- as.dist(1 - cor_mat_W2)


# Run PERMDISP
bd_W2 <- betadisper(dist_W2, meta_W2$Outcome, type = "centroid", add = TRUE)
bd_W2
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W2, group = meta_W2$Outcome, type = "centroid", add
#                  = TRUE)
# 
# No. of Positive Eigenvalues: 12
# No. of Negative Eigenvalues: 0
# 
# Average distance to centroid:
#   Cure Relapse 
# 0.3830  0.4845 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 12 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 0.4837 0.3327 0.2763 0.2592 0.2239 0.2135 0.2071 0.1938 

set.seed(42)
permutest(bd_W2)
set.seed(42)
permutest_bd_W2 <- permutest(bd_W2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups     1 0.031673 0.031673 10.065    999  0.009 **
#   Residuals 11 0.034615 0.003147                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1  

boxplot(bd_W2, main = "W2 Distance to Centroid")
plot(bd_W2)

anova(bd_W2)
# Analysis of Variance Table
# 
# Response: Distances
# Df   Sum Sq  Mean Sq F value   Pr(>F)   
# Groups     1 0.031673 0.031673  10.065 0.008879 **
#   Residuals 11 0.034615 0.003147                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# PERMANOVA
adonis2(dist_W2 ~ meta_W2$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_W2 ~ meta_W2$Outcome)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1  0.27644 0.10399 1.2766  0.043 *
#   Residual 11  2.38194 0.89601                
# Total    12  2.65838 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Of distances?
adonis2(dist(bd_W2$distances) ~ meta_W2$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist(bd_W2$distances) ~ meta_W2$Outcome)
# Df SumOfSqs     R2      F Pr(>F)   
# Model     1 0.031673 0.4778 10.065  0.005 **
#   Residual 11 0.034615 0.5222                 
# Total    12 0.066288 1.0000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Trying with the equation from https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering ##
# dist2_W2 <- as.dist(sqrt(2*(1 - cor_mat_W2)))
dist2_W2 <- as.dist(sqrt(1 - cor_mat_W2))
bd2_W2 <- betadisper(dist2_W2, meta_W2$Outcome, type = "centroid", add = FALSE)
bd2_W2
permutest(bd2_W2)
# This also works?

### CAP ###
# Not sure what is going on here...... Think should try this for prediction
# cap_W2 <- capscale(dist_W2 ~ Outcome, data = meta_W2)
# plot(cap_W2)
# anova(cap_W2, permutations = 999)
# # Permutation test for capscale under reduced model
# # Permutation: free
# # Number of permutations: 999
# # 
# # Model: capscale(formula = dist_W2 ~ Outcome, data = meta_W2)
# # Df SumOfSqs      F Pr(>F)   
# # Model     1  0.33942 2.1232   0.01 **
# #   Residual 11  1.75852                 
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# # get site scores
# scores_df <- as.data.frame(scores(cap_W2, display = "sites"))
# 
# # add labels
# scores_df$Outcome <- meta_W2$Outcome
# 
# centroids <- scores_df %>%
#   group_by(Outcome) %>%
#   summarise(across(starts_with("CAP"), mean))
# 
# threshold <- mean(centroids$CAP1)
# 
# pred_labels <- ifelse(scores_df$CAP1 > threshold, "Relapse", "Cure")

# simple nearest-centroid classification
# library(FNN)
# # Use only CAP1
# centroid_mat <- centroids[, "CAP1", drop = FALSE]
# score_mat <- scores_df[, "CAP1", drop = FALSE]
# pred <- get.knnx(data = centroid_mat, query = score_mat, k = 1)
# pred_labels <- centroids$Outcome[pred$nn.index]


###########################################################
########## W2 Log2(TPM) CORRELATION - PCoA FIG ############

# Pull out the distances of each sample
scores_W2 <- as.data.frame(bd_W2$vectors)
scores_W2$SampleID <- rownames(scores_W2)
scores_W2$Outcome <- meta_W2$Outcome

# % variation explained = proportion of distance variation captured
eig <- bd_W2$eig
var_exp <- eig / sum(eig)
percent_var <- round(var_exp * 100, 1)
percent_var[1:5]

# Determine where the centroids are 
centroids_W2 <- scores_W2 %>%
  dplyr::group_by(Outcome) %>%
  dplyr::summarise(PCoA1 = mean(PCoA1), 
                   PCoA2 = mean(PCoA2))

scores_W2 <- merge(scores_W2, centroids_W2, by = "Outcome", suffixes = c("", "_centroid"))

# Make hull df for polygon on plot
hulls <- scores_W2 %>%
  group_by(Outcome) %>%
  slice(chull(PCoA1, PCoA2))

# Make the figure
PermDisp_PCoA_fig1 <- scores_W2 %>% 
  ggplot(aes(x = PCoA1, y = PCoA2, fill = Outcome, shape = Outcome)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8) +
  geom_polygon(data = hulls, aes(color = Outcome, group = Outcome), 
               fill = NA, linewidth = 0.6, show.legend = F) +
  geom_segment(data = scores_W2, aes(xend = PCoA1_centroid, yend = PCoA2_centroid, color = "black"), alpha = 0.3) +
  geom_point(data = centroids_W2, aes(x = PCoA1, y = PCoA2, color = Outcome), 
             size = 4, shape = 4, stroke = 1.5, show.legend = F) +
  scale_fill_manual(values = c(`Cure` = "#0072B2", `Relapse` = "#bc5300")) +  
  scale_color_manual(values = c(`Cure` = "#0072B2", `Relapse` = "#bc5300")) +  
  scale_shape_manual(values = c(`Cure` = 21, `Relapse` = 21)) + 
  labs(title = "W2: Log2(TPM+1) → Pearson Corr → PERMDISP",
       subtitle = paste0("Avg dist to centroid: Cure=", round(bd_W2$group.distances[[1]], 3), "; Relapse=", round(bd_W2$group.distances[[2]], 3), "; p=", permutest_bd_W2$tab$`Pr(>F)`[1]),
       x = "PCoA1",
       y = "PCoA2") +
  my_plot_themes
PermDisp_PCoA_fig1
ggsave(PermDisp_PCoA_fig1,
       file = paste0("W2_Corr_PCoA_v2.png"),
       path = "Figures/PERMDISP",
       dpi = 600,
       width = 7, height = 5.8, units = "in")


###########################################################
##################### W2 TPM euclidean ####################

meta_W2 <- GoodSputum60_Metadata %>% filter(Week == "Week 2")

samples_W2 <- meta_W2$SampleID
TPM_W2 <- GoodSputum60_tpmf[, samples_W2]

# remove genes with zero variance
TPM_W2_f <- TPM_W2[apply(TPM_W2, 1, var) != 0, ]
# Now only 4025 genes

TPM_W2_t <- t(TPM_W2_f)

# Make a distance matrix
dist_W2 <- vegdist(TPM_W2_t, method = "euclidean")

# Run PERMDISP
bd_W2 <- betadisper(dist_W2, meta_W2$Outcome)
bd_W2
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W2, group = meta_W2$Outcome)
# 
# No. of Positive Eigenvalues: 12
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Cure Relapse 
# 19327   28147 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 12 eigenvalues)
# PCoA1     PCoA2     PCoA3     PCoA4     PCoA5     PCoA6     PCoA7     PCoA8 
# 2.957e+09 1.421e+09 9.221e+08 8.060e+08 5.787e+08 5.134e+08 4.084e+08 3.792e+08 

permutest(bd_W2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 239348770 239348770 3.1218    999  0.117
# Residuals 11 843383991  76671272  

boxplot(bd_W2, main = "W2 Distance to Centroid")

plot_df_W2 <- data.frame(
  distance = bd_W2$distances,
  Outcome = meta_W2$Outcome,
  Timepoint = "W2")

###########################################################
##################### W2 TPM CORRELATION ##################

meta_W2 <- GoodSputum60_Metadata %>% filter(Week == "Week 2")

samples_W2 <- meta_W2$SampleID
TPM_W2 <- GoodSputum60_tpmf[, samples_W2]

# remove genes with zero variance
TPM_W2_f <- TPM_W2[apply(TPM_W2, 1, var) != 0, ]
# Now only 4025 genes

# Make a correlation distance matrix
cor_mat_W2 <- cor(TPM_W2_f, method = "pearson")
dist_W2 <- as.dist(1 - cor_mat_W2)

# Run PERMDISP
bd_W2 <- betadisper(dist_W2, meta_W2$Outcome, type = "centroid", add = TRUE)
bd_W2
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_W2, group = meta_W2$Outcome, type = "centroid", add = TRUE)
# 
# No. of Positive Eigenvalues: 11
# No. of Negative Eigenvalues: 0
# 
# Average distance to centroid:
#   Cure Relapse 
# 0.2380  0.4931 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 11 eigenvalues)
# PCoA1   PCoA2   PCoA3   PCoA4   PCoA5   PCoA6   PCoA7   PCoA8 
# 0.73637 0.38647 0.33338 0.18069 0.14243 0.12274 0.09213 0.05871 

permutest(bd_W2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.20013 0.200126 21.352    999  0.001 ***
#   Residuals 11 0.10310 0.009373                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1        

boxplot(bd_W2, main = "W2 Distance to Centroid")
plot(bd_W2)

anova(bd_W2)
# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     1 0.20013 0.200126  21.352 0.0007394 ***
#   Residuals 11 0.10310 0.009373                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# PERMANOVA
adonis2(dist_W2 ~ meta_W2$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_W2 ~ meta_W2$Outcome)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1  0.33941 0.16189 2.1247  0.013 *
#   Residual 11  1.75719 0.83811                
# Total    12  2.09660 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Of distances?
adonis2(dist(bd_W2$distances) ~ meta_W2$Outcome)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist(bd_W2$distances) ~ meta_W2$Outcome)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     1  0.20013 0.65999 21.352  0.001 ***
#   Residual 11  0.10310 0.34001                  
# Total    12  0.30323 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot_df_W2 <- data.frame(
  distance = bd_W2$distances,
  Outcome = meta_W2$Outcome,
  Timepoint = "W2")



###########################################################
################# ∆ Log2(TPM) euclidean ###################

# Adjust dataframe
delta_mat <- BothWk_log2TPMf_delta %>%
  dplyr::select(Gene, Patient, delta) %>%
  pivot_wider(names_from = Patient, values_from = delta) %>%
  as.data.frame()
rownames(delta_mat) <- delta_mat$Gene
delta_mat <- delta_mat[, -1]  # remove Gene column
delta_mat_f <- delta_mat[apply(delta_mat, 1, var) != 0, ] # Now 4025 genes
delta_mat_t <- t(delta_mat_f)

# Create metadata
meta_delta <- BothWk_log2TPMf_delta %>%
  dplyr::select(Patient, Outcome) %>%
  distinct()

# Make a distance matrix
dist_delta <- vegdist(delta_mat_t, method = "euclidean")
# dist_W2 <- vegdist(log2TPM_W2_t, method = "bray")

# Run PERMDISP
bd_delta <- betadisper(dist_delta, meta_delta$Outcome)
bd_delta
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_delta, group = meta_delta$Outcome)
# 
# No. of Positive Eigenvalues: 12
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Cure Relapse 
# 157.3   186.1 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 12 eigenvalues)
# PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
# 87343 61636 46936 41577 37685 33926 29767 23390 

permutest(bd_delta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1   2560  2559.9 1.7565    999  0.224
# Residuals 11  16032  1457.5      

boxplot(bd_delta, main = "W2 Distance to Centroid")
plot(bd_delta)



###########################################################
################# ∆ Log2(TPM) CORRELATION #################
# W0->W2 change BothWk_log2TPMf_delta

# Adjust dataframe
delta_mat <- BothWk_log2TPMf_delta %>%
  dplyr::select(Gene, Patient, delta) %>%
  pivot_wider(names_from = Patient, values_from = delta) %>%
  as.data.frame()
rownames(delta_mat) <- delta_mat$Gene
delta_mat <- delta_mat[, -1]  # remove Gene column
delta_mat_f <- delta_mat[apply(delta_mat, 1, var) != 0, ] # Now 4025 genes
delta_mat_t <- t(delta_mat_f)

# Create metadata
meta_delta <- BothWk_log2TPMf_delta %>%
  dplyr::select(Patient, Outcome) %>%
  distinct()

# Make a correlation distance matrix
cor_mat_delta <- cor(delta_mat_f, method = "pearson")
dist_delta <- as.dist(1 - cor_mat_delta)
# dist_delta <- as.dist(sqrt(2 * (1 - cor_mat_delta))) # If we do it like this we don't have to add=T

# Run PERMDISP
bd_delta <- betadisper(dist_delta, meta_delta$Outcome, type = "centroid", add = T)
bd_delta
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dist_delta, group = meta_delta$Outcome, type = "centroid", add =
#                    T)
# 
# No. of Positive Eigenvalues: 12
# No. of Negative Eigenvalues: 0
# 
# Average distance to centroid:
#   Cure Relapse 
# 0.4361  0.4882 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 12 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 0.4671 0.4014 0.3465 0.3196 0.2728 0.2605 0.2081 0.1914 

permutest(bd_delta)
permutest_bd_delta <- permutest(bd_delta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.008341 0.0083407 2.2042    999  0.176
# Residuals 11 0.041623 0.0037839  

boxplot(bd_delta, main = "W2 Distance to Centroid")
plot(bd_delta)









