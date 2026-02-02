# Make graphs showing distribution of lineages
# E. Lamont 
# 1/28/26

# 2/2/26: re-doing after getting the lineage info from all the samples

source("Import_data.R")

# Using Lineage_info_all from Import_SampleMetadata.R

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10))# , 
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank())


###########################################################
###################### ORGANIZE DATA ######################

# Not sure what all the timepoints mean, what does S mean?

# Some of the samples are duplicated (multiple timepoints). Can see that here
Lineage_info[duplicated(Lineage_info$SUBJID),]

my_lineages <- metadata_4 %>% 
  select(SUBJID, outcome, sample, Collection_TimePoint, main_lineage, sub_lineage, arm) %>%
  mutate(Outcome = ifelse(outcome == "Prob Relap", "Relapse", outcome)) %>%
  mutate(Outcome = ifelse(Outcome == "success", "Cure", Outcome)) %>%
  mutate(main_lineage2 = recode(main_lineage, "lineage2" = "L2", "lineage3" = "L3", "lineage4" = "L4", "lineage1" = "L1")) %>%
  mutate(main_lineage2 = ifelse(grepl(";", main_lineage2), "mixed", main_lineage2)) %>% # recode the mixed lineages
  filter(main_lineage2 != "NA") %>%
  filter(!is.na(arm)) %>% # Remove the NA sample 14004 (where did this come from?)
  filter(Collection_TimePoint %in% c("D0", "S", "D0_2")) %>% # Remove duplicates by keeping just the D0 and S (don't know the difference)
  distinct(SUBJID, .keep_all = T) # Make sure am only keeping one SUBJID row

# Now there are not duplicates
my_lineages[duplicated(my_lineages$SUBJID),]

# See the arm fall out
my_lineages %>% 
  filter(Outcome %in% c("Relapse", "Cure")) %>% 
  group_by(arm, Outcome) %>%
  summarize(N_samples = n())

###########################################################
##################### MAKE A BARPLOT ######################

Lineage_plot <- my_lineages %>%
  filter(Outcome %in% c("Relapse", "Cure")) %>% 
  # mutate(Type = factor(Type)) %>% 
  ggplot(aes(x = Outcome, fill = main_lineage2)) + 
  geom_bar(position = position_dodge(), color = "black") + 
  geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.3) + 
  scale_fill_manual(values = c("L2" = "dodgerblue2", "L3" = "purple3", "L4" = "red", "L1" = "magenta2", "mixed" = "grey")) + 
  # scale_y_continuous(limits = c(0,62), breaks = seq(0, 62, 10), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "# of samples", fill = "Lineage",
       title = "subMIC subset lineages") + 
  my_plot_themes
Lineage_plot
# ggsave(Lineage_plot,
#        file = "subMIC_Subset_Lineages_Bar1.pdf",
#        path = "Figures/Lineages",
#        width = 4, height = 4, units = "in")


###########################################################
################# MAKE A STACKED BARPLOT ##################

my_lineages_Percents <- my_lineages %>%
  filter(Outcome %in% c("Relapse", "Cure")) %>%
  count(Outcome, main_lineage2) %>%
  group_by(Outcome) %>%
  mutate(Percent = round(n/sum(n) * 100, 1)) %>%
  ungroup()
my_lineages_Percents

Lineage_plot2 <- my_lineages_Percents %>%
  ggplot(aes(x = Outcome, y = Percent, fill = main_lineage2)) + 
  geom_col(color = "black") + 
  # geom_text(aes(label = scales::percent(after_stat(prop), accuracy = 1)), stat = "count", position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("L2" = "dodgerblue2", "L3" = "purple3", "L4" = "red", "L1" = "magenta2", "mixed" = "grey")) + 
  scale_y_continuous(limits = c(0,100.1), breaks = seq(0, 100.1, 10), expand = expansion(mult = c(0, 0))) +
  annotate("text", x = "Cure", y = 80, label = "41.0% L2", size = 3) +
  annotate("text", x = "Cure", y = 30, label = "55.9% L4", size = 3) +
  annotate("text", x = "Relapse", y = 80, label = "52.4% L2", size = 3) +
  annotate("text", x = "Relapse", y = 30, label = "42.9% L4", size = 3) +
  labs(x = NULL, y = "% of samples", fill = "Lineage",
       title = NULL) + 
  my_plot_themes
Lineage_plot2
ggsave(Lineage_plot2,
       file = "All_Lineages_Bar1.pdf",
       path = "Figures/Lineages",
       width = 4, height = 4, units = "in")

###########################################################
#################### CHI-SQUARED TEST #####################

my_lineages_testing <- my_lineages %>%
  filter(Outcome %in% c("Relapse", "Cure")) %>%
  filter(main_lineage2 %in% c("L2", "L4")) %>% # Just keep L2 and L4 to make it easy
  select(SUBJID, Outcome, main_lineage2)
  # count(Outcome, main_lineage2) # %>%
  # group_by(Outcome)

# Are lineage and outcome independent?
# What are the expected counts if they are independent?

# Create a table of the data
my_lineages_table <- table(my_lineages_testing$main_lineage2, my_lineages_testing$Outcome)

# Checked by hand that all expected counts are >5

# Run the Chi-square test
chisq.test(my_lineages_testing$main_lineage2, my_lineages_testing$Outcome, correct = F)
# X-squared = 1.2283, df = 1, p-value = 0.2677

# Run a more complicated Chi-square test that shows expected counts
library(gmodels)
CrossTable(my_lineages_testing$main_lineage2, my_lineages_testing$Outcome,
           expected = F, prop.r = T, prop.c = T, prop.chisq = F, chisq = T, format = "SPSS")


###########################################################
################ TWO-SAMPLE BINOMIAL TEST #################

# Same as Chi-squared but gives confidence interval as well

# Create a table of the data
my_lineages_table <- table(my_lineages_testing$main_lineage2, my_lineages_testing$Outcome)

prop.test(my_lineages_table, correct = F)
# X-squared = 1.2283, df = 1, p-value = 0.2677
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.09111547  0.02645518
# sample estimates:
#   prop 1    prop 2 
# 0.9147287 0.9470588 






