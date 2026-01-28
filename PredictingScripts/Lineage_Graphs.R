# Make graphs showing distribution of lineages
# E. Lamont 
# 1/28/26

source("Import_data.R")

# Using Lineage_info from Import_SampleMetadata.R

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
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
  select(SUBJID, outcome, sample, Collection_TimePoint, main_lineage, sub_lineage) %>%
  mutate(Outcome = ifelse(outcome == "Prob Relap", "Relapse", outcome)) %>%
  mutate(Outcome = ifelse(Outcome == "success", "Cure", Outcome)) %>%
  mutate(main_lineage = recode(main_lineage, "lineage2" = "L2", "lineage3" = "L3", "lineage4" = "L4"))

# Remove duplicates by keeping just the D0 and S (don't know the difference)
my_lineages <- my_lineages %>%
  filter(Collection_TimePoint %in% c("D0", "S", "D0_2")) %>%
  distinct(SUBJID, .keep_all = T)

# Now there are not duplicates
my_lineages[duplicated(my_lineages$SUBJID),]

###########################################################
##################### MAKE A BARPLOT ######################

Lineage_plot <- my_lineages %>%
  filter(Outcome %in% c("Relapse", "Cure")) %>% 
  # mutate(Type = factor(Type)) %>% 
  ggplot(aes(x = Outcome, fill = main_lineage)) + 
  geom_bar(position = position_dodge(), color = "black") + 
  geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.3) + 
  scale_fill_manual(values = c("L2" = "dodgerblue2", "L3" = "purple3", "L4" = "red")) + 
  scale_y_continuous(limits = c(0,62), breaks = seq(0, 62, 10), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "# of samples", fill = "Lineage") + 
  my_plot_themes
Lineage_plot


