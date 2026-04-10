# GSEA
# E. Lamont
# 4/6/26

# Lets just look at W2 cure vs relapse for now 


source("Import_data.R")


###########################################################
################ GET DEG DATA (DUFFYTOOLS) ################

# Remove NAs 
res <- list_dfs_f2$W2.cure_vs_W2.relapse %>%
  # filter(!is.na(padj)) %>% 
  arrange(desc(LOG2FOLD))

# Create a ranked vector (using stat or Log2FoldChange)
ranks <- res$LOG2FOLD
names(ranks) <- res$GENE_ID
ranks <- sort(ranks, decreasing = TRUE)

# Filter out NAs or duplicated genes
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[!duplicated(names(ranks))]


###########################################################
################ CUSTOM GENE SETS ################

# Load custom gene sets. Needs to have columns Gene, GeneSet. Will take from allGeneSetList
custom_list <- allGeneSetList$EllaGeneSets_2026.03.19
EllaGeneSets_2026.03.19 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2026.03.19.csv") # Also need the csv to get the grouping for plotting with facets

# Run fgsea
fgsea_res <- fgsea(
  pathways = custom_list,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

###########################################################
####################### BUBBLE PLOT #######################

my_plot_themes <- theme_bw() +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10)
  )
facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))




# Add groupings
fgsea_res2 <- fgsea_res %>%
  mutate(pathway = str_wrap(pathway, width = 50)) %>% 
  mutate(FDR_Significance = ifelse(padj < 0.05, "significant", "not significant")) %>%
  left_join(EllaGeneSets_2026.03.19 %>% # Add the Group names
              rename(pathway = GeneSet) %>%
              dplyr::select(pathway, Group), by = "pathway")

# Make the bubble plot
my_bubblePlot <- fgsea_res2 %>%
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>%
  mutate(Group_wrapped = case_when(Group_wrapped == "Ribosomal proteins" ~ "Ribosomal\nproteins", Group_wrapped == "Hypoxia related" ~ "Hypoxia\nrelated", TRUE ~ Group_wrapped)) %>%
  filter(!Group %in% c("Toxin/Antitoxin", "ESX genes", "Metal", "Nucleic Acid")) %>% # Remove this because I don't think its interesting
  mutate(pathway2 = paste0(pathway, " (n=", size, ")")) %>%
  ggplot(aes(x = NES, y = pathway2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & NES>0 ~ "pos",
                                  FDR_Significance == "significant" & NES<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.9) + 
  scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
  facet_grid(rows = vars(Group_wrapped), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "W2 Relapse vs Cure (Run1-4) DuffyTools DEG -> GSEA", y = NULL, x = "Normalized Enrichment Score") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot




###########################################################
################ TRY TO GET GENE ONTOLOGY #################

# Found Gene Ontology download here:
# https://geneontology.org/docs/download-go-annotations/#2-all-other-organisms
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000277735.2/
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/277/735/GCF_000277735.2_ASM27773v2/


library(readr)

# gaf <- read_tsv("Data/GCF_000277735.2_ASM27773v2_gene_ontology.gaf", comment = "!")


gaf <- readr::read_tsv("Data/GCF_000277735.2_ASM27773v2_gene_ontology.gaf",
  comment = "!", col_names = FALSE, show_col_types = FALSE)

colnames(gaf) <- c("DB", "DB_Object_ID", "DB_Object_Symbol",
  "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", "With_or_From", "Aspect",
  "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID")

mtb <- gaf %>%
  filter(Taxon == "taxon:83332") %>%
  select(DB_Object_Symbol, GO_ID)

# Convert to gene sets
gaf_gene_sets <- gaf %>%
  group_by(GO_ID) %>%
  summarise(genes = list(unique(DB_Object_Symbol)))


###########################################################
################ NEW DOWNLOAD: TRY TO GET GENE ONTOLOGY #################

### Another thing to download here ###
# https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/

gaf_2 <- read_tsv("Data/30.M_tuberculosis_ATCC_25618.goa.txt",
  comment = "!", col_names = FALSE, show_col_types = FALSE)
colnames(gaf_2) <- c("DB", "DB_Object_ID", "GENE_NAME",
                   "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", "With_or_From", "Aspect",
                   "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID")

gaf_clean <- gaf_2 %>%
  filter(
    Taxon == "taxon:83332",  # Mycobacterium tuberculosis
    Qualifier != "NOT"       # remove negated annotations
  )

gaf_clean <- gaf_clean %>%
  filter(Evidence_Code != "IEA") # Remove electronic annotations

gaf_clean <- gaf_clean %>%
  mutate(GENE_ID = str_extract(DB_Object_Synonym, "Rv\\d+\\w*"))


gaf_GeneSet <- gaf_clean %>%
  filter(!is.na(GENE_ID)) %>%
  select(GO_ID, GENE_ID) %>%
  distinct() %>%
  group_by(GO_ID) %>%
  summarise(genes = list(GENE_ID), .groups = "drop")

gaf_GeneSetList <- setNames(gaf_GeneSet$genes, gaf_GeneSet$GO_ID)

# Run fgsea
fgsea_res <- fgsea(
  pathways = gaf_GeneSetList,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)
