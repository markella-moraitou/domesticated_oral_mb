##### NORMALISATIONS #####

#### Access OTU tree from Open Tree of Life, perform philr and clr normalisations
#### and create subsets for further analysis

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(tidyr)
library(tibble)
library(readr)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(microbiome)
library(stringr)
library(philr)
library(microbiome)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

## Set up for plotting
source(file.path("plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))

# Links to download GTDB tree
url <- "https://data.gtdb.ecogenomic.org/releases/release220/220.0/"

# Download tree and metadata
options(timeout = 600)
download.file(paste0(url, "bac120_r220.tree.gz"), file.path(outdir, "bac120_r220.tree.gz"))
download.file(paste0(url, "bac120_taxonomy_r220.tsv.gz"), file.path(outdir, "bac120_taxonomy_r220.tsv.gz"))

bac_tree <- read.tree(gzfile(file.path(outdir, "bac120_r220.tree.gz")))
bac_meta <- read_tsv(file.path(outdir, "bac120_taxonomy_r220.tsv.gz"), col_names = FALSE)

#######################
#### GET TAXA TREE ####
#######################

# Keep only metadata in tree
bac_meta_f <- bac_meta %>% filter(X1 %in% bac_tree$tip.label)

# Split taxonomy column
bac_meta_f <- bac_meta_f %>% separate(X2, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(s = str_remove(s, "s__"))

# Filter metadata to only include taxa in phyloseq object
bac_taxa <- data.frame(taxa_names = taxa_names(subset_taxa(phy_sp_f, superkingdom != "Archaea")))
bac_taxa$fix_names <- str_remove(bac_taxa$taxa_names, "\\*")

cat('Taxa that are not found in the bac120_taxonomy_r220 table:\n')
bac_taxa[!bac_taxa$fix_names %in% bac_meta_f$s]

# Extract tip labels from the species in the dataset
bac_meta_f <- filter(bac_meta_f, s %in% bac_taxa$fix_names) %>% select(X1, p, s) %>% ungroup()

# Add phyname
bac_meta_f$phyname <- bac_taxa$taxa_names[match(bac_meta_f$s, bac_taxa$fix_names)]

# Subset tree
tree <- drop.tip(bac_tree, setdiff(bac_tree$tip.label, bac_meta_f$X1))
tree$tip.label <- bac_meta_f$phyname[match(tree$tip.label, bac_meta_f$X1)]
tree$node.label <- paste0("N", 1:tree$Nnode) # Not very informative for now

write.tree(tree, file = file.path(phydir, "phy_tree.tree"))

#################################
#### CREATE PHYLOSEQ OBJECTS ####
#################################

# CLR transformation
phy_sp_f_clr <- phy_sp_f %>% microbiome::transform("clr")

# PhiLR transformation
phy_sp_philr <- phyloseq(otu_table(phy_sp_f), sample_data(phy_sp_f), tax_table(phy_sp_f), phy_tree(tree))
philr_otu <- philr(phy_sp_philr, pseudocount=10^-5)

write.table(philr_otu, file = file.path(phydir, "phy_sp_philr_OTU.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)

phy_sp_philr <- phyloseq(otu_table(philr_otu, taxa_are_rows = FALSE), sample_data(phy_sp_f))

# Save phyloseq objects
saveRDS(phy_sp_f, file.path(phydir, "phy_sp_f.RDS"))
saveRDS(phy_sp_f_clr, file.path(phydir, "phy_sp_f_clr.RDS"))
saveRDS(phy_sp_philr, file.path(phydir, "phy_sp_philr.RDS"))

###################
#### PLOT TREE ####
###################

# Plot microbial tree and colour tips by phylum
tree_meta <- bac_meta_f %>% filter(s %in% tree$tip.label) %>% rename(tip.label = s) %>%
  mutate(p = str_remove(p, "p__")) %>% select(tip.label, p)

p <- ggtree(tree)  %<+% tree_meta +
  geom_tiplab(size = 1, aes(colour = p), name = "Phylum") +
  scale_color_manual(values = phylum_palette) +
  ggtitle("Bacterial phylogenetic tree for OTUs in dataset") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size=3)))

# Add average abundance as barplot 
abund_df <- as.data.frame(otu_table(phy_sp_f)) %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance") %>%
  group_by(OTU) %>%
  summarise(Mean_Abundance = mean(Abundance)) %>%
  rename(tip.label = OTU)

p <- p + geom_fruit(data = abund_df, geom = geom_bar,
                    mapping = aes(x = Mean_Abundance, y = tip.label),
                    orientation = "y", stat = "identity",
                    pwidth = 0.3, size = 0.1, fill = "darkgrey")

ggsave(p, filename = file.path(outdir, "phy_tree.png"), width = 5, height = 30, units = "in", dpi = 300)

# Remove downloaded files to save space
file.remove(file.path(outdir, "bac120_r220.tree.gz"))
file.remove(file.path(outdir, "bac120_taxonomy_r220.tsv.gz"))
