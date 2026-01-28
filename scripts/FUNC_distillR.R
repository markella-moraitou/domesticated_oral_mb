##### INPUT FOR FUNCTIONAL ANALYSIS #####

#### Prepares tables and phyloseq objects for analysis of functional annotations

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(ggplot2)
library(ggExtra)
library(microViz)
library(rphylopic)
library(RColorBrewer)
library(distillR)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
taxdir <- normalizePath(file.path("..", "output", "community_analysis")) # Directory with taxonomy analysis
datadir <- normalizePath(file.path("..", "output", "function", "data")) # Directory with functional data files
subdir <- normalizePath(file.path("..", "output", "function", "distillR")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))

# Stratified sample data
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

################################
#### DISTILL COMMUNITY WIDE ####
################################

# Identify functions in entire communities

if(file.exists(file.path(subdir, "GIFTs_community.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs_community.csv"), row.names = 1)
} else {
  # Prep input for distillR
  data <- psmelt(phy_gene_f) %>% select(OTU, Sample, gene_name, Abundance, Total_abundance) %>%
      # Remove zero abundances
      filter(Abundance > 0)
  cat("Running distillR to get GIFTs\n")
  GIFTs <- distill(data, GIFT_db, genomecol=2, annotcol=c(1, 3))
  write.csv(GIFTs, file.path(subdir, "GIFTs_community.csv"), row.names = TRUE)
}

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)
GIFTs_elements <- GIFTs_elements[, colSums(GIFTs_elements) > 0] # Remove empty columns

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_element", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Common.name, Group, Genus, Domestication))) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

p <- ggplot(GIFTs_elements_long, aes(y=Sample, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(Group), cols = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements_community.png"), width = 25, height = 20)

write.csv(GIFTs_elements_long, file.path(subdir, "GIFTs_elements_community.csv"), row.names = FALSE)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_functions <- GIFTs_functions[, colSums(GIFTs_functions) > 0] # Remove empty columns

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_function", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Common.name, Group, Genus, Domestication))) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain"))))

p <- ggplot(GIFTs_functions_long, aes(y=Sample, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(Group), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions_community.png"), width = 10, height = 15)

write.csv(GIFTs_functions_long, file.path(subdir, "GIFTs_functions_long.csv"), row.names = FALSE)

#############
#### RDA ####
#############

# Add GIFT elemens to a phyloseq and get host groups as separate variables
phy_distillr <- phyloseq(sample_data(phy_gene_f), otu_table(GIFTs_elements, taxa_are_rows = FALSE))

phy_distillr <- phy_distillr %>%
        ps_mutate(Domestic_sheep = (Genus == "Ovis" & Domestication == "domestic"),
                  Feral_sheep = (Genus == "Ovis" & Domestication == "feral"),
                  Wild_argali = (Genus == "Ovis" & Domestication == "wild"),
                  Domestic_horse = (Genus == "Equus" & Domestication == "domestic"),
                  Wild_zebra = (Genus == "Equus" & Domestication == "wild"),
                  Domestic_pig = (Genus == "Sus" & Domestication == "domestic"),
                  Wild_boar = (Genus == "Sus" & Domestication == "wild"))

# Species traits to use as constraints
species_traits <- c("Domestic_sheep", "Feral_sheep", "Wild_argali",
                    "Domestic_horse", "Wild_zebra",
                    "Domestic_pig", "Wild_boar")

# Ordinate using all data
ord <- ord_calc(phy_distillr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_gifts.png"), p, width=3, height=3)

# Plot ordination

dom_shape_palette <- c("domestic" = 1, "wild" = 16, "feral" = 6, "human" = 8)

p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
  scale_color_manual(values=species_palette, name = "Species") +
  geom_phylopic(data = centroids(ord@ord, phy_distillr), aes(colour = Species), uuid = centroids(ord@ord, phy_distillr)$uid, width = 0.1, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))
  
p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "gifts_ordination.png"), width=6, height=6)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get info for for plotting
arrows <- GIFT_db %>% select(Code_element, Element, Function, Domain) %>%
  mutate(Element = paste(Element, Domain)) %>% unique %>%
  right_join(rownames_to_column(data.frame(arrows), "Code_element"), by = "Code_element") %>%
  column_to_rownames("Code_element") %>%
  arrange(desc(distance))

# label 8 strongest associations plus the 5 largest effects of RDA1 and RDA2
arrows$to_label <- (rownames(arrows) %in% head(rownames(arrows), 8) | 
                    rownames(arrows) %in% head(rownames(arrange(arrows, desc(abs(RDA1)))), 5) |
                    rownames(arrows) %in% head(rownames(arrange(arrows, desc(abs(RDA2)))), 5))

# Save
write.csv(rownames_to_column(arrows, "Code_element"), file.path(subdir, "gift_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Group uncommon categories
common_categories <- table(arrows$Function) %>% sort(decreasing = TRUE) %>% head(6) %>% names

arrows <- arrows %>%
    mutate(category_grouped = factor(case_when(Function %in% common_categories ~ Function,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(arrows$category_grouped)[-length(unique(arrows$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey50" # Set "Other" to grey

p <- ggplot(data = arrows) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Module") +
  geom_label(aes(label = Element, x = RDA1, y = RDA2, colour = category_grouped, hjust = ifelse(RDA1 < 0, 0.2, 0.8)), size = 2, data = filter(arrows, to_label)) +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

ggsave(p, filename = file.path(subdir, "gifts_ordination_arrows.png"), width=6, height=6)

###################
#### PERMANOVA ####
###################

# Explanatory variables
otus <- as.data.frame(subset_samples(phy_distillr, Species != "Homo sapiens")@otu_table)

sample_data <- as.data.frame(subset_samples(phy_distillr, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ genus * dom + reads,
        permutations = 1000, by = "term", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_gifts.csv"), row.names = TRUE, quote = TRUE)

##########################
#### DISTILL BY TAXON ####
##########################

if(file.exists(file.path(subdir, "GIFTs_by_taxon.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs_by_taxon.csv"), row.names = 1)
} else {
  cat("Running distillR to get GIFTs\n")
  data <- gene_str %>% select(gene_id, database, species, Sample, mapped_reads) %>%
    # Keep entries with a totalAvgDepth > 1
    filter(mapped_reads > 1) %>%
    # Keep only taxa identified at the species level
    filter(species != "no support") %>%
    #filter(species %in% diff_taxa$species) %>%
    # Keep genes that are present in phy_gene_f and are KEGG
    filter(gene_id %in% taxa_names(phy_gene_f)) %>%
    filter(database == "KEGG") %>%
    # Combine species and sample into a "genome" column
    mutate(genome = paste(species, Sample, sep = " - ")) %>%
    # Remove genomes with less than 200 genes
    group_by(genome) %>% mutate(n_genes = n_distinct(gene_id)) %>% filter(n_genes > 400) %>%
    # Keep microbial species found in at least 5 samples
    group_by(species) %>% filter(n_distinct(Sample) > 5) %>% ungroup %>%
    select(genome, gene_id)
  GIFTs <- distill(data, GIFT_db, genomecol=1, annotcol=2)
  write.csv(GIFTs, file.path(subdir, "GIFTs_by_taxon.csv"), row.names = TRUE)
}

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)
GIFTs_elements <- GIFTs_elements[, colSums(GIFTs_elements) > 0] # Remove empty columns

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("genome") %>%
  separate(genome, into = c("species", "Sample.ID"), sep = " - ", remove = FALSE) %>%
  pivot_longer(cols = -c(species, Sample.ID, genome), names_to = "Code_element", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Sample.ID, Common.name, Group, Genus, Domestication))) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  mutate(genus = str_remove(species, " .*")) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

write.csv(GIFTs_elements_long, file.path(subdir, "GIFTs_elements_by_taxon.csv"), row.names = FALSE)

# Plot genera with at least 10 taxa
GIFTs_elements_long <- GIFTs_elements_long %>% group_by(genus) %>%
  filter(n_distinct(genome) > 10) %>% ungroup()

p <- ggplot(GIFTs_elements_long, aes(x=genome, y=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90)) +
  labs(y="Compound", x="Sample", fill="Completeness") +
  facet_grid(cols = vars(genus), rows = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements_by_taxon.png"), width = 25, height = 20)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_functions <- GIFTs_functions[, colSums(GIFTs_functions) > 0] # Remove empty columns

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("genome") %>%
  separate(genome, into = c("species", "Sample.ID"), sep = " - ", remove = FALSE) %>%
  pivot_longer(cols = -c(species, Sample.ID, genome), names_to = "Code_function", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Sample.ID, Species, Common.name, Group, Genus, Domestication))) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain")))) %>%
  mutate(genus = str_remove(species, " .*"))

p <- ggplot(GIFTs_functions_long, aes(y=genome, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(genus), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions_by_taxon.png"), width = 10, height = 15)

write.csv(GIFTs_functions_long, file.path(subdir, "GIFTs_functions_by_taxon.csv"), row.names = FALSE)

#### PCA ####

ord <- prcomp(GIFTs_elements)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings <- data.frame(Code_element = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")]) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  mutate(Variables = paste(Element, Domain, sep = " ")) %>%
  # Keep the longest arrows
  arrange(desc(sqrt(PC1^2 + PC2^2))) %>%
  slice(1:12) %>% select(PC1, PC2, Variables, Function)

metadata <- GIFTs_functions_long[match(rownames(ord$x), GIFTs_functions_long$genome),] %>%
  select(Sample, Species, Common.name, Group, Genus, Domestication, genome, genus) %>%
  # Group genomes to highlight most abundant genera in plot
  group_by(genus) %>%
  mutate(n_genomes = n_distinct(genome))

top_genera <- metadata %>% arrange(desc(n_genomes)) %>% pull(genus) %>% unique %>% head(8)

metadata <- metadata %>%
  mutate(genus_grouped = case_when(genus %in% top_genera ~ genus,
                                   TRUE ~ "Other")) %>%
  mutate(genus_grouped = factor(genus_grouped, levels = c(top_genera, "Other")))

# Get palette with RColorBrewer
genus_palette <- brewer.pal(n = length(unique(metadata$genus_grouped))-1, name = "Dark2")
names(genus_palette) <- unique(metadata$genus_grouped)[-which(unique(metadata$genus_grouped) == "Other")]
genus_palette["Other"] <- "grey"

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=metadata$genus_grouped), data = data.frame(ord$x)) +
  geom_point(size=2, alpha = 0.8) +
  scale_colour_manual(values = genus_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*6),
                                       yend = (PC2*6)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  geom_label(data = loadings, aes(x = (PC1*6), y = (PC2*6), label = Variables),
            size = 2, hjust = 0.5, vjust = -0.5, color = "black", alpha = 0.7) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(pca, filename = file.path(subdir, "PCA_GIFTs_by_taxon.png"), width = 10, height = 10)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=metadata$Species, shape = metadata$Domestication), data = data.frame(ord$x)) +
  geom_point(size=1, alpha = 0.8) +
  scale_colour_manual(values = species_palette, name = "") +
  scale_shape_manual(values = dom_shape_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  facet_wrap(~metadata$genus_grouped, ncol = 2) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3))

ggsave(pca, filename = file.path(subdir, "PCA_GIFTs_by_taxon_faceted.png"), width = 10, height = 10)
