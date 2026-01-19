##### EXPLORATORY PLOTS #####

#### Plots exploratory plots such as ordinations and heatmaps

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(microViz)
library(microbiome)
library(rphylopic)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(ggExtra)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
datadir <- normalizePath(file.path("..", "output", "function", "data"))
pathdir <- normalizePath(file.path("..", "output", "function", "pathway_completeness")) # Directory with pathway analysis output
subdir <- normalizePath(file.path("..", "output", "function", "exploratory_plots")) # subdirectory for the output of this script

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

# Phyloseq objects
phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(datadir, "phy_gene_f_clr.RDS"))

phy_pathway <- readRDS(file.path(pathdir, "phy_pathway.RDS"))
phy_pathway_clr <- readRDS(file.path(pathdir, "phy_pathway_clr.RDS"))

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#########################
#####  BETA DISPER  #####
#########################

# Calculate beta disper for species, order and diet, and run tukey's test

#### GENE COMPOSITION (CLR) ####

tukey_results <- data.frame()

# Compare species
disp <- betadisper(vegdist(t(otu_table(phy_gene_f_clr)), method = "euclidean"), group = phy_gene_f_clr@sam_data$Group)

disp_tukey <- TukeyHSD(disp, which = "group", ordered = FALSE)$group %>% data.frame %>% rownames_to_column("Comparison") %>%
  separate(Comparison, into = c("Group1", "Group2"), sep = "-")

tukey_results <- rbind(tukey_results, disp_tukey)

write.csv(tukey_results, file = file.path(subdir, "betadisper_gene_tukey_results.csv"), row.names = FALSE, quote = FALSE)

disp_df <- data.frame(Sample = sample_names(phy_gene_f_clr),
                      Group = phy_gene_f_clr@sam_data$Group,
                      Species = phy_gene_f_clr@sam_data$Species,
                      Genus = phy_gene_f_clr@sam_data$Genus,
                      Distance = disp$distances)

p <- ggplot(data = disp_df, aes(x = Group, y = Distance, fill = Species)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Species)) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    scale_fill_manual(values = species_palette, name = "Species") +
    theme(legend.position = "none", axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "betadisper_gene.png"), p, width = 8, height = 12)

#### PATHWAY COMPOSITION (CLR) ####

tukey_results <- data.frame()

# Compare species
disp <- betadisper(vegdist(t(otu_table(phy_pathway_clr)), method = "euclidean"), group = phy_pathway_clr@sam_data$Group)

disp_tukey <- TukeyHSD(disp, which = "group", ordered = FALSE)$group %>% data.frame %>% rownames_to_column("Comparison") %>%
  separate(Comparison, into = c("Group1", "Group2"), sep = "-")

tukey_results <- rbind(tukey_results, disp_tukey)

write.csv(tukey_results, file = file.path(subdir, "betadisper_path_tukey_results.csv"), row.names = FALSE, quote = FALSE)

disp_df <- data.frame(Sample = sample_names(phy_pathway_clr),
                      Group = phy_pathway_clr@sam_data$Group,
                      Species = phy_pathway_clr@sam_data$Species,
                      Genus = phy_pathway_clr@sam_data$Genus,
                      Distance = disp$distances)

p <- ggplot(data = disp_df, aes(x = Group, y = Distance, fill = Species)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Species)) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    scale_fill_manual(values = species_palette, name = "Species") +
    theme(legend.position = "none", axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "betadisper_path.png"), p, width = 8, height = 12)

###################
#### RDA PLOTS ####
###################

#### GENE ABUNDANCE (CLR) ####

# Recode order and habitat as TRUE and FALSE
phy_gene_f_clr <- phy_gene_f_clr %>%
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
ord <- ord_calc(phy_gene_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_genes.png"), p, width=3, height=3)

# Plot ordination

dom_shape_palette <- c("domestic" = 1, "wild" = 16, "feral" = 6, "human" = 8)

p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
  scale_color_manual(values=species_palette, name = "Species") +
  geom_phylopic(data = centroids(ord@ord, phy_gene_f), aes(colour = Species), uuid = centroids(ord@ord, phy_gene_f)$uid, width = 0.3, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))
  
p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "gene_ordination.png"), width=6, height=6)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get gene category
arrows$category <- as.character(phy_gene_f_clr@tax_table[match(rownames(arrows),  rownames(phy_gene_f_clr@tax_table)), "category"])

arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), 500))

# Save
write.csv(rownames_to_column(arrows, "gene"), file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2")), category)

# Group uncommon categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(6) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows_filt$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(arrows_filt$category_grouped)[-length(unique(arrows_filt$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey90" # Set "Other" to grey

p <- ggplot(data = arrows_filt) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Module") +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

ggsave(p, filename = file.path(subdir, "gene_ordination_arrows.png"), width=6, height=6)

#### PATH ABUNDANCE (CLR) ####

phy_pathway_clr <- phy_pathway_clr %>% 
         ps_mutate(Domestic_sheep = (Genus == "Ovis" & Domestication == "domestic"),
                  Feral_sheep = (Genus == "Ovis" & Domestication == "feral"),
                  Wild_argali = (Genus == "Ovis" & Domestication == "wild"),
                  Domestic_horse = (Genus == "Equus" & Domestication == "domestic"),
                  Wild_zebra = (Genus == "Equus" & Domestication == "wild"),
                  Domestic_pig = (Genus == "Sus" & Domestication == "domestic"),
                  Wild_boar = (Genus == "Sus" & Domestication == "wild"))

# Ordinate using all data
ord <- ord_calc(phy_pathway_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_pathways.png"), p, width=3, height=3)

# Color by order
p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
  scale_color_manual(values=species_palette, name = "Species") +
  geom_phylopic(data = centroids(ord@ord, phy_pathway_clr), aes(colour = Species), uuid = centroids(ord@ord, phy_pathway_clr)$uid, width = 0.2, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "pathway_ordination.png"), width=6, height=6)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get gene category
arrows$name <- as.character(phy_pathway_clr@tax_table[match(rownames(arrows), rownames(phy_pathway_clr@tax_table)), "path_name"])
arrows$category <- as.character(phy_pathway_clr@tax_table[match(rownames(arrows), rownames(phy_pathway_clr@tax_table)), "path_class"])
arrows$category <- str_remove(arrows$category, ".*; ")

arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), nrow(arrows)))

# Save
write.csv(rownames_to_column(arrows, "gene"), file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2")), name, category)

# Group uncommon categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(8) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows_filt$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(arrows_filt$category_grouped)[-length(unique(arrows_filt$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey90" # Set "Other" to grey

p <- ggplot(data = arrows_filt) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  #geom_text(aes(x = RDA1, y = RDA2, label = name, hjust = ifelse(RDA1 < 0, 1, 0)), size = 1.5, vjust = 1) +
  scale_color_manual(values = arrow_colours, name = "Pathway class") +
  xlim(c(min(arrows_filt$RDA1)*2, max(arrows_filt$RDA1)*2)) +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(size = 7), legend.title = element_text(size = 10)) +
  guides(colour = guide_legend(ncol = 2))

ggsave(p, filename = file.path(subdir, "pathway_ordination_arrows.png"), width=6, height=5)

###################
#### PERMANOVA ####
###################

# Remove humans for permanova analyses

#### PATHWAYS ####

# Explanatory variables
otus <- t(as.data.frame(subset_samples(phy_pathway_clr, Species != "Homo sapiens")@otu_table))

sample_data <- as.data.frame(subset_samples(phy_pathway_clr, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ genus * dom + reads,
        permutations = 1000, by = "term", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_pathways.csv"), row.names = TRUE, quote = TRUE)

### GENES ###

# Explanatory variables
otus <- t(as.data.frame(subset_samples(phy_gene_f_clr, Species != "Homo sapiens")@otu_table))

sample_data <- as.data.frame(subset_samples(phy_gene_f_clr, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ genus * dom + reads,
        permutations = 1000, by = "term", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_genes.csv"), row.names = TRUE, quote = TRUE)

############################
#### DISTANCES TO HUMAN ####
############################

#### GENE CLR ABUNDANCES ####

# Calculate human centroid

human_centroid <- rowMeans(otu_table(subset_samples(phy_gene_f_clr, Species == "Homo sapiens")))

# Calculate Aitchison distances of each sample to the human centroid
distances = data.frame(sample = character(), distance = numeric())

nonhuman_samples <- subset_samples(phy_gene_f_clr, Species != "Homo sapiens") %>% sample_names

for (sample in nonhuman_samples) {
  values <- otu_table(phy_gene_f_clr)[,sample]
  df <- cbind(values, human_centroid) %>% t
  dist <- vegdist(df, method = "euclidean")
  distances <- rbind(distances, data.frame(sample = sample, distance = as.numeric(dist)))
}

# Add metadata
sample_meta <- data.frame(phy_gene_f_clr@sam_data) %>%
      select(Species, Genus, Common.name, Domestication) %>% rownames_to_column("sample")

distances <- distances %>% left_join(sample_meta, by = "sample")

write.csv(distances, file = file.path(subdir, "distances_to_human_gene_clr.csv"), quote = FALSE, row.names = FALSE)

# Plot
p <- ggviolin(data = distances, x = "Domestication", y = "distance", fill = "Species", facet.by = "Genus") +
  scale_fill_manual(values = species_palette) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Aitchison distances")

# Run Kruskal Wallis tests
stat.test <- distances %>%
  group_by(Genus) %>%
  wilcox_test(distance ~ Domestication) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()

write.csv(stat.test, file = file.path(subdir, "distances_to_human_gene_clr_test.csv"), quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>% add_xy_position(x = "Domestication")

p <- p +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = -2, hide.ns = TRUE,
    label = "{p.adj.signif}")

ggsave(file.path(subdir, "distances_to_human_gene_clr.png"), p, width=5, height=4)

#### PATHWAY CLR ABUNDANCES ####

# Calculate human centroid

human_centroid <- rowMeans(otu_table(subset_samples(phy_pathway_clr, Species == "Homo sapiens")))

# Calculate Aitchison distances of each sample to the human centroid
distances = data.frame(sample = character(), distance = numeric())

nonhuman_samples <- subset_samples(phy_pathway_clr, Species != "Homo sapiens") %>% sample_names

for (sample in nonhuman_samples) {
  values <- otu_table(phy_pathway_clr)[,sample]
  df <- cbind(values, human_centroid) %>% t
  dist <- vegdist(df, method = "euclidean")
  distances <- rbind(distances, data.frame(sample = sample, distance = as.numeric(dist)))
}

# Add metadata
sample_meta <- data.frame(phy_pathway_clr@sam_data) %>%
      select(Species, Genus, Common.name, Domestication) %>% rownames_to_column("sample")

distances <- distances %>% left_join(sample_meta, by = "sample")

write.csv(distances, file = file.path(subdir, "distances_to_human_path_clr.csv"), quote = FALSE, row.names = FALSE)

# Plot
p <- ggviolin(data = distances, x = "Domestication", y = "distance", fill = "Species", facet.by = "Genus") +
  scale_fill_manual(values = species_palette) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Aitchison distances")

# Run Kruskal Wallis tests
stat.test <- distances %>%
  group_by(Genus) %>%
  wilcox_test(distance ~ Domestication) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()

write.csv(stat.test, file = file.path(subdir, "distances_to_human_path_clr_test.csv"), quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>% add_xy_position(x = "Domestication")

p <- p +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = -2, hide.ns = TRUE,
    label = "{p.adj.signif}")

ggsave(file.path(subdir, "distances_to_human_path_clr.png"), p, width=5, height=4)
