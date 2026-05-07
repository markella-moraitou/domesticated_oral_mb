##### COMPARE DISTANCE TO HUMAN MICROBIOM #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggnewscale)
library(ggExtra)
library(ggpubr)
library(rstatix)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output", "community_analysis")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "human_dist")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

# Get presence absence of taxa (0 and 1)
phy_sp_f_pa <- microbiome::transform(phy_sp_f, "pa")

############################
#### DISTANCES TO HUMAN ####
############################

#### CLR ABUNDANCES ####

# Calculate human centroid

human_centroid <- rowMeans(otu_table(subset_samples(phy_sp_f_clr, Species == "Homo sapiens")))

# Calculate Aitchison distances of each sample to the human centroid
distances = data.frame(sample = character(), distance = numeric())

nonhuman_samples <- subset_samples(phy_sp_f_clr, Species != "Homo sapiens") %>% sample_names

for (sample in nonhuman_samples) {
  values <- otu_table(phy_sp_f_clr)[,sample]
  df <- cbind(values, human_centroid) %>% t
  dist <- vegdist(df, method = "euclidean")
  distances <- rbind(distances, data.frame(sample = sample, distance = as.numeric(dist)))
}

# Add metadata
sample_meta <- data.frame(phy_sp_f_clr@sam_data) %>%
      select(Species, Genus, Common.name, Domestication) %>% rownames_to_column("sample")

distances <- distances %>% left_join(sample_meta, by = "sample")

write.csv(distances, file = file.path(subdir, "distances_to_human_clr.csv"), quote = FALSE, row.names = FALSE)

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

write.csv(stat.test, file = file.path(subdir, "distances_to_human_clr_test.csv"), quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>% add_xy_position(x = "Domestication")

p_clr <- p +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = -2, hide.ns = TRUE,
    label = "{p.adj.signif}")

ggsave(file.path(subdir, "distances_to_human_clr.png"), p_clr, width=5, height=4)

#### PhilR ####

human_centroid <- colMeans(otu_table(subset_samples(phy_sp_philr, Species == "Homo sapiens")))

# Calculate Jaccrd distances of each sample to the human centroid
distances = data.frame(sample = character(), distance = numeric())

nonhuman_samples <- subset_samples(phy_sp_philr, Species != "Homo sapiens") %>% sample_names

for (sample in nonhuman_samples) {
  values <- t(otu_table(phy_sp_philr))[,sample]
  df <- cbind(values, human_centroid) %>% t
  dist <- vegdist(df, method = "euclidean")
  distances <- rbind(distances, data.frame(sample = sample, distance = as.numeric(dist)))
}

# Add metadata
sample_meta <- data.frame(phy_sp_f_pa@sam_data) %>%
      select(Species, Genus, Common.name, Domestication) %>% rownames_to_column("sample")

distances <- distances %>% left_join(sample_meta, by = "sample")

write.csv(distances, file = file.path(subdir, "distances_to_human_philr.csv"), quote = FALSE, row.names = FALSE)

# Plot
p <- ggviolin(data = distances, x = "Domestication", y = "distance", fill = "Species", facet.by = "Genus") +
  scale_fill_manual(values = species_palette) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("PhILR distances")

# Run Kruskal Wallis tests
stat.test <- distances %>%
  group_by(Genus) %>%
  wilcox_test(distance ~ Domestication) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()

write.csv(stat.test, file = file.path(subdir, "distances_to_human_philr_test.csv"), quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>% add_xy_position(x = "Domestication")

p_philr <- p +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = 0, hide.ns = TRUE,
    label = "{p.adj.signif}")

ggsave(file.path(subdir, "distances_to_human_philr.png"), p_philr, width=5, height=4)

#### PRESENCE ABSENCE ####

human_centroid <- rowMeans(otu_table(subset_samples(phy_sp_f_pa, Species == "Homo sapiens")))

# Calculate Jaccrd distances of each sample to the human centroid
distances = data.frame(sample = character(), distance = numeric())

nonhuman_samples <- subset_samples(phy_sp_f_pa, Species != "Homo sapiens") %>% sample_names

for (sample in nonhuman_samples) {
  values <- otu_table(phy_sp_f_pa)[,sample]
  df <- cbind(values, human_centroid) %>% t
  dist <- vegdist(df, method = "jaccard")
  distances <- rbind(distances, data.frame(sample = sample, distance = as.numeric(dist)))
}

# Add metadata
sample_meta <- data.frame(phy_sp_f_pa@sam_data) %>%
      select(Species, Genus, Common.name, Domestication) %>% rownames_to_column("sample")

distances <- distances %>% left_join(sample_meta, by = "sample")

write.csv(distances, file = file.path(subdir, "distances_to_human_pa.csv"), quote = FALSE, row.names = FALSE)

# Plot
p <- ggviolin(data = distances, x = "Domestication", y = "distance", fill = "Species", facet.by = "Genus") +
  scale_fill_manual(values = species_palette) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Jaccard distances")

# Run Kruskal Wallis tests
stat.test <- distances %>%
  group_by(Genus) %>%
  wilcox_test(distance ~ Domestication) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()

write.csv(stat.test, file = file.path(subdir, "distances_to_human_pa_test.csv"), quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>% add_xy_position(x = "Domestication")

p_pa <- p +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = 0, hide.ns = TRUE,
    label = "{p.adj.signif}")

ggsave(file.path(subdir, "distances_to_human_pa.png"), p_pa, width=5, height=4)

## Combine

p_combined <- plot_grid(p_pa + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
                        p_clr + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
                                        strip.background = element_blank(), strip.text = element_blank()),
                        p_philr + theme(strip.background = element_blank(), strip.text = element_blank()),
                        ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1, 1.2))

ggsave(file.path(subdir, "distances_to_human_combined.png"), p_combined, width=5, height=8)
