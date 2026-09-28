##### ALPHA DIVERSITY #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(stringr)
library(vegan)
library(ape)
library(picante)
library(parallel)
library(ggpubr)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "alpha_diversity")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

source(file.path("modules", "phylo_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

# Bacterial phylogeny (GTDB)
bac_tree <- read.tree(file.path(phydir, "phy_tree.tree"))

# Host phylogeny
host_consensus <- read.tree(file.path(indir, "host_consensus.tre"))

############################
#### RAREFACTION CURVES ####
############################

set.seed(1)

#### Filtered dataset ####

rare_results <- data.frame(Sample = character(),
                           S = numeric(),
                           se = numeric(),
                           subsample = numeric(),
                           stringsAsFactors = FALSE)

max <- 25000
step <- max / 20

for (s in seq(0, max, by=step)) {
  cat("Rarefying to", s, "sequences per sample...\n")
  rare <- rarefy(data.frame(phy_sp_f@otu_table), MARGIN = 2, sample = s, se = FALSE) %>% data.frame()
  colnames(rare) <- "S"
  rare <- rare %>% rownames_to_column(var = "Sample")
  rare$subsample <- s
  # Append to results
  rare_results <- rbind(rare_results, rare)
}

# Add metadata
metadata <- data.frame(Sample = sample_names(phy_sp_f),
                        lib_size = sample_sums(phy_sp_f),
                        Species = phy_sp_f@sam_data$Species,
                        Domestication = phy_sp_f@sam_data$Domestication,
                        stringsAsFactors = FALSE) %>%
            mutate(Species_short = str_replace(Species, "[a-z]+ ", ". ")) %>%
            arrange(Domestication, Species_short) %>%
            mutate(Species_short = factor(Species_short, levels = unique(Species_short)))

# Don't show subsamples larger than the achieved library size
rare_results_filt <- rare_results %>% left_join(metadata, by = "Sample") %>%
  filter(subsample <= lib_size)

write.csv(rare_results_filt, file.path(subdir, "rarefaction_results_filt.csv"), row.names = FALSE)

# Plot rarefaction curves
p <- ggplot(rare_results_filt) +
  geom_line(aes(x = subsample, y = S, group = Sample, colour = Species)) +
  scale_color_manual(values=species_palette, name = "Species") +
  facet_grid(~ Species_short) +
  theme(legend.position = "none", strip.text = element_text(face = "italic")) +
  xlab("Number of sequences sampled") +
  ylab("Observed species richness") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "rarefaction_curves_filt.png"), p, width=8, height=6)

#### Raw dataset ####

set.seed(1)

rare_results <- data.frame(Sample = character(),
                           S = numeric(),
                           se = numeric(),
                           subsample = numeric(),
                           stringsAsFactors = FALSE)

max <- 10^floor(log10(mean(sample_sums(phy_sp))))
step <- max / 20

for (s in seq(0, max, by=step)) {
  cat("Rarefying to", s, "sequences per sample...\n")
  rare <- rarefy(data.frame(phy_sp@otu_table), MARGIN = 2, sample = s, se = FALSE) %>% data.frame()
  colnames(rare) <- "S"
  rare <- rare %>% rownames_to_column(var = "Sample")
  rare$subsample <- s
  # Append to results
  rare_results <- rbind(rare_results, rare)
}

# Add metadata
metadata <- data.frame(Sample = sample_names(phy_sp),
                        lib_size = sample_sums(phy_sp),
                        Species = phy_sp@sam_data$Species,
                        Domestication = phy_sp@sam_data$Domestication,
                        stringsAsFactors = FALSE) %>%
            mutate(Species_short = str_replace(Species, "[a-z]+ ", ". ")) %>%
            arrange(Domestication, Species_short) %>%
            mutate(Species_short = factor(Species_short, levels = unique(Species_short)))

# Don't show subsamples larger than the achieved library size
rare_results_filt <- rare_results %>% left_join(metadata, by = "Sample") %>%
  filter(subsample <= lib_size)

write.csv(rare_results_filt, file.path(subdir, "rarefaction_results_raw.csv"), row.names = FALSE)

# Plot rarefaction curves
p <- ggplot(rare_results_filt) +
  geom_line(aes(x = subsample, y = S, group = Sample, colour = Species)) +
  scale_color_manual(values=species_palette, name = "Species") +
  facet_grid(~ Species_short) +
  theme(legend.position = "none", strip.text = element_text(face = "italic")) +
  xlab("Number of sequences sampled") +
  ylab("Observed species richness") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "rarefaction_curves_raw.png"), p, width=8, height=6)

#########################
#### ALPHA DIVERSITY ####
#########################

set.seed(1)

# Rarefy to the depth suggested by rarefaction curves
rar_level <- 20000
phy_sp_rarefied <- rarefy_even_depth(subset_samples(phy_sp_f, sample_sums(phy_sp_f) > rar_level), sample.size = rar_level, rngseed = 1)

alpha_div <- data.frame(estimate_richness(phy_sp_f, measures = c("Observed"))) %>%
  rownames_to_column(var = "Sample") %>% rename(filt = Observed) %>%
  left_join(data.frame(estimate_richness(phy_sp_rarefied, measures = c("Observed"))) %>%
              rownames_to_column(var = "Sample") %>% rename(filt_rarefied = Observed),
            by = "Sample")

# Also raw data
rar_level <- 500000
phy_sp_raw_rarefied <- rarefy_even_depth(subset_samples(phy_sp, sample_sums(phy_sp) > rar_level), sample.size = rar_level, rngseed = 1)

alpha_div <- left_join(alpha_div,
                      data.frame(estimate_richness(phy_sp, measures = c("Observed"))) %>%
                      rownames_to_column(var = "Sample") %>% rename(raw = Observed) %>%
              left_join(data.frame(estimate_richness(phy_sp_raw_rarefied, measures = c("Observed"))) %>%
                      rownames_to_column(var = "Sample") %>% rename(raw_rarefied = Observed),
                      by = "Sample"))

# Add metadata
alpha_div <- alpha_div %>%
  left_join(data.frame(phy_sp_f@sam_data) %>%
              rownames_to_column(var = "Sample") %>%
              select(Sample, Species, Genus, Group, Common.name, Domestication, contig_reads_count),
            by = c("Sample")) %>%
            mutate(Group = str_to_lower(Group)) %>%
            arrange(Domestication, Species) %>%
            mutate(Group = factor(Group, levels = unique(Group)))

write.csv(alpha_div, file = file.path(subdir, "alpha_diversity.csv"), quote = FALSE, row.names = FALSE)

# Filtered

p <- ggplot(alpha_div, aes(x=Group, y=filt)) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  facet_grid(Genus ~ ., scales = "free_y", space = "free_y") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Observed species richness") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_filt.png"), p, width=4, height=5)

# Filtered & Rarefied
p <- ggplot(alpha_div, aes(x=Group, y=filt_rarefied)) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  facet_grid(Genus ~ ., scales = "free_y", space = "free_y") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Observed species richness\n(after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_filt_rarefied.png"), p, width=4, height=5)

# Raw
p <- ggplot(alpha_div, aes(x=Group, y=raw)) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  facet_grid(Genus ~ ., scales = "free_y", space = "free_y") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Observed species richness") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_raw.png"), p, width=4, height=5)

# Raw & Rarefied
p <- ggplot(alpha_div, aes(x=Group, y=raw_rarefied)) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  facet_grid(Genus ~ ., scales = "free_y", space = "free_y") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Observed species richness\n(after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_raw_rarefied.png"), p, width=4, height=5)

####################
#### FAITH'S PD ####
####################

bac_tree$tip.label <- gsub("_", " ", bac_tree$tip.label)

# Calculate Faith's Phylogenetic Diversity
otu_table <- t(as.matrix(subset_taxa(phy_sp_rarefied, superkingdom == "Bacteria")@otu_table))

# Get phylogenetic diversity
phy_div <- pd(otu_table, bac_tree) %>% rownames_to_column(var = "Sample") %>%
  left_join(data.frame(phy_sp_f@sam_data) %>%
              rownames_to_column(var = "Sample") %>%
              select(Sample, Species, Genus, Group, Common.name, Domestication),
            by = c("Sample"))

write.csv(phy_div, file = file.path(subdir, "phylogenetic_diversity.csv"), quote = FALSE, row.names = FALSE)

# Plot phylogenetic diversity
p <- ggplot(phy_div, aes(x=Group, y=PD)) +
  geom_boxplot(aes(fill=Species)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=species_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp@sam_data$Group, phy_sp@sam_data$Species)) +
  facet_grid(Genus ~ ., scales = "free_y", space = "free_y") +
  theme(legend.position = "none", axis.title.y = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Faith's PD") +
  coord_flip()

ggsave(file.path(subdir, "phylogenetic_diversity.png"), p, width=4, height=5)

# Plot relationship between PD and species richness
p <- ggplot(phy_div, aes(x=SR, y=PD)) +
  geom_point(aes(colour = Species)) +
  #geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_colour_manual(values=species_palette, name = "Species") +
  #scale_x_discrete(labels = setNames(phy_sp@sam_data$Common.name, phy_sp@sam_data$Species)) +
  theme(legend.position = "bottom", axis.text = element_text(size = 8),
        strip.text = element_text(size = 10)) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  xlab("Observed species richness\n(after rarefaction)") + ylab("Faith's PD\n(after rarefaction)")

ggsave(file.path(subdir, "alpha_vs_pd.png"), width=4, height=5)

###################
#### RUN TESTS ####
###################

# Combine with alpha diversity
div <- full_join(select(alpha_div, c(Sample, filt, filt_rarefied, raw_rarefied, contig_reads_count)),
                  select(phy_div, c(-SR)), by = "Sample")

write.csv(div, file = file.path(subdir, "diversity.csv"), quote = FALSE, row.names = FALSE)

div_filt <- div %>% filter(!is.na(filt_rarefied) & !is.na(PD))

cor.test(div_filt$filt_rarefied, div_filt$PD, method = "pearson")

# Linear model
model <- aov(filt_rarefied ~ contig_reads_count + Genus*Domestication, data = div_filt)
res <- summary(model)[[1]]

shapiro.test(residuals(model))
summary(residuals(model))

write.csv(res, file = file.path(subdir, "anova_alpha_diversity_filt.csv"), quote = FALSE)

model_signif <- aov(filt_rarefied ~ Genus + Genus:Domestication, data = div_filt)

# Tukey's HSD
tukey <- do.call("rbind", TukeyHSD(model_signif)) %>% data.frame %>% filter(!is.na(diff))

# For the interactions, keep only those comparing domestication within genera
# As well as comparisons with humans
tukey <-
  tukey %>%
  rownames_to_column(var = "Comparison") %>%
  filter(str_detect(Comparison, "[A-Z][a-z]+-[A-Z][a-z]+$") |
        str_detect(Comparison, "Ovis:.*-Ovis:.*") |
        str_detect(Comparison, "Sus:.*-Sus:.*") |
        str_detect(Comparison, "Equus:.*-Equus:.*") |
        str_detect(Comparison, "Homo"))

write.csv(tukey, file = file.path(subdir, "tukey_hsd_alpha_diversity_filt.csv"), quote = FALSE, row.names = FALSE)