##### Plot MAG tree #####

#### Plot MAG trees for Bacteria and Archaea and annotate using relevant metadata

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output")) # Directory with output
subdir <- normalizePath(file.path(outdir, "mags")) # subdirectory for the output of this script

source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(subdir, "bac_meta.tsv"), sep="\t", header=TRUE, quote = "", comment = "")
ar_meta <- read.table(file.path(subdir, "ar_meta.tsv"), sep="\t", header=TRUE, quote = "", comment = "")

# MAG trees
bac_tree <- read.tree(file = file.path(subdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(subdir, "ar_tree.tree"))

# List of dereplicated bins
drep_bins <- read.table(file.path(outdir, "M4_dereplicated_bins", "dereplicated_bins_list.txt"), header=FALSE, sep="\t") %>% pull(V1)

###########################
#### KEEP ONLY HQ MAGS ####
###########################

drep_bins <- gsub("dereplicated_genomes/", "", drep_bins)

## Keep only HQ MAGs
hq_bacs <- bac_meta %>% filter(Completeness >= 90 & Contamination <= 5) %>% filter(bin %in% drep_bins) %>% pull(label) 

hq_ars <- ar_meta %>% filter(Completeness >= 90 & Contamination <= 5) %>% filter(bin %in% drep_bins) %>% pull(label)

## Subset trees to only include HQ MAGs
bac_tree <- drop.tip(bac_tree, setdiff(bac_tree$tip.label, hq_bacs))
ar_tree <- drop.tip(ar_tree, setdiff(ar_tree$tip.label, hq_ars))

## Subset metadata to only include HQ MAGs
bac_meta <- bac_meta %>% filter(label %in% c(bac_tree$tip.label, bac_tree$node.label))
ar_meta <- ar_meta %>% filter(label %in% c(ar_tree$tip.label, ar_tree$node.label))

# Save trees and tables
write.tree(bac_tree, file = file.path(subdir, "bac_tree_drep.tree"))
write.tree(ar_tree, file = file.path(subdir, "ar_tree_drep.tree"))

write.table(bac_meta, file = file.path(subdir, "bac_meta_drep.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(ar_meta, file = file.path(subdir, "ar_meta_drep.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

########################
#### PROCESS TABLES ####
########################

# Get specimen ages and damage for plotting
ages_damage <- bac_meta %>% rbind(ar_meta) %>% select(label, bin, domain, Sample, host_genus, host_species, Period, min_damage_model_p, q1_damage_model_p, mean_damage_model_p, median_damage_model_p, q3_damage_model_p, max_damage_model_p) %>%
  filter(!is.na(bin))

####################
#### PLOT TREES ####
####################

#### Bacteria tree ####
# Colour by order and habitat
bac_p <- ggtree(bac_tree, layout = "circular", aes(color=phylum), size = 1.5) %<+%
    select(bac_meta, c(label, phylum, host_species)) +
  scale_colour_manual(values = phylum_palette, name = "MAG phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=2, aes(colour=host_species)) +
  scale_colour_manual(values = species_palette, name = "Host species", na.value = "black") +
  new_scale_color() +
  geom_tippoint(size = 2, aes(color=host_species)) +
  scale_colour_manual(values = species_palette, name = "Host species", na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-6, -6, -6, 1), "cm"), # Remove margins
        legend.position=c(0.06, 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
  color = guide_legend(override.aes = list(size = 5)))

bac_py <- filter(select(bac_meta, c(label, bin, median_damage_model_p)), !is.na(bin))

# Add info
bac_p <- bac_p +
  # Add boxplot with damage patterns
  geom_fruit(data=bac_py, geom=geom_bar, mapping = aes(y=label, x = -log10(median_damage_model_p + 0.001)),
             stat = "identity", axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.3, pwidth = 0.2, alpha = 0.3) +
  new_scale_color() +
  # Add point with collection year
  geom_fruit(data=ages_damage, geom=geom_point, mapping = aes(y=label, colour = Period, shape = Period), size = 2,
             offset = -0.2) +
  #scale_colour_viridis_c(option = "magma", name = "Collection year", na.value = "transparent") +
  #scale_shape(name = "", labels = c("Yes" = "Approximated", "No" = "From records")) +
  new_scale_colour() +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

ggsave(bac_p, file=file.path(subdir, "bac_genome_tree.png"), width = 12, height = 10)

#### Archaea tree ####
# Colour by order and habitat
ar_p <- ggtree(ar_tree) %<+%
    select(ar_meta, c(label, host_species)) +
  new_scale_color() +
  geom_tiplab(size=2, aes(colour=host_species)) +
  scale_colour_manual(values = species_palette, name = "Host species", na.value = "black") +
  new_scale_color() +
  geom_tippoint(size = 2, aes(color=host_species)) +
  scale_colour_manual(values = species_palette, name = "Host species", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5))) 

ar_py <- filter(select(ar_meta, c(label, bin, median_damage_model_p)), !is.na(bin))

# Add info
ar_p <- ar_p +
  # Add boxplot with damage patterns
  geom_fruit(data=ar_py, geom=geom_bar, mapping = aes(y=label, x = -log10(median_damage_model_p + 0.001)),
             stat = "identity", axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.39, pwidth = 0.2, alpha = 0.3) +
  new_scale_color() +
  # Add point with collection year
  geom_fruit(data=ages_damage, geom=geom_point, mapping = aes(y=label, colour = Period, shape = Period), size = 2,
             offset = -0.19) +
  #scale_colour_viridis_c(option = "magma", name = "Collection year", na.value = "transparent") +
  #scale_shape(name = "", labels = c("Yes" = "Approximated", "No" = "From records")) +
  new_scale_colour() +
  theme(legend.position = "none")

ggsave(ar_p, file=file.path(subdir, "ar_genome_tree.png"), width = 6, height = 4)

###########################
#### PLOT AGE V DAMAGE ####
###########################

p <- ggplot(aes(x = Period, y = median_damage_model_p, colour = host_species), data = ages_damage) +
        geom_jitter(alpha = 0.5, width = 0.1) +
        scale_colour_manual(values = species_palette, name = "Host species")

ggsave(p, file=file.path(subdir, "age_v_damage.png"), width = 8, height = 6)
