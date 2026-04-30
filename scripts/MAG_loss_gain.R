##### LOSS & GAIN OF TAXA #####

#### Identify consistent loss and gain of taxa during domestication

################
#### SET UP ####
################

library(dplyr)
library(tidyr)
library(microbiome)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggVennDiagram)
library(scales)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "output", "mags"))
subdir <- normalizePath(file.path(outdir, "mag_loss_gain"))

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Metadata
sample_metadata <- read.csv(file.path(indir, "sample_metadata.csv"), header=TRUE)

# Bin metadata
bac_meta <- read.table(file.path(outdir, "bac_meta_drep.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(outdir, "ar_meta_drep.tsv"), sep="\t", header=TRUE)

# Presence absence of MAGs
presabs <- read.csv(file.path(outdir, "mag_mapping_stats", "hq_mag_presence_per_host.csv"))

###################
#### PREP DATA ####
###################

presabs <- mutate(presabs, present = TRUE)

metadata <- sample_metadata %>% select(Species, Domestication) %>%
        filter(Domestication != "feral") %>% distinct

##############################
#### IDENTIFY LOSS & GAIN ####
##############################

# Keep taxa that are not always present
presabs_expanded <- presabs %>% select(host_species, label) %>% 
  expand(host_species, label) %>%
  left_join(select(presabs, c(host_species, label, present)), by = c("host_species", "label")) %>%
  mutate(present = case_when(is.na(present) ~ FALSE, TRUE ~ present)) %>%
  rename(Species = host_species, MAG = label) %>%
  # Add metadata
  left_join(metadata, by = "Species", relationship = "many-to-one") %>%
  # Get host genus
  mutate(Genus = str_extract(Species, "^[A-Za-z]+")) %>%
  # Remove humans
  filter(Genus != "Homo")

loss_gain <- 
  presabs_expanded %>%
  group_by(MAG, Genus) %>%
  # Label as loss and gain
  summarise(loss_gain = case_when(
    present[Domestication == "wild"] & !present[Domestication == "domestic"] ~ "loss",
    present[Domestication == "domestic"] & !present[Domestication == "wild"] ~ "gain"
  )) %>% group_by(MAG) %>%
  mutate(
    losses = sum(loss_gain == "loss", na.rm = TRUE),
    gains = sum(loss_gain == "gain", na.rm = TRUE)
  ) %>%
  filter(!is.na(loss_gain))

write.csv(loss_gain, file = file.path(subdir, "loss_gain_taxa.csv"), quote = FALSE, row.names = FALSE)

# Summarise the behaviour of each taxon
lg_summary <-
  loss_gain %>% group_by(MAG) %>%
  summarise(
    losses = sum(loss_gain == "loss", na.rm = TRUE),
    gains = sum(loss_gain == "gain", na.rm = TRUE)
  )

write.csv(lg_summary, file = file.path(subdir, "loss_gain_taxa_summary.csv"), quote = FALSE, row.names = FALSE)

p <- lg_summary %>% group_by(losses, gains) %>% summarise(n_MAGs = n_distinct(MAG)) %>%
        ggplot(aes(x = losses, y = gains, fill = n_MAGs)) +
        geom_tile() +
        scale_fill_viridis_c(option = "mako", direction = -1, name = "Number of MAGs") +
        geom_text(aes(label = n_MAGs), color = "grey40", size = 8) +
        annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -0.5, ymax = 0.5, linewidth = 2, colour = "black", fill = "transparent") +
        annotate("rect", xmin = -0.5, xmax = 0.5, ymin = 0.5, ymax = 3.5, linewidth = 2, colour = "black", fill = "transparent") +
        xlab("MAG losses") + ylab("MAG gains") + theme(legend.position = "none")

ggsave(filename = file.path(subdir, "loss_gain_taxa_summary_heatmap.png"), plot = p, width = 4, height = 4)

############################
#### Plot Venn Diagrams ####
############################

#### Losses ####

# Keep taxa that are consistently lost during domestication
loss_taxa <- lg_summary %>% filter(losses >= 1 & gains == 0) %>% pull(MAG)

# Collect loss taxa in a list
loss_taxa_list <- list()

for (g in c("Equus", "Ovis", "Sus")) {
    loss_taxa_list[[g]] <-
        loss_gain %>% filter(Genus == g & loss_gain == "loss" & MAG %in% loss_taxa) %>%
        pull(MAG) %>% unique
}

genus_palette <- species_palette[c("Equus caballus", "Ovis aries", "Sus domesticus")]
names(genus_palette) <- c("Equus", "Ovis", "Sus")

loss_taxa_venn <- 
  ggVennDiagram(loss_taxa_list, set_color = genus_palette[names(loss_taxa_list)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. losses") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")

ggsave(loss_taxa_venn, file=file.path(subdir, "loss_taxa_venn.png"),
       device = "png", width = 4, height = 4)

#### Gains ####

# Keep taxa that are consistently gained during domestication
gain_taxa <- lg_summary %>% filter(gains >= 1 & losses == 0) %>% pull(MAG)

# Collect gain taxa in a list
gain_taxa_list <- list()

for (g in c("Equus", "Ovis", "Sus")) {
    gain_taxa_list[[g]] <-
        loss_gain %>% filter(Genus == g & loss_gain == "gain" & MAG %in% gain_taxa) %>%
        pull(MAG) %>% unique
}

gain_taxa_venn <- 
  ggVennDiagram(gain_taxa_list, set_color = genus_palette[names(loss_taxa_list)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. gains") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")

ggsave(gain_taxa_venn, file=file.path(subdir, "gain_taxa_venn.png"),
       device = "png", width = 4, height = 4)