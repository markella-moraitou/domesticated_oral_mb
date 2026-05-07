##### LOSS & GAIN OF TAXA #####

#### Identify consistent loss and gain of taxa during domestication

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(microbiome)
library(tidyr)
library(tibble)
library(stringr)
library(ANCOMBC)
library(ggplot2)
library(ggVennDiagram)
library(scales)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "loss_gain"))
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

#########################
#### PREP TEST INPUT ####
#########################

phy_genus <- phy_sp_f %>% tax_glom("genus")

phy_genus@tax_table[, "genus"] <- make.names(phy_genus@tax_table[, "genus"], unique = TRUE)
taxa_names(phy_genus) <- phy_genus@tax_table[, "genus"]

phy_genus <- phy_genus %>% 
      subset_samples(Genus != "Homo" & Domestication != "feral") # Remove humans and feral samples

#################
#### ANCOMBC ####
#################

# Run ANCOMBC with group = Species to identify structural zeroes

ancom <- ancombc2(data = phy_genus,
               fix_formula = "Species",
               tax_level = "genus", 
               p_adj_method = "holm", prv_cut = 0, 
               group="Species",
               struc_zero = TRUE,
               lib_cut = 0,
               verbose = TRUE)
     
str_zero <- ancom$zero_ind %>% rename_with(~ gsub(")$", "", str_remove(.x, "structural_zero.*Species = ")), .cols = everything())

write.csv(str_zero, file = file.path(subdir, "ancombc_structural_zeroes.csv"), quote = FALSE, row.names = FALSE)

##############################
#### IDENTIFY LOSS & GAIN ####
##############################

metadata <- as_tibble(sample_data(phy_genus)) %>%
  select(Species, Genus, Domestication) %>% unique

# Keep taxa that are not always present
str_zero_filt <- str_zero %>% pivot_longer(cols = -taxon, names_to = "Species", values_to = "structural_zero") %>%
  # Add metadata
  left_join(metadata, by = "Species", relationship = "many-to-one")

loss_gain <- 
  str_zero_filt %>%
  group_by(taxon, Genus) %>%
  mutate(present = !structural_zero, structural_zero = NULL) %>%
  # Label as loss and gain
  summarise(loss_gain = case_when(
    present[Domestication == "wild"] & !present[Domestication == "domestic"] ~ "loss",
    present[Domestication == "domestic"] & !present[Domestication == "wild"] ~ "gain"
  )) %>% group_by(taxon) %>%
  mutate(
    losses = sum(loss_gain == "loss", na.rm = TRUE),
    gains = sum(loss_gain == "gain", na.rm = TRUE)
  ) %>%
  filter(!is.na(loss_gain))

write.csv(loss_gain, file = file.path(subdir, "loss_gain_taxa.csv"), quote = FALSE, row.names = FALSE)

# Summarise the behaviour of each taxon
lg_summary <-
  loss_gain %>% group_by(taxon) %>%
  summarise(
    losses = sum(loss_gain == "loss", na.rm = TRUE),
    gains = sum(loss_gain == "gain", na.rm = TRUE)
  )

write.csv(lg_summary, file = file.path(subdir, "loss_gain_taxa_summary.csv"), quote = FALSE, row.names = FALSE)

p <- lg_summary %>% group_by(losses, gains) %>% summarise(n_taxa = n_distinct(taxon)) %>%
        ggplot(aes(x = losses, y = gains, fill = n_taxa)) +
        geom_tile() +
        scale_fill_viridis_c(option = "mako", direction = -1, name = "Number of taxa") +
        geom_text(aes(label = n_taxa), color = "grey40", size = 8) +
        annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -0.5, ymax = 0.5, linewidth = 2, colour = "black", fill = "transparent") +
        annotate("rect", xmin = -0.5, xmax = 0.5, ymin = 0.5, ymax = 3.5, linewidth = 2, colour = "black", fill = "transparent") +
        xlab("Taxon losses") + ylab("Taxon gains") + theme(legend.position = "none")

ggsave(filename = file.path(subdir, "loss_gain_taxa_summary_heatmap.png"), plot = p, width = 4, height = 4)

############################
#### Plot Venn Diagrams ####
############################

#### Losses ####

# Keep taxa that are consistently lost during domestication
loss_taxa <- lg_summary %>% filter(losses >= 1 & gains == 0) %>% pull(taxon)

# Collect loss taxa in a list
loss_taxa_list <- list()

for (g in c("Equus", "Ovis", "Sus")) {
    loss_taxa_list[[g]] <-
        loss_gain %>% filter(Genus == g & loss_gain == "loss" & taxon %in% loss_taxa) %>%
        pull(taxon) %>% unique
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
gain_taxa <- lg_summary %>% filter(gains >= 1 & losses == 0) %>% pull(taxon)

# Collect gain taxa in a list
gain_taxa_list <- list()

for (g in c("Equus", "Ovis", "Sus")) {
    gain_taxa_list[[g]] <-
        loss_gain %>% filter(Genus == g & loss_gain == "gain" & taxon %in% gain_taxa) %>%
        pull(taxon) %>% unique
}

gain_taxa_venn <- 
  ggVennDiagram(gain_taxa_list, set_color = genus_palette[names(loss_taxa_list)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. gains") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")

ggsave(gain_taxa_venn, file=file.path(subdir, "gain_taxa_venn.png"),
       device = "png", width = 4, height = 4)

##################################
#### PLOT RELATIVE ABUNDANCES ####
##################################

lg_summary_filtered <- 
  lg_summary %>%
  filter((losses == 3 | gains == 3)) %>%
  mutate(
    loss_gain = case_when(
      losses == 3 ~ "loss",
      gains == 3 ~ "gain"
    )
  )

lg_abundances <- transform(phy_genus, "compositional")@otu_table[lg_summary_filtered$taxon, ] %>% data.frame %>% rownames_to_column("taxon") %>%
      pivot_longer(cols = -taxon, names_to = "Sample", values_to = "Relative abundance") %>%
      right_join(
        as_tibble(sample_data(phy_genus)) %>% select(new_name, Species),
        by = c("Sample" = "new_name")
      )

p <- ggplot(aes(x = Species, colour = Species, y = `Relative abundance`), data = lg_abundances) +
      geom_jitter(alpha = 0.7, width = 0.1, height = 0) +
      scale_colour_manual(values = species_palette) +
      scale_y_continuous(labels = label_percent()) +
      facet_wrap(~ paste0(taxon, " (", lg_summary_filtered$loss_gain[match(lg_abundances$taxon, lg_summary_filtered$taxon)], ")"), scales = "free_y") +
      theme(legend.position = "none", axis.text.x = element_text(hjust = 1))

ggsave(p, file = file.path(subdir, "loss_gain_taxa_abundances.png"), width = 8, height = 8)      
