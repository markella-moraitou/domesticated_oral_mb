##### DIFFERENCES IN ABUNDANCE OF STREPTOCOCCUS CLADES #####

#### Check for differences in the abundance of the clades Mitis, Sanguinis, Salivarius, Anginosus, Mutans and Pyogenes

################
#### SET UP ####
################

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "strep_clades"))
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

phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))

# Define the clades based on Richards et al. (2014) Phylogenomics and the Dynamic Genome Evolution of the Genus Streptococcus
strep_df <- list("Mitis" = c("Streptococcus mitis", "Streptococcus oralis", "Streptococcus infantis", "Streptococcus parasanguinis"),
                   "Sanguinis" = c("Streptococcus gordonii", "Streptococcus sanguinis", "Streptococcus cristatus"),
                   "Salivarius" = c("Streptococcus salivarius", "Streptococcus thermophilus"),
                   "Anguinosus" = c("Streptococcus intermedius", "Streptococcus constellatus"),
                   "Mutans" = c("Streptococcus ratti", "Streptococcus mutans", "Streptococcus macacae")) %>% unlist %>% as.data.frame %>%
                   rownames_to_column("clade") %>% mutate(clade = str_remove(clade, "[0-9]")) %>%
                   rename("species" = ".")

#################
#### PROCESS ####
#################

strep_abundances <- psmelt(subset_taxa(phy_sp_f, genus == "Streptococcus")) %>%
    # Remove suffixes to match 
    mutate(species = str_remove(species, "_.*")) %>%
    # Add clade info
    left_join(strep_df) %>%
    mutate(clade = case_when(is.na(clade) ~ "Unknown", TRUE ~ clade)) %>%
    group_by(Sample) %>% mutate(rel_abundance = Abundance/sum(Abundance))


write.csv(strep_abundances, file = file.path(subdir, "strep_clade_abundance.tsv"), quote = TRUE, row.names = FALSE)

# By amylase type
p <- ggplot(strep_abundances, aes(y = Sample, x = rel_abundance, fill =  forcats::fct_rev(clade))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(y = "Sample",
       x = "Relative abundance of streptococcus clades",
       fill = "Streptococcus clade") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))

ggsave(filename = file.path(subdir, "strep_clade_abundance.png"), width = 5, height = 8)
