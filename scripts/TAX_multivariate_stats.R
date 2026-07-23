##### EXPLORE FILTERED & DECONTAMINATED DATA #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(reshape2)
library(microViz)
library(rphylopic)
library(ggplot2)
library(ggExtra)
library(ggnewscale)
library(vegan)
library(microbiomeutilities)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output", "community_analysis")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "multivariate_stats")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

# Get functions
source(file.path("modules", "ordination_functions.R"))
source(file.path("modules", "phylo_functions.R"))

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

#########################
#####  BETA DISPER  #####
#########################

#### For host group ####

tukey_results <- data.frame()

# Compare species
disp <- betadisper(vegdist(t(otu_table(phy_sp_f_clr)), method = "euclidean"), group = phy_sp_f_clr@sam_data$Group)

disp_anova <- anova(disp)

disp_tukey <- TukeyHSD(disp, which = "group", ordered = FALSE)$group %>% data.frame %>% rownames_to_column("Comparison") %>%
  separate(Comparison, into = c("Group1", "Group2"), sep = "-")

tukey_results <- rbind(tukey_results, disp_tukey)

write.csv(tukey_results, file = file.path(subdir, "betadisper_tukey_results.csv"), row.names = FALSE, quote = FALSE)

disp_df <- data.frame(Sample = sample_names(phy_sp_f_clr),
                      Group = phy_sp_f_clr@sam_data$Group,
                      Species = phy_sp_f_clr@sam_data$Species,
                      Genus = phy_sp_f_clr@sam_data$Genus,
                      Distance = disp$distances)

p <- ggplot(data = disp_df, aes(x = Group, y = Distance, fill = Species)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Species)) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    scale_fill_manual(values = species_palette, name = "Species") +
    theme(legend.position = "none", axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "betadisper.png"), p, width = 5, height = 5)

#### For domesticated versus wild (exclude feral) ####

phy_dw <- subset_samples(phy_sp_f_clr, Domestication %in% c("wild", "domestic"))

disp <- betadisper(vegdist(t(otu_table(phy_dw)), method = "euclidean"), group = phy_dw@sam_data$Domestication)

disp_anova <- anova(disp)

disp_df <- data.frame(Sample = sample_names(phy_dw),
                      Domestication = phy_dw@sam_data$Domestication,
                      Distance = disp$distances)

p <- ggplot(data = disp_df, aes(x = Domestication, y = Distance)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.5, width = 0.2) +
    #scale_fill_manual(values = species_palette, name = "Species") +
    theme(legend.position = "none", axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "betadisper_dom.png"), p, width = 5, height = 5)

#########################
#####  COMPOSITION  #####
#########################

#### PHYLUM COMPOSITION ####
phy_phylum <- tax_glom(phy_sp_f, taxrank = "phylum")
taxa_names(phy_phylum) <- as.vector(phy_phylum@tax_table[,"phylum"])

# Get only 10 most common phyla and turn rest to other
phylum_grouped = data.frame(abundance = taxa_sums(phy_phylum), superkingdom = phy_phylum@tax_table[,"superkingdom"]) %>% rownames_to_column("phylum") %>% arrange(-abundance) %>%
  mutate(phylum_grouped = ifelse(row_number() > 5, paste("Other", superkingdom, sep = " "), phylum))

# Change phylum names to grouped names and reaggregate
phy_phylum@tax_table[,"phylum"] <- phylum_grouped$phylum_grouped[match(phy_phylum@tax_table[,"phylum"], phylum_grouped$phylum)]

# Aggregate again
phy_phylum <- tax_glom(phy_phylum, taxrank = "phylum")
taxa_names(phy_phylum) <- phy_phylum@tax_table[,"phylum"] 

# Melt and turn phyla into a factor and reorder
phy_phylum_melt <- psmelt(transform(phy_phylum, "compositional"))
phy_phylum_melt$OTU <- factor(phy_phylum_melt$OTU , levels=names(phylum_palette))

# Order plots facets
group_order <- phy_phylum_melt %>% select(Group, Species, Genus, Domestication) %>% unique %>%
  arrange(Genus, desc(Domestication)) %>% pull(Group) %>% unique

phy_phylum_melt$Group <- factor(phy_phylum_melt$Group, levels=group_order)

# Order by Pseudomonadota
sample_levels <- select(phy_phylum_melt, c(Sample, Species, Genus, Domestication, OTU, Abundance)) %>% filter(OTU == "Pseudomonadota") %>%
  arrange(Genus, desc(Domestication), desc(Abundance)) %>% pull(Sample)

phy_phylum_melt$Sample <- factor(phy_phylum_melt$Sample, levels=sample_levels)

p = ggplot(data = phy_phylum_melt, aes(x = Abundance, y = Sample, fill = OTU)) +
  geom_bar(stat = "identity") +
  facet_grid(Group~., space = "free_y", scales = "free_y", switch = "y") +
  scale_fill_manual(values=phylum_palette, name = "Phylum") +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "left", legend.title.position = "top", legend.key.spacing.x = unit(0.5, "cm"),
        legend.direction = "vertical", 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        strip.background = element_blank(), strip.text = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE)) +
  xlab("")

## Get sample_metadata
species_bar <- phy_phylum_melt %>% select(Sample, Species, Group, Common.name, Domestication) %>% unique %>%
               arrange(Sample) %>%
               # Choose one label per species (in the middle of the species group)
               group_by(Group, Species) %>% mutate(order = row_number()) %>%
               mutate(label = ifelse(order == floor(mean(order)), as.character(Group), "")) %>%
               mutate(label = ifelse(is.na(label), as.character(Group), label))

p_bar <-
  ggplot(data = species_bar, aes(y = Sample, x=1, fill = Species, group = Group)) +
  geom_tile() + scale_fill_manual(values = species_palette, name = "") +
  facet_grid(rows = vars(Group), scales = "free", space = "free", switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", 
        strip.background = element_blank(),
        strip.text = element_blank()) + scale_y_discrete(position = "right", label = setNames(species_bar$label, species_bar$Sample)) +
    xlab("") + ylab("")

ggsave(filename = file.path(subdir, "phy_sp_f_composition.png"), device="png", width=6, height=6,
       plot_grid(p, p_bar, ncol = 2, align = "h", rel_widths = c(3.5, 1.5)))

#############
#### PCA ####
#############

#### Get shape scales for plotting ####

dom_shape_palette <- c("domestic" = 1, "wild" = 16, "feral" = 6, "human" = 8)

#### CLR ABUNDANCES ####

ord <- ord_calc(phy_sp_f_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_clr_screeplot.png"), p, width=3, height=3)

# Color by species
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="Species", shape="Domestication", arrows_scaling = 1, type = "PCA") +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication")

ggsave(file.path(subdir, "PCA_clr_1_2.png"), p, width=10, height=6)

#### PHILR ABUNDANCES ####

ord <- ord_calc(phy_sp_f_pa, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_pa_screeplot.png"), p, width=3, height=3)

# Color by species
p <- custom_ord_plot(phy_sp_f_pa, ord, colour="Species", shape="Domestication", arrows_scaling = 1, type = "PCA") +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication")

ggsave(file.path(subdir, "PCA_pa_1_2.png"), p, width=10, height=6)

#### All data ####
ord <- ord_calc(phy_sp_philr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_philr.png"), p, width=3, height=3)

# Color by species
p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0, auto_caption = NA) +
  geom_point(size = 2, alpha = 0.8, aes(colour = Species, shape = Domestication)) +
  custom_theme() +
  scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
  scale_color_manual(values=species_palette, name = "Species") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_philr), aes(colour = Species), uuid = centroids(ord@ord, phy_sp_philr)$uid, width = 0.5, alpha = 0.8) +
  guides(shape = guide_legend(ncol = 1),
        colour = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(shape = c(16, 1, 16, 1, 16, 1, 8))))

ggsave(file.path(subdir, "PCA_philr_1_2.png"), p, width=10, height=6)

#############
#### RDA ####
#############

#### CLR ABUNDANCES #####

phy_sp_f_clr <- phy_sp_f_clr %>%
        ps_mutate(sheep = (Genus == "Ovis" & Domestication %in% c("domestic", "feral")),
                  argali = (Genus == "Ovis" & Domestication == "wild"),
                  horse = (Genus == "Equus" & Domestication == "domestic"),
                  zebra = (Genus == "Equus" & Domestication == "wild"),
                  pig = (Genus == "Sus" & Domestication == "domestic"),
                  boar = (Genus == "Sus" & Domestication == "wild"))

# Species traits to use as constraints
species_traits <- c("sheep", "argali",
                    "horse", "zebra",
                    "pig", "boar")

# Ordinate using all data
ord <- ord_calc(phy_sp_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_clr_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

# Color by species
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="Species", shape="Domestication", type = "RDA")

ggsave(file.path(subdir, "RDA_clr_1_2.png"), p, width=10, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_sp_f_clr)[["plot"]]
ggsave(file.path(subdir, "RDA_clr_taxa_1_2.png"), p, width=10, height=6)

write.csv(taxa_plot(ord, phy_sp_f_clr)[["data"]], file = file.path(subdir, "RDA_clr_taxa_scores.csv"), row.names = FALSE, quote = TRUE)

## Without humans

phy_sp_f_clr_nohuman <- phy_sp_f_clr %>% subset_samples(Species != "Homo sapiens")

# Ordinate using all data
ord <- ord_calc(phy_sp_f_clr_nohuman, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_clr_nohuman_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

# Color by species
p <- custom_ord_plot(phy_sp_f_clr_nohuman, ord, colour="Species", shape="Domestication", type = "RDA")

ggsave(file.path(subdir, "RDA_clr_nohuman_1_2.png"), p, width=10, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_sp_f_clr)[["plot"]]
ggsave(file.path(subdir, "RDA_clr_nohuman_taxa_1_2.png"), p, width=10, height=6)

write.csv(taxa_plot(ord, phy_sp_f_clr)[["data"]], file = file.path(subdir, "RDA_clr_nohuman_taxa_scores.csv"), row.names = FALSE, quote = TRUE)

#### PRESENCE-ABSENCE ####

phy_sp_f_pa <- phy_sp_f_pa %>%
        ps_mutate(sheep = (Genus == "Ovis" & Domestication %in% c("domestic", "feral")),
                  argali = (Genus == "Ovis" & Domestication == "wild"),
                  horse = (Genus == "Equus" & Domestication == "domestic"),
                  zebra = (Genus == "Equus" & Domestication == "wild"),
                  pig = (Genus == "Sus" & Domestication == "domestic"),
                  boar = (Genus == "Sus" & Domestication == "wild"))

# Ordinate using all data
ord <- ord_calc(phy_sp_f_pa, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_pa_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

# Color by species
p <- custom_ord_plot(phy_sp_f_pa, ord, colour="Species", shape="Domestication", type = "RDA")

ggsave(file.path(subdir, "RDA_pa_1_2.png"), p, width=10, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_sp_f_clr)[["plot"]]
ggsave(file.path(subdir, "RDA_pa_taxa_1_2.png"), p, width=10, height=6)

write.csv(taxa_plot(ord, phy_sp_f_clr)[["data"]], file = file.path(subdir, "RDA_pa_taxa_scores.csv"), row.names = FALSE, quote = TRUE)

## Without humans

phy_sp_f_pa_nohuman <- phy_sp_f_pa %>% subset_samples(Species != "Homo sapiens")

# Ordinate using all data
ord <- ord_calc(phy_sp_f_pa_nohuman, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_pa_nohuman_screeplot.png"), p, width=3, height=3)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_pa_nohuman_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

# Color by species
p <- custom_ord_plot(phy_sp_f_pa_nohuman, ord, colour="Species", shape="Domestication", type = "RDA")

ggsave(file.path(subdir, "RDA_pa_nohuman_1_2.png"), p, width=10, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_sp_f_pa)[["plot"]]
ggsave(file.path(subdir, "RDA_pa_nohuman_taxa_1_2.png"), p, width=10, height=6)

write.csv(taxa_plot(ord, phy_sp_f_pa)[["data"]], file = file.path(subdir, "RDA_pa_nohuman_taxa_scores.csv"), row.names = FALSE, quote = TRUE)

#### PHILR ABUNDANCES #####

phy_sp_philr <- phy_sp_philr %>%
        ps_mutate(sheep = (Genus == "Ovis" & Domestication %in% c("domestic", "feral")),
                  argali = (Genus == "Ovis" & Domestication == "wild"),
                  horse = (Genus == "Equus" & Domestication == "domestic"),
                  zebra = (Genus == "Equus" & Domestication == "wild"),
                  pig = (Genus == "Sus" & Domestication == "domestic"),
                  boar = (Genus == "Sus" & Domestication == "wild"))

# Ordinate using all data
ord <- ord_calc(phy_sp_philr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_philr_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

centroids_ord <- centroids(ord@ord, phy_sp_philr)

pp_height <- diff(range(vegan::scores(ord@ord, display="sites", choices=2)[,1]))/15

# Color by species
p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0, auto_caption = NA,
              constraint_vec_style = vec_constraint(alpha = 0, linewidth = 0), constraint_lab_style = list(alpha = 0, size = 0, linewidth = 0)) +
    geom_point(size = 2, alpha = 0.8, aes(colour = Species, shape = Domestication)) +
    custom_theme() +
    geom_phylopic(data = centroids_ord, aes(colour = Species), uuid = centroids_ord$uid, fill = "transparent", height = pp_height) +
    scale_colour_manual(values=species_palette, name = "Species") +
    scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
    guides(shape = guide_legend(ncol = 1), colour = guide_legend(ncol = 3, byrow = FALSE))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_philr_1_2.png"), p, width=10, height=6)

## Without humans

phy_sp_philr_nohuman <- phy_sp_philr %>% subset_samples(Species != "Homo sapiens")

# Ordinate using all data
ord <- ord_calc(phy_sp_philr_nohuman, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_philr_nohuman_screeplot.png"), p, width=3, height=3)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "RDA_philr_nohuman_screeplot.png"), p, width=3, height=3)

## SAMPLE PLOTS

centroids_ord <- centroids(ord@ord, phy_sp_philr_nohuman)

# Color by species
p <- ord_plot(ord, colour="Species", shape="Domestication", alpha = 0, auto_caption = NA,
              constraint_vec_style = vec_constraint(alpha = 0, linewidth = 0), constraint_lab_style = list(alpha = 0, size = 0, linewidth = 0)) +
    geom_point(size = 2, alpha = 0.8, aes(colour = Species, shape = Domestication)) +
    custom_theme() +
    geom_phylopic(data = centroids_ord, aes(colour = Species), uuid = centroids_ord$uid, fill = "transparent", height = pp_height) +
    scale_colour_manual(values=species_palette, name = "Species") +
    scale_shape_manual(values=dom_shape_palette, name = "Domestication") +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
    guides(shape = guide_legend(ncol = 1), colour = guide_legend(ncol = 3, byrow = FALSE))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_philr_nohuman_1_2.png"), p, width=10, height=6)

###################
#### PERMANOVA ####
###################

# Remove humans for permanova analyses

#### CLR ABUNDANCES ####

# Explanatory variables
otus <- t(as.data.frame(subset_samples(phy_sp_f_clr, Species != "Homo sapiens")@otu_table))

sample_data <- as.data.frame(subset_samples(phy_sp_f_clr, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ reads + genus * dom,
        permutations = 1000, by = "term", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr.csv"), row.names = TRUE, quote = TRUE)

#### PhILR ####
# Explanatory variables

otus <- as.data.frame(subset_samples(phy_sp_philr, Species != "Homo sapiens")@otu_table)

sample_data <- as.data.frame(subset_samples(phy_sp_philr, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ reads + genus * dom,
        permutations = 1000, by = "term", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_philr.csv"), row.names = TRUE, quote = TRUE)

#### PRESENCE-ABSENCE ####
# Explanatory variables

otus <- t(as.data.frame(subset_samples(phy_sp_f, Species != "Homo sapiens")@otu_table))

sample_data <- as.data.frame(subset_samples(phy_sp_f, Species != "Homo sapiens")@sam_data)

genus <- sample_data$Genus
dom <- sample_data$Domestication
reads <- sample_data$contig_reads_count

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(otus ~ reads + genus * dom,
        permutations = 1000, by = "term", method = "jaccard")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_pa.csv"), row.names = TRUE, quote = TRUE)

##########################
#### PHYLOPICS LEGEND ####
##########################

phylopics$Common.name <- phy_sp_f@sam_data$Common.name[match(phylopics$Species, as.vector(phy_sp_f@sam_data$Species))]

p <- ggplot(phylopics[!is.na(phylopics$Common.name),], aes(y = Common.name, x = 1, colour = Species)) +
    geom_phylopic(aes(uuid = uid), height = 0.5) +
    scale_color_manual(values = species_palette, name = "Species") +
    theme_void() + theme(axis.title = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_text(size = 5, hjust = 1),
                         legend.position = "none",
                         plot.background = element_rect(fill = "transparent", color = "transparent"),
                         panel.background = element_rect(fill = "transparent", color = "transparent"))

ggsave(file.path(subdir, "phylopics_legend.png"), p, width=1.5, height=1.5)

#################
#### HEATMAP ####
#################

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

set.seed(123)
png(filename = file.path(subdir, "heatmap_all_taxa.png"), width=16, height=20, units="in", res=300)
plot_taxa_heatmap(phy_sp_f, subset.top=ntaxa(phy_sp_f), transformation="clr",
                  VariableA=c("Species", "unmapped_count"),
                  annotation_colors = list("Species" = species_palette),
                  show_rownames = FALSE, cluster_cols = TRUE,
                  show_colnames = FALSE, heatcolors = grad_palette)$plot
dev.off()
