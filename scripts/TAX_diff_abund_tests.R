##### DIFFERENTIAL ABUNDANCE TESTS #####

#### Run differential abundance analysis on filtered data

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
library(colorspace)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "differential_abundance"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

# Extract OTU table and sample data
otu_table <- as.data.frame(phy_sp_f@otu_table)
sample_data <- as.data.frame(phy_sp_f@sam_data)

# Transpose the OTU table
otu_table_t <- t(otu_table)

#########################
#### PREP TEST INPUT ####
#########################

phy_genus <- phy_sp_f %>% tax_glom("genus")

phy_genus@tax_table[, "genus"] <- make.names(phy_genus@tax_table[, "genus"], unique = TRUE)
taxa_names(phy_genus) <- phy_genus@tax_table[, "genus"]

phy_genus_clr <- phy_genus %>% transform("clr")

#################
#### ANCOMBC ####
#################

# Run for all mammal samples and then for each genus separately

#### ALL WILD VS DOMESTIC ####

# Remove feral sheep and Humans because they are an ambiguous category
# If diff abund genera are found, I will plot the abundances to see
# if they match the domesticated or wild patterns better
phy_ancom <- phy_genus %>% subset_samples(Domestication != "feral" & Domestication != "human")
phy_ancom <- phy_ancom %>% prune_taxa(taxa_sums(phy_ancom) > 0, .)

# Compare wild and domestic animals to humans
phy_ancom@sam_data$Domestication <- factor(phy_ancom@sam_data$Domestication, levels = c("wild", "domestic"))

ancom_all <- ancombc2(data = phy_ancom,
               fix_formula = "Domestication + Genus",
               tax_level = "genus", 
               p_adj_method = "holm", prv_cut = 0.1, 
               group="Domestication",
               struc_zero = FALSE,
               lib_cut = 0,
               verbose = TRUE)

write.csv(ancom_all$res, file = file.path(subdir, "ancombc_all.csv"), quote = FALSE, row.names = FALSE)

#### Within Equus ####
phy_equus <- phy_genus %>% subset_samples(Genus == "Equus")
phy_equus <- phy_equus %>% prune_taxa(taxa_sums(phy_equus) > 0, .)

ancom_equus <- ancombc2(data = phy_equus,
               fix_formula = "Domestication", 
               tax_level = "genus", 
               p_adj_method = "holm", prv_cut = 0.1,
               group="Domestication",
               struc_zero = FALSE,
               lib_cut = 0,
               verbose = TRUE)

write.csv(ancom_equus$res, file = file.path(subdir, "ancom_equus.csv"), quote = FALSE, row.names = FALSE)

#### Within Ovis ####
phy_ovis <- phy_genus %>% subset_samples(Genus == "Ovis" & Domestication != "feral")
phy_ovis <- phy_ovis %>% prune_taxa(taxa_sums(phy_ovis) > 0, .)

ancom_ovis <- ancombc2(data = phy_ovis,
               fix_formula = "Domestication",
               tax_level = "genus",
               p_adj_method = "holm", prv_cut = 0.1,
               struc_zero = FALSE,
               lib_cut = 0,
               verbose = TRUE)

write.csv(ancom_ovis$res, file = file.path(subdir, "ancom_ovis.csv"), quote = FALSE, row.names = FALSE)

#### Within Sus ####

phy_sus <- phy_genus %>% subset_samples(Genus == "Sus")
phy_sus <- phy_sus %>% prune_taxa(taxa_sums(phy_sus) > 0, .)

ancom_sus <- ancombc2(data = phy_sus,
               fix_formula = "Domestication",
               tax_level = "genus",
               p_adj_method = "holm", prv_cut = 0,
               group="Domestication",
               struc_zero = FALSE,
               lib_cut = 0,
               verbose = TRUE)

write.csv(ancom_sus$res, file = file.path(subdir, "ancom_sus.csv"), quote = FALSE, row.names = FALSE)

#### Human vs Wild ####

phy_human <- phy_genus %>% subset_samples(Domestication %in% c("human", "wild"))
phy_human <- phy_human %>% prune_taxa(taxa_sums(phy_human) > 0, .)

phy_human@sam_data$Domestication <- ifelse(phy_human@sam_data$Domestication == "human", "human", "animal")
phy_human@sam_data$Domestication <- factor(phy_human@sam_data$Domestication, levels = c("animal", "human"))

ancom_human <- ancombc2(data = phy_human,
               fix_formula = "Domestication",
               tax_level = "genus",
               p_adj_method = "holm", prv_cut = 0.1,
               group="Domestication",
               struc_zero = FALSE,
               lib_cut = 0,
               verbose = TRUE)

write.csv(ancom_human$res, file = file.path(subdir, "ancom_human.csv"), quote = FALSE, row.names = FALSE)

#########################
#### COMBINE RESULTS ####
#########################

ancom_all$res$dataset <- "Domestic vs Wild (All)"
ancom_equus$res$dataset <- "Horse vs Zebra"
ancom_ovis$res$dataset <- "Sheep vs Argali"
ancom_sus$res$dataset <- "Pig vs Boar"
ancom_human$res$dataset <- "Human vs Wild"

# Combine results
res <- lapply(list(ancom_all$res, ancom_equus$res, ancom_ovis$res, ancom_sus$res, ancom_human$res), 
              function(x) {
                select(x, taxon, contains("Domestication"), dataset) %>%
                rename_with(., ~str_remove_all(., "_Domesticationdomestic") %>% str_remove_all(., "_Domesticationhuman"))
                }) %>% bind_rows()

write.csv(res, file.path(subdir, "ancom_res.csv"), quote = FALSE, row.names = FALSE)

# Keep only taxa that are differentially abundant in any of the analyses
signif_taxa <- res %>% filter(q < 0.05 & passed_ss) %>% pull(taxon) %>% unique

# Order by lfc in the largest subset
taxa_order <- res %>% filter(taxon %in% signif_taxa) %>%
               group_by(taxon) %>% select(taxon, dataset, lfc) %>%
               pivot_wider(names_from = dataset, values_from = lfc) %>%
               arrange(`Domestic vs Wild (All)`, `Human vs Wild`) %>% pull(taxon)

res_filt <- res %>% filter(taxon %in% signif_taxa) %>%
            mutate(taxon = factor(taxon, levels = taxa_order)) %>%
            mutate(signif = case_when(q < 0.001 & passed_ss ~ "***",
                                      q < 0.01 & passed_ss ~ "**",
                                      q < 0.05 & passed_ss ~ "*",
                                      q < 0.1 & passed_ss ~ ".",
                                      TRUE ~ "")) %>%
            # Make x labels look a bit nicer
            mutate(dataset = gsub(" vs ", "\nvs ", dataset)) %>%
            mutate(dataset = factor(dataset, levels = c("Domestic\nvs Wild (All)", "Horse\nvs Zebra", "Sheep\nvs Argali", "Pig\nvs Boar",  "Human\nvs Wild")))

# Plot heatmap
p <- ggplot(data = res_filt, aes(x = dataset, y = taxon, fill = lfc)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Log-fold change", na.value = "transparent") +
    geom_text(aes(label = signif), color = "black", size = 3) +
    theme(legend.position = "top", legend.text = element_text(angle = 45, vjust = 0.5),
          panel.background = element_rect(fill = "grey90"), panel.grid = element_blank(),
          axis.text.y = element_text(size = 8))

ggsave(p, filename = file.path(subdir, "ancom_res_heatmap.png"), width = 6, height = 10)

# Also plot lfc's as scatterplots to show if they correlate
data_palette <- darken(species_palette[c("Equus quagga", "Ovis ammon", "Sus scrofa", "Homo sapiens")])
names(data_palette) <- c("Horse\nvs Zebra", "Sheep\nvs Argali", "Pig\nvs Boar",  "Human\nvs Wild")

p <- filter(res_filt, c(!dataset %in% c("Domestic\nvs Wild (All)", "Human\nvs Wild"))) %>%
    ggplot(aes(y = taxon, x = lfc)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = dataset, shape = dataset), size = 3, alpha = 0.8) +
    geom_line(aes(group = taxon), linetype = "dotted", size = 0.3) +
    scale_colour_manual(values = data_palette, name = "Dataset") +
    scale_shape_manual(values = c(0, 2, 4), name = "Dataset") +
    theme(legend.position = "top", axis.text.y = element_text(size = 8)) +
    xlab("Log-fold change in domestic animals")

ggsave(p, filename = file.path(subdir, "ancom_lfc_comparison.png"), width = 6, height = 10)

#########################
#### PLOT ABUNDANCES ####
#########################

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_genus_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Abundance") %>%
              left_join(rownames_to_column(select(data.frame(phy_genus_clr@sam_data), Common.name, Species, Genus, Domestication), "Sample"), by = "Sample") %>%
              mutate(Common.name = case_when(Domestication == "feral" ~ "Feral sheep",
                                          TRUE ~ Common.name))

ancom_abund <- res_filt %>% filter(q < 0.05 & passed_ss) %>%
            mutate(association = case_when(lfc < 0 ~ "wild+",
                                           lfc > 0 ~ "wild-")) %>%
            mutate(dataset = paste("comparison with", str_remove_all(dataset, "\n.*"))) %>%
            mutate(label = paste(association, dataset)) %>%
            # Combine labels for different results from the same OTU
            group_by(taxon) %>%
            summarise(label = paste(label, collapse="\n")) %>%
            # Add abundances
            left_join(select(abundances, OTU, Abundance, Species, Common.name, Genus, Domestication), by = c("taxon" ="OTU"))

# Reorder host species
species_levels <- abundances %>% data.frame %>% arrange(Genus, Domestication) %>% select(Species, Common.name) %>% unique
ancom_abund$Common.name <- factor(ancom_abund$Common.name, levels = species_levels$Common.name)

p <- ggplot(ancom_abund, aes(x = Common.name, y = Abundance, fill = Species, colour = Species)) +
    geom_boxplot(alpha = 0.8, size = 0.5, outliers = FALSE) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.8) +    scale_fill_manual(values = species_palette, name = "Species") +
    scale_colour_manual(values = darken(species_palette), name = "Species") +
    facet_wrap(~ paste(as.character(taxon), label, sep = "\n"), ncol = 4, scales = "free_y") +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "ancom_abundances.png"), width = 10, height = 30)
