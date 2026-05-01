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
library(MCMCglmm)
library(coda)
library(parallel)
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

# Source MCMCglmm wrapper and output processing functions
source(file.path("modules", "mcmcglmm_functions.R"))

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

# The analysis will be only run on the most abundant taxa (top 100 per host genus)
# to reduce the number of tests and make the MCMCglmm run in a reasonable time frame.
# The analysis will also be run between humans and wild animals, only for any taxa found in the top abundant taxa
# To infer if they reflect 'humanisation' of the microbiome

phy_genus <- phy_sp_f %>% tax_glom("genus")

phy_genus@tax_table[, "genus"] <- make.names(phy_genus@tax_table[, "genus"], unique = TRUE)
taxa_names(phy_genus) <- phy_genus@tax_table[, "genus"]

phy_genus_clr <- phy_genus %>% transform("clr")

# Remove feral sheep and Humans because they are an ambiguous category
# If diff abund genera are found, I will plot the abundances to see
# if they match the domesticated or wild patterns better
phy_da <- phy_genus_clr %>% subset_samples(Domestication != "feral" & Domestication != "human")
phy_da <- phy_da %>% prune_taxa(taxa_sums(phy_da) > 0, .)

phy_da@sam_data$Domestication <- factor(phy_da@sam_data$Domestication, levels = c("wild", "domestic"))

# Get as data frame
data <- psmelt(phy_da) %>%
        select(OTU, Abundance, Sample, Species, Genus, Domestication, unmapped_count, contig_reads_count)

## Create subsets per genus for the within-genus analyses

# Equus
data_equus <- data %>% filter(Genus == "Equus") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

# Get most abundant taxa
top_equus <- data_equus %>% group_by(OTU) %>%
  summarise(av_abundance = mean(Abundance)) %>%
  arrange(desc(av_abundance)) %>%
  slice_head(n = 100) %>% pull(OTU)

# Ovis
data_ovis <- data %>% filter(Genus == "Ovis") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

top_ovis <- data_ovis %>% group_by(OTU) %>%
  summarise(av_abundance = mean(Abundance)) %>%
  arrange(desc(av_abundance)) %>%
  slice_head(n = 100) %>% pull(OTU)

# Sus
data_sus <- data %>% filter(Genus == "Sus") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

top_sus <- data_sus %>% group_by(OTU) %>%
  summarise(av_abundance = mean(Abundance)) %>%
  arrange(desc(av_abundance)) %>%
  slice_head(n = 100) %>% pull(OTU)

# Top taxa
top_taxa <- c(top_equus, top_ovis, top_sus) %>% unique

# Filter data frame to top taxa for MCMCglmm
data_top <- data %>% filter(OTU %in% top_taxa)

data_equus <- data_equus %>% filter(OTU %in% top_taxa)
data_ovis <- data_ovis %>% filter(OTU %in% top_taxa)
data_sus <- data_sus %>% filter(OTU %in% top_taxa)

## Human vs wild comparison
data_human <- phy_genus_clr %>% subset_samples(Domestication %in% c("human", "wild")) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    psmelt %>% select(OTU, Abundance, Sample, Species, Genus, Domestication) %>%
    filter(OTU %in% top_taxa) %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

data_human$Domestication <- factor(data_human$Domestication, levels = c("wild", "human"))

######################
#### RUN ANALYSES ####
######################

# Run for all mammal samples and then for each genus separately

#### ALL WILD VS DOMESTIC ####

formula <- as.formula("Abundance ~ OTU + Genus:OTU + Domestication:OTU")

mcmc_out_all <- mcmcglmm_wrapper(data_top, formula, "mcmcglmm_all", 1)

mcmc_res_all <- process_mcmcglmm_out(mcmc_out_all, "mcmcglmm_all")

#### Within Equus ####

formula <- as.formula("Abundance ~ OTU + Domestication:OTU")

mcmc_out_equus <- mcmcglmm_wrapper(data_equus, formula, "mcmcglmm_equus", 1)

mcmc_res_equus <- process_mcmcglmm_out(mcmc_out_equus, "mcmcglmm_equus")

#### Within Ovis ####

mcmc_out_ovis <- mcmcglmm_wrapper(data_ovis, formula, "mcmcglmm_ovis", 1)

mcmc_res_ovis <- process_mcmcglmm_out(mcmc_out_ovis, "mcmcglmm_ovis")

#### Within Sus ####

mcmc_out_sus <- mcmcglmm_wrapper(data_sus, formula, "mcmcglmm_sus", 1)

mcmc_res_sus <- process_mcmcglmm_out(mcmc_out_sus, "mcmcglmm_sus")

#### Human vs Wild ####

formula <- as.formula("Abundance ~ OTU + Domestication:OTU")

mcmc_out_human <- mcmcglmm_wrapper(data_human, formula, "mcmcglmm_human", 1)

mcmc_res_human <- process_mcmcglmm_out(mcmc_out_human, "mcmcglmm_human")

#########################
#### COMBINE RESULTS ####
#########################

mcmc_res_all$dataset <- "Domestic vs Wild (All)"
mcmc_res_equus$dataset <- "Horse vs Zebra"
mcmc_res_ovis$dataset <- "Sheep vs Argali"
mcmc_res_sus$dataset <- "Pig vs Boar"
mcmc_res_human$dataset <- "Human vs Wild"

# Combine results
res <- bind_rows(mcmc_res_all, mcmc_res_equus, mcmc_res_ovis, mcmc_res_sus, mcmc_res_human) %>%
        filter(!grepl("Genus.*", term)) %>%
        group_by(dataset) %>%
        # Adjust pMCMC for multiple testing per dataset
        mutate(pMCMC_adj = p.adjust(pMCMC, method = "BH")) %>% ungroup %>%
        mutate(term = str_remove(term, "Domestication")) %>%
        arrange(dataset, term) %>%
        rename(coeff = post.mean)

write.csv(res, file.path(subdir, "mcmcglmm_res_combined.csv"), quote = FALSE, row.names = FALSE)

# Keep only taxa that are differentially abundant in any of the analyses
signif_taxa <- res %>% filter(pMCMC_adj < 0.05) %>% pull(OTU) %>% unique

# Order by difference of abundance in the largest subset
taxa_order <- res %>% filter(OTU %in% signif_taxa) %>% filter(dataset == "Domestic vs Wild (All)") %>%
               arrange(coeff) %>% pull(OTU)

res_filt <- res %>% filter(OTU %in% signif_taxa) %>%
            mutate(OTU = factor(OTU, levels = taxa_order)) %>%
            mutate(signif = case_when(pMCMC_adj < 0.001 ~ "***",
                                      pMCMC_adj < 0.01 ~ "**",
                                      pMCMC_adj < 0.05 ~ "*",
                                      pMCMC_adj < 0.1 ~ ".",
                                      TRUE ~ "")) %>%
            # Make x labels look a bit nicer
            mutate(dataset = gsub(" vs ", "\nvs ", dataset)) %>%
            mutate(dataset = factor(dataset, levels = c("Domestic\nvs Wild (All)", "Horse\nvs Zebra", "Sheep\nvs Argali", "Pig\nvs Boar",  "Human\nvs Wild")))

# Plot heatmap
p <- ggplot(data = res_filt, aes(x = dataset, y = OTU, fill = coeff)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "taxon domestic coefficient", na.value = "transparent") +
    geom_text(aes(label = signif), color = "black", size = 3) +
    theme(legend.position = "top", legend.text = element_text(angle = 90, vjust = 0.5, size = 10),
          legend.title = element_text(size = 10),
          panel.background = element_rect(fill = "grey90"), panel.grid = element_blank(),
          axis.text.y = element_text(size = 6), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) +
    guides(fill = guide_colorbar(barwidth = unit(2, "cm"), barheight = unit(0.5, "cm")))

ggsave(p, filename = file.path(subdir, "mcmcglmm_res_heatmap.png"), width = 4, height = 7)

# Also plot coeffs as scatterplots to show if they correlate
res_filt$dataset <- factor(gsub("\n", " ", res_filt$dataset), level = gsub("\n", " ", levels(res_filt$dataset)))

data_palette <- darken(species_palette[c("Equus quagga", "Ovis ammon", "Sus scrofa", "Homo sapiens")])
names(data_palette) <- c("Horse vs Zebra", "Sheep vs Argali", "Pig vs Boar",  "Human vs Wild")

# Maintain taxon order
taxon_order = levels(res_filt$OTU)

p <- filter(res_filt, c(!dataset %in% c("Domestic vs Wild (All)", "Human vs Wild"))) %>%
    mutate(taxon = factor(OTU, levels = taxon_order)) %>%
    ggplot(aes(y = taxon, x = coeff)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = dataset, shape = dataset), size = 3, alpha = 0.8) +
    geom_line(aes(group = taxon), linetype = "dotted", size = 0.3) +
    scale_colour_manual(values = data_palette, name = "") +
    scale_shape_manual(values = c(0, 2, 4), name = "") +
    theme(legend.position = "top", legend.direction = "vertical", axis.text.y = element_text(size = 6)) +
    xlab("OTU - domestic\ncoefficient")

ggsave(p, filename = file.path(subdir, "mcmcglmm_coeff_comparison.png"), width = 4, height = 7)

# Correlation test of coeffs between datasets
res_wide <- res %>%
    mutate(dataset = str_replace(dataset, "\n", " ")) %>%
    pivot_wider(id_cols = "OTU", names_from = dataset, values_from = coeff)

coeff_correlations = data.frame(
    Dataset1 = character(),
    Dataset2 = character(),
    Correlation = numeric(),
    P_value = numeric(),
    samples = integer()
)

datasets = res_filt$dataset %>% levels %>% str_replace_all(., "\n", " ")
datasets = datasets[datasets != "Domestic vs Wild (All)"]

for (d1 in datasets) {
  for (d2 in datasets) {
      # Keep only data for which there are coeff estimates in both datasets
      data <- res_wide[ , c("OTU", d1, d2)]
      data <- data[complete.cases(data), ]
      cases = data$OTU %>% unique %>% length
      cor_test <- cor.test(x = data[[d1]], y = data[[d2]])
      coeff_correlations <- rbind(coeff_correlations, data.frame(
          Dataset1 = d1,
          Dataset2 = d2,
          Correlation = cor_test$estimate,
          P_value = cor_test$p.value,
          samples = cases)
      )
  }
}

coeff_correlations$signif <- case_when(coeff_correlations$P_value < 0.001 ~ "***",
                                      coeff_correlations$P_value < 0.01 ~ "**",
                                      coeff_correlations$P_value < 0.05 ~ "*",
                                      coeff_correlations$P_value < 0.1 ~ ".",
                                      TRUE ~ "")

coeff_correlations$Dataset1 <- factor(coeff_correlations$Dataset1, levels = datasets)
coeff_correlations$Dataset2 <- factor(coeff_correlations$Dataset2, levels = rev(datasets))

write.csv(coeff_correlations, file.path(subdir, "mcmcglmm_coeff_correlations.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(aes(x = Dataset1, y = Dataset2, fill = Correlation), data = filter(coeff_correlations, Dataset1!=Dataset2)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Coefficient Correlation", na.value = "transparent") +
    geom_text(aes(label = signif), color = "black", size = 5) +
    theme(axis.title = element_blank(), legend.position = "top")

ggsave(p, filename = file.path(subdir, "mcmcglmm_coeff_correlation_heatmap.png"), width = 6, height = 6)

#########################
#### PLOT ABUNDANCES ####
#########################

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_genus_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Abundance") %>%
              left_join(rownames_to_column(select(data.frame(phy_genus_clr@sam_data), Common.name, Species, Genus, Domestication), "Sample"), by = "Sample") %>%
              mutate(Common.name = case_when(Domestication == "feral" ~ "Feral sheep",
                                          TRUE ~ Common.name))

res_abund <- res_filt %>% filter(pMCMC_adj < 0.05) %>%
            mutate(association = case_when(coeff < 0 ~ "wild+",
                                           coeff > 0 ~ "dom+")) %>%
            mutate(dataset = paste("comparison with", str_remove_all(dataset, "\n.*"))) %>%
            mutate(label = paste(association, dataset)) %>%
            # Combine labels for different results from the same OTU
            group_by(OTU) %>%
            summarise(label = paste(label, collapse="\n")) %>%
            # Add abundances
            left_join(select(abundances, OTU, Abundance, Species, Common.name, Genus, Domestication), by = c("OTU"))

# Reorder host species
species_levels <- abundances %>% data.frame %>% arrange(Genus, Domestication) %>% select(Species, Common.name) %>% unique
res_abund$Common.name <- factor(res_abund$Common.name, levels = species_levels$Common.name)

p <- ggplot(res_abund, aes(x = Common.name, y = Abundance, fill = Species, colour = Species)) +
    geom_boxplot(alpha = 0.8, size = 0.5, outliers = FALSE) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.8) +    scale_fill_manual(values = species_palette, name = "Species") +
    scale_colour_manual(values = darken(species_palette), name = "Species") +
    facet_wrap(~ paste(as.character(OTU), label, sep = "\n"), ncol = 5, scales = "free_y") +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "mcmcglmm_abundances.png"), width = 15, height = 30)
