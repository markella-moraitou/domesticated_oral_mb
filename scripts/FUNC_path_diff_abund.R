##### DIFFERENTIAL ABUNDANCE OF PATHWAYS #####

#### Run differential abundance analysis on pathway data

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
datadir <- normalizePath(file.path("..", "output", "function", "data"))
pathdir <- normalizePath(file.path("..", "output", "function", "pathway_completeness")) # Directory with pathway analysis output
subdir <- normalizePath(file.path("..", "output", "function", "path_diff_abundance")) # subdirectory for the output of this script

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

# Source MCMCglmm wrapper and output processing functions
source(file.path("modules", "mcmcglmm_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Load pathway phyloseq
for (phy_file in list.files(pathdir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(pathdir, phy_file)))
}

# Load gene phyloseq
for (phy_file in list.files(datadir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(datadir, phy_file)))
}

# pathway to ko
path_to_ko <- read.csv(file.path(pathdir, "pathways_to_kos.csv"))

#########################
#### PREP TEST INPUT ####
#########################

# Remove feral sheep and Humans because they are an ambiguous category
# If diff abund genera are found, I will plot the abundances to see
# if they match the domesticated or wild patterns better
phy_da <- phy_pathway_clr %>% subset_samples(Domestication != "feral" & Domestication != "human")
phy_da <- phy_da %>% prune_taxa(taxa_sums(phy_da) > 0, .)

phy_da@sam_data$Domestication <- factor(phy_da@sam_data$Domestication, levels = c("wild", "domestic"))

# Get as data frame
data <- psmelt(phy_da) %>%
        select(OTU, Abundance, Sample, Species, Genus, Domestication, unmapped_count, contig_reads_count)

## Create subsets per genus for the within-genus analyses

# Equus
data_equus <- data %>% filter(Genus == "Equus") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

# Ovis
data_ovis <- data %>% filter(Genus == "Ovis") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

# Sus
data_sus <- data %>% filter(Genus == "Sus") %>% group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

## Human vs wild comparison
data_human <- phy_pathway_clr %>% subset_samples(Domestication %in% c("human", "wild")) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    psmelt %>% select(OTU, Abundance, Sample, Species, Genus, Domestication) %>%
    group_by(OTU) %>% filter(sum(Abundance) > 0) %>% ungroup

data_human$Domestication <- factor(data_human$Domestication, levels = c("wild", "human"))

######################
#### RUN ANALYSES ####
######################

# Run for all mammal samples and then for each genus separately

#### ALL WILD VS DOMESTIC ####

formula <- as.formula("Abundance ~ OTU + Genus:OTU + Domestication:OTU")

mcmc_out_all <- mcmcglmm_wrapper(data, formula, "mcmcglmm_all", 1)

mcmc_res_all <- process_mcmcglmm_out(mcmc_out_all, "mcmcglmm_all")

mcmc_res_all <- mcmc_res_all %>% rename(pathway = OTU)

#### Within Equus ####

formula <- as.formula("Abundance ~ OTU + Domestication:OTU")

mcmc_out_equus <- mcmcglmm_wrapper(data_equus, formula, "mcmcglmm_equus", 1)

mcmc_res_equus <- process_mcmcglmm_out(mcmc_out_equus, "mcmcglmm_equus")

mcmc_res_equus <- mcmc_res_equus %>% rename(pathway = OTU)

#### Within Ovis ####

mcmc_out_ovis <- mcmcglmm_wrapper(data_ovis, formula, "mcmcglmm_ovis", 1)

mcmc_res_ovis <- process_mcmcglmm_out(mcmc_out_ovis, "mcmcglmm_ovis")

mcmc_res_ovis <- mcmc_res_ovis %>% rename(pathway = OTU)

#### Within Sus ####

mcmc_out_sus <- mcmcglmm_wrapper(data_sus, formula, "mcmcglmm_sus", 1)

mcmc_res_sus <- process_mcmcglmm_out(mcmc_out_sus, "mcmcglmm_sus")

mcmc_res_sus <- mcmc_res_sus %>% rename(pathway = OTU)

#### Human vs Wild ####

formula <- as.formula("Abundance ~ OTU + Domestication:OTU")

mcmc_out_human <- mcmcglmm_wrapper(data_human, formula, "mcmcglmm_human", 1)

mcmc_res_human <- process_mcmcglmm_out(mcmc_out_human, "mcmcglmm_human")

mcmc_res_human <- mcmc_res_human %>% rename(pathway = OTU)

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
        rename(coeff = post.mean) %>%
        left_join(rownames_to_column(data.frame(phy_pathway@tax_table), "pathway"))

write.csv(res, file.path(subdir, "mcmcglmm_res_func.csv"), quote = TRUE, row.names = FALSE)

# Keep only taxa that are differentially abundant in any of the wild vs dom analyses
signif_paths <- res %>% filter(pMCMC_adj < 0.05 & dataset != "Human vs Wild") %>% pull(path_name) %>% unique

# Order by difference of abundance in the largest subset
path_order <- res %>% filter(path_name %in% signif_paths) %>% filter(dataset == "Domestic vs Wild (All)") %>%
               arrange(coeff) %>% pull(path_name)

res_filt <- res %>% filter(path_name %in% signif_paths) %>%
            mutate(path_name = factor(path_name, levels = path_order)) %>%
            mutate(signif = case_when(pMCMC_adj < 0.001 ~ "***",
                                      pMCMC_adj < 0.01 ~ "**",
                                      pMCMC_adj < 0.05 ~ "*",
                                      pMCMC_adj < 0.1 ~ ".",
                                      TRUE ~ "")) %>%
            # Make x labels look a bit nicer
            mutate(dataset = gsub(" vs ", "\nvs ", dataset)) %>%
            mutate(dataset = factor(dataset, levels = c("Domestic\nvs Wild (All)", "Horse\nvs Zebra", "Sheep\nvs Argali", "Pig\nvs Boar",  "Human\nvs Wild")))

# Plot heatmap
p <- ggplot(data = res_filt, aes(x = dataset, y = path_name, fill = coeff)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "domestic\ncoefficient", na.value = "transparent") +
    geom_text(aes(label = signif), color = "black", size = 3) +
    theme(legend.position = "top", legend.text = element_text(angle = 90, vjust = 0.5, size = 10),
          legend.title = element_text(size = 10),
          panel.background = element_rect(fill = "grey90"), panel.grid = element_blank(),
          axis.text.y = element_text(size = 8), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) +
    guides(fill = guide_colorbar(barwidth = unit(2, "cm"), barheight = unit(0.5, "cm")))

ggsave(p, filename = file.path(subdir, "mcmcglmm_res_heatmap.png"), width = 5, height = 5)

# Also plot coeffs as scatterplots to show if they correlate
res_filt$dataset <- factor(gsub("\n", " ", res_filt$dataset), level = gsub("\n", " ", levels(res_filt$dataset)))
res_filt$signif <- ifelse(res_filt$pMCMC_adj < 0.05, "TRUE", "FALSE")

data_palette <- darken(species_palette[c("Equus quagga", "Ovis ammon", "Sus scrofa", "Homo sapiens")])
names(data_palette) <- c("Horse vs Zebra", "Sheep vs Argali", "Pig vs Boar",  "Human vs Wild")

p <- filter(res_filt, c(!dataset %in% c("Domestic vs Wild (All)", "Human vs Wild"))) %>%
    mutate(taxon = factor(path_name, levels = path_order)) %>%
    ggplot(aes(x = path_name, y = coeff)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = dataset, shape = dataset, fill = dataset, alpha = signif), size = 3) +
    geom_line(aes(group = pathway), linetype = "dotted", size = 0.3) +
    scale_colour_manual(values = data_palette, name = "") +
    scale_fill_manual(values = data_palette, name = "") +
    scale_shape_manual(values = c(22, 23, 25), name = "") +
    scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 0.8), name = "adj. pMCMC < 0.05") +
    theme(legend.position = "right", legend.direction = "vertical",
          axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) +
    ylab("domestic coefficient")

ggsave(p, filename = file.path(subdir, "mcmcglmm_coeff_comparison.png"), width = 10, height = 6)

#########################
#### PLOT ABUNDANCES ####
#########################

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_pathway_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Abundance") %>%
              # Fix pathway names
              mutate(OTU = str_replace(pattern = "path.", replacement = "path:", OTU)) %>%
              # Add sample and pathway info
              left_join(rownames_to_column(select(data.frame(phy_pathway_clr@sam_data), Common.name, Species, Group, Genus, Domestication), "Sample"), by = "Sample") %>%
              left_join(rownames_to_column(select(data.frame(phy_pathway@tax_table), path_name), "OTU"), by = "OTU")

mcmc_abund <- res_filt %>% filter(pMCMC_adj < 0.05) %>%
            mutate(association = case_when(coeff < 0 ~ "wild+",
                                           coeff > 0 ~ "wild-")) %>%
            mutate(dataset = str_remove_all(dataset, "\n.*")) %>%
            mutate(label = paste(association, dataset)) %>%
            # Combine labels for different results from the same OTU
            group_by(path_name) %>%
            summarise(label = paste(label, collapse="\n")) %>%
            # Add abundances
            left_join(select(abundances, path_name, Abundance, Species, Group, Common.name, Genus, Domestication))

# Reorder host species
group_levels <- levels(phy_pathway@sam_data$Group)
mcmc_abund$Group <- factor(mcmc_abund$Group, levels = group_levels)

p <- ggplot(mcmc_abund, aes(x = Group, y = Abundance, fill = Species, colour = Species)) +
    geom_boxplot(alpha = 0.8, size = 0.5, outliers = FALSE) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.8) +    scale_fill_manual(values = species_palette, name = "Species") +
    scale_colour_manual(values = darken(species_palette), name = "Species") +
    facet_wrap(~ paste(as.character(path_name), label, sep = "\n"), ncol = 4, scales = "free_y") +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "mcmcglmm_abundances.png"), width = 10, height = 15)

#################################
#### CHECK GENES & TAXONOMY  ####
#################################

# For differentially abundant pathways, check which KOs are involved
# And which taxa they're encoded by

# Find differentially abundant pathways between domestic and wild animals
da_kos <- res_filt %>% filter(dataset != "Human\nvs Wild") %>%
                  filter(pMCMC < 0.05) %>% pull(pathway) %>% unique

# Get genes for these KOs
da_genes <- path_to_ko %>% filter(pathway %in% da_kos) %>%
            left_join(rownames_to_column(data.frame(phy_pathway@tax_table), "pathway"))

# Get gene info and abundances
gene_abundances <- phy_gene_f %>%
  subset_taxa(taxa_names(phy_gene_f) %in% str_remove(da_genes$kos, "ko:")) %>%
  psmelt %>% 
  select(OTU, Sample.ID, Abundance, Species, Domestication, gene_name, gene_description) %>%
  mutate(kos = paste0("ko:", OTU)) %>%
  left_join(da_genes, by = "kos", relationship = "many-to-many")

# Plot abundances of genes per pathway and species
p <- gene_abundances %>% group_by(Species, Domestication, OTU, gene_name, gene_description, pathway, path_name) %>%
      summarise(median_abundance = median(Abundance) + 1) %>%
      ggplot(aes(y = OTU, x = median_abundance)) +
      geom_bar(stat = "identity", aes(fill = Species)) +
      scale_fill_manual(values = species_palette) +
      facet_grid(rows = vars(path_name), cols = vars(Species), scales = "free", space = "free_y", labeller = label_wrap_gen(width=20)) +
      scale_x_continuous(trans = 'log10') +
      theme(strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_blank(), legend.position = "none")

ggsave(p, filename = file.path(subdir, "da_pathway_genes.png"), width = 10, height = 20)
