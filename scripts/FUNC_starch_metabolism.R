##### Amylase gene abundance #####

#### Test if amylase gene abundance differs between domesticates and wild

################
#### SET UP ####
################

#### LOAD PACKAGES ####
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
indir <- normalizePath(file.path("..", "input")) 
outdir <- normalizePath(file.path("..", "output", "function"))
datadir <- normalizePath(file.path(outdir, "data"))
subdir <- normalizePath(file.path(outdir, "starch_metabolism")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Stratified sample data
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))

# Select sample metadata
meta <- data.frame(phy_gene_f@sam_data) %>% rownames_to_column("Sample") %>%
        select(Sample, Sample.ID, Species, Genus, Common.name, Group, Domestication, Total_abundance, contig_reads_count)

#################################
#### AMYLASE GENE ABUNDANCES ####
#################################

# Identify differences in amylase gene abundances
# From Mehta and Satyanarayana 2016:
amylases_df <- list("alpha-amylases" = "3.2.1.1",
                      "beta-amylases" = "3.2.1.2",
                      "gamma-amylases" = "3.2.1.3",
                      "maltogenetic amylases" = "3.2.1.133") %>% do.call(rbind, .) %>% data.frame %>%
                rownames_to_column("amylase_type") %>% rename("EC" = ".")

amy_abundances <- gene_str %>%
  filter(str_detect(gene_name, "3\\.2\\.1\\.[123]\\b")) %>%
  # Extract EC numbers from gene_description or name
  mutate(EC = str_extract(paste(gene_description, gene_name), "3.2.[0-9]+.[0-9]+")) %>%
  # Match the codes
  filter(EC %in% amylases_df$EC) %>%
  # Add info on category
  left_join(amylases_df) %>%
  # Add sample metadata
  left_join(meta, by = c("Sample" = "Sample.ID")) %>%
  filter(!is.na(Total_abundance)) %>%
  # Calculate relative abundance
  mutate(rel_abundance = mapped_reads / Total_abundance)

write.csv(amy_abundances,
          file = file.path(subdir, "amylase_gene_abundance_stratified.tsv"),
          quote = TRUE, row.names = FALSE)

# Get median abundances per amylase gene
amy_median_abund <- amy_abundances %>%
  group_by(gene_name) %>%
  summarize(median_abundance = median(rel_abundance)) %>%
  ungroup() %>%
  arrange(desc(median_abundance)) %>% unique

write.csv(amy_median_abund, file = file.path(subdir, "amylase_gene_median_abundances.tsv"), quote = TRUE, row.names = FALSE)

# Summarize amylase gene abundance per sample and microbial genus
amy_summ <- amy_abundances %>%
  group_by(Sample, genus, Species, Group, Common.name, Domestication) %>%
  # Sum all amylase mapped reads per sample and genus
  summarize(rel_abundance = sum(rel_abundance))

# Keep genera with largest total amylase mapped reads across all samples
# Group the rest as 'other'

top_genera <- amy_summ %>%
  group_by(genus) %>% filter(genus != "no support") %>%
  summarize(total_rel_abundance = sum(rel_abundance)) %>%
  ungroup() %>%
  arrange(desc(total_rel_abundance)) %>%
  slice_head(n = 8) %>%
  pull(genus)

amy_summ <- amy_summ %>%
  mutate(genus = ifelse(genus %in% top_genera | genus == "no support", genus, "other"))

genus_palette <- brewer.pal(n = length(top_genera), name = "Set1")
names(genus_palette) <- top_genera

genus_palette["no support"] <- "black"
genus_palette["other"] <- "grey50"

amy_summ$genus <- factor(amy_summ$genus,
                                    levels = c(top_genera, "other", "no support"))

#### Plot ####

# By amylase type
p_t <- ggplot(amy_abundances, aes(y = Sample, x = rel_abundance, fill =  forcats::fct_rev(amylase_type))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(y = "Sample",
       x = "Relative abundance of amylase genes",
       fill = "Amylase type") +
  scale_fill_manual(values = c("alpha-amylases" = "#D36306", "gamma-amylases" = "#FF9F4F")) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))

ggsave(filename = file.path(subdir, "amylase_gene_abundance_by_type.png"), width = 5, height = 8)

# By microbial genus
p_g <- ggplot(amy_abundances, aes(y = Sample, x = rel_abundance, fill = genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  scale_fill_manual(values = genus_palette) +
  labs(y = "Sample",
       x = "Relative abundance of amylase genes",
       fill = "Microbial genus") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 2, byrow = FALSE))

ggsave(filename = file.path(subdir, "amylase_gene_abundance_by_genus.png"), width = 5, height = 9)

#### TEST ####

# Function for lm diagnostics
diagnose_lm <- function(model) {
  diag_plots <- list() # Store all diagnostic plots here
 
  # Assumption: residuals shouldn't be correlated to fitted values
  tbl <- data.frame(fitted.values = fitted.values(model), residuals = residuals(model))
  ranf_fitted <- ggplot(data = tbl, aes(x = fitted.values, y = residuals)) + geom_point() + theme_bw() +
    labs(title = "Model residuals to fitted values")
  
  diag_plots[["ranef_fitted"]] <- ranf_fitted
  
  # Assumption: residuals are normal
  res_qqplot <- ggplot(data = data.frame(residuals = residuals(model)),
                         aes(sample = residuals)) +
    geom_qq() + geom_qq_line() + theme_bw() + xlab("Theoretical quantiles") + ylab("Observed quantiles") +
    labs(title = paste("Q-Q Plot for model residuals"))
  
  diag_plots[["res_qqplot"]] <- res_qqplot
  return(diag_plots)
}

# Sum abundances of all genes
amy_abundances_total <- amy_abundances %>% group_by(Sample, Genus, Domestication, contig_reads_count) %>%
        summarise(rel_abundance = sum(rel_abundance))

# Linear model
model <- aov(rel_abundance ~ Genus*Domestication + contig_reads_count, data = amy_abundances_total)
res <- summary(model)[[1]]

shapiro.test(residuals(model))
summary(residuals(model))

diagn <- diagnose_lm(model)

ggsave(filename  =  file.path(subdir, "diagnostics_anova.png"), plot_grid(plotlist = diagn),
       width  =  14, height = 14)

write.csv(res, file = file.path(subdir, "anova_amylase_relabund.csv"), quote = FALSE)

model_signif <- aov(rel_abundance ~ Domestication, data = amy_abundances_total)

# Tukey's HSD
tukey <- do.call("rbind", TukeyHSD(model_signif)) %>% data.frame %>% filter(!is.na(diff))

write.csv(tukey, file = file.path(subdir, "tukey_amylase_relabund.csv"), quote = FALSE)
