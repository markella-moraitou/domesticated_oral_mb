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
library(cowplot)
library(grid)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) 
abpdir <- normalizePath(file.path("..", "output", "F2_abp_proteins")) 
statsdir <- normalizePath(file.path("..", "output", "function"))
datadir <- normalizePath(file.path(statsdir, "data"))
subdir <- normalizePath(file.path(statsdir, "starch_metabolism")) # subdirectory for the output of this script

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
        select(Sample, Sample.ID, Species, Genus, Common.name, Group, Domestication, Total_abundance, contig_reads_count) %>%
        mutate(Group = factor(str_to_lower(Group), levels = c(str_to_lower(levels(Group)))))

#  Abp finding from HMMER and relevant contig mapping stats
abpA_contigs <- read.table(file.path(abpdir, "abpA_mapping_stats.txt"), header = TRUE, sep = "\t")

abpA_tax <- read.table(file.path(abpdir, "abpA_tax.txt"), header = TRUE, sep = "\t", quote = "", comment = "")

# This file is a bit annoying to load
abpA_hmmer_file <- readLines(file.path(abpdir, "abpA_hmmer_table.txt"))
abpA_hmmer_filtered <- abpA_hmmer_file[!grepl("#", abpA_hmmer_file)]
abpA_hmmer_filtered <- gsub("\\s\\s+", "\t", abpA_hmmer_filtered)
abpA_hmmer <- read.table(text = abpA_hmmer_filtered, header = FALSE, sep = "\t", quote = "", comment = "")
abpA_hmmer <- abpA_hmmer[,c(1, 5, 6, 7, 18)]

colnames <- gsub("# ", "", abpA_hmmer_file[2]) %>% gsub("\\s\\s+", "\t", .) %>% str_split_1("\t")
colnames(abpA_hmmer) <- c("target.name", "E.value", "score", "bias", "description")

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
  mutate(rel_abundance = mapped_reads / Total_abundance) %>%
  mutate(genus = str_remove(genus, "\\*$"))

write.table(amy_abundances,
          file = file.path(subdir, "amylase_gene_abundance_stratified.tsv"),
          sep = "\t", quote = TRUE, row.names = FALSE)

# Get median abundances per amylase gene
amy_median_abund <- amy_abundances %>%
  group_by(gene_name, amylase_type, database) %>%
  summarize(median_abundance = median(rel_abundance)) %>%
  ungroup() %>%
  arrange(desc(median_abundance)) %>% unique

write.table(amy_median_abund, file = file.path(subdir, "amylase_gene_median_abundances.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

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

amy_abundances <- amy_abundances %>%
  mutate(genus_grouped = ifelse(genus %in% top_genera | genus == "no support", genus, "other"))

genus_palette <- brewer.pal(n = length(top_genera), name = "Set1")
names(genus_palette) <- top_genera

genus_palette["no support"] <- "black"
genus_palette["other"] <- "grey50"

amy_abundances$genus_grouped <- factor(amy_abundances$genus_grouped,
                                        levels = c(top_genera, "other", "no support"))

# Get mean abundance per host group
amy_summ %>% 
  # Sum by sample
  group_by(Sample, Group) %>% summarise(sum = sum(rel_abundance*100)) %>%
  group_by(Group) %>% summarise(mean_ra = mean(sum))

#### Plot ####

# By amylase type
p_t <- ggplot(amy_abundances, aes(y = Sample, x = rel_abundance, fill =  forcats::fct_rev(amylase_type))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01), expand=c(0,0), limits = c(0, 0.006)) +
  labs(y = "Sample",
       x = "",
       fill = "Amylase type") +
  scale_fill_manual(values = c("alpha-amylases" = "#D36306", "gamma-amylases" = "#FF9F4F")) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        strip.text.y = element_blank(),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# By microbial genus
p_g <- ggplot(amy_abundances, aes(y = Sample, x = rel_abundance, fill = genus_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01), expand=c(0,0), limits = c(0, 0.006)) +
  scale_fill_manual(values = genus_palette) +
  labs(y = "",
       x = "",
       fill = "Microbial genus") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))

p <- plot_grid(p_t, p_g, align = "h", axis = "tb", ncol = 2, rel_widths = c(1, 1.2))

p <- ggdraw(p) + 
  draw_label("Relative abundance of amylase genes", fontface = "bold", size = 15, x = 0.5, y = 0.27)
             
ggsave(filename = file.path(subdir, "amylase_gene_abundance.png"), width = 10, height = 8)

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
        summarise(rel_abundance = sum(rel_abundance)*100) # Get as percentage

# Encode diet as herbivore or omnivore
amy_abundances_total <- amy_abundances_total %>%
      # Remove humans
      filter(Genus != "Homo") %>%
      mutate(Diet = case_when(Genus == "Sus" ~ "omnivore",
                              Genus %in% c("Equus", "Ovis") ~ "herbivore"))

# Linear model
model <- aov(rel_abundance ~ Diet*Domestication + contig_reads_count, data = amy_abundances_total)
res <- summary(model)[[1]]

shapiro.test(residuals(model))
summary(residuals(model))

diagn <- diagnose_lm(model)

ggsave(filename  =  file.path(subdir, "diagnostics_anova.png"), plot_grid(plotlist = diagn),
       width  =  5, height = 5)

write.csv(res, file = file.path(subdir, "anova_amylase_relabund.csv"), quote = FALSE)

model_signif <- aov(rel_abundance ~ Domestication, data = amy_abundances_total)

# Tukey's HSD
tukey <- do.call("rbind", TukeyHSD(model_signif)) %>% data.frame %>% filter(!is.na(diff))

write.csv(tukey, file = file.path(subdir, "tukey_amylase_relabund.csv"), quote = FALSE)

########################
#### ABP ABUNDANCES ####
########################

# Fix tax table 
abpA_tax <- abpA_tax %>% mutate(across(superkingdom:species, ~ str_remove(str_remove(.x, "^.__"), ": .*")))

abpA_table <- abpA_hmmer %>% select(target.name, E.value, description) %>%
  # Get contig name
  mutate(contigName = str_remove(str_remove(target.name, "clustered_contigs_rep_seq_"), "_[0-9]+$")) %>%
  # Combine with mapping stats
  left_join(abpA_contigs) %>% filter(!is.na(coverage)) %>%
  # Filter for identity and coverage as we do for the general gene table
  mutate(mapped_reads = case_when(identity > 98 & coverage > 50 ~ mapped_reads,
                                  TRUE ~ 0)) %>%
  # Add taxonomy of contig
  left_join(select(abpA_tax, c(X..contig, phylum, class, order, family, genus, species)), by = c("contigName" = "X..contig")) %>%
  # When classification are missing, use 'no support'
  mutate(across(phylum:species, ~ replace_na(.x, "no support"))) %>%
  # Combine with sample metadata
  rename(Sample.ID = sample) %>%
  full_join(meta) %>% 
  # Remove missing species (these are controls)
  filter(!is.na(Species)) %>% 
  mutate(mapped_reads = case_when(is.na(mapped_reads) ~ 0, TRUE ~ mapped_reads)) %>%
  # Calculate rel abundance
  mutate(rel_abundance = mapped_reads/Total_abundance) %>%
  # Simplify description
  mutate(description = case_when(grepl("rank: E", description) ~ "no annotation",
                                 TRUE ~ str_remove(str_remove(description, ".*rank: .; "), " \\(db=.*\\)")))

write.csv(abpA_table, file.path(subdir, "abpA_stats.csv"), quote = FALSE, row.names = FALSE)

# First plot relatively loose hits: E.value < 0.05
abpA_table_loose <- 
    abpA_table %>%
    mutate(rel_abundance = case_when(E.value < 0.05 ~ rel_abundance, TRUE ~ 0))

# Get top genera for plotting
top_genera2 <- abpA_table_loose %>% group_by(genus) %>% summarise(rel_abundance = sum(rel_abundance)) %>%
        arrange(desc(rel_abundance)) %>% filter(genus != "no support" & rel_abundance > 0) %>% head(3) %>% pull(genus)

# Match with previous plots
genus_palette2 <- genus_palette[top_genera2]
genus_palette2 <- genus_palette2[!is.na(genus_palette2)]

# Add remaining genera
missing_gen <- top_genera2[!top_genera2 %in% names(genus_palette2) & !top_genera2 %in% c("no support", "other")]
missing_col <- brewer.pal(n = length(missing_gen), name = "Set2")

genus_palette2 <- c(genus_palette2, setNames(missing_col, missing_gen))

genus_palette2["no support"] <- "black"
genus_palette2["other"] <- "grey50"

abpA_table_loose <- 
    abpA_table_loose %>%
    mutate(genus_grouped = case_when(genus %in% c(top_genera2, "no support") ~ genus, TRUE ~ "other")) %>%
    mutate(genus_grouped = factor(genus_grouped, levels = names(genus_palette2)))

# Plot
p1 <- ggplot(abpA_table_loose, aes(y = Sample.ID, x = rel_abundance, fill = genus_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.001)) +
  labs(y = "Sample",
       x = "Rel. abundance of loose abpA matches\n(E-value < 0.05)",
       fill = "Microbial genus") +
  scale_fill_manual(values = genus_palette2) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 1, byrow = FALSE))

ggsave(filename = file.path(subdir, "abpA_abundance_loose.png"), width = 5, height = 9)

# Stricter hits
abpA_table_filt <- abpA_table_loose %>%
  # Keep only hits with E-value < 10^-3
  mutate(rel_abundance = case_when(E.value < 10^-3 ~ rel_abundance, TRUE ~ 0))

# Plot
p <- ggplot(abpA_table_filt, aes(y = Sample.ID, x = rel_abundance, fill = genus_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.001)) +
  labs(y = "Sample",
       x = "Rel. abundance of strict abpA matches\n(E-value < 10^-3)",
       fill = "Microbial genus") +
  scale_fill_manual(values = genus_palette2) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 1, byrow = FALSE))

ggsave(filename = file.path(subdir, "abpA_abundance_strict.png"), width = 5, height = 9)
