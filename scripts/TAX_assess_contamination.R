##### ASSESS CONTAMINANTION #####

#### Use various sources of information to assess if a taxon is a contaminant or not
#### and if a sample is contaminated or not

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(tidyr)
library(stringr)
library(tibble)
library(scales)
library(cuperdec)
library(cowplot)
library(microbiomeutilities)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
tabledir <- normalizePath(file.path("..", "output", "Final_tables"))
outdir <- normalizePath(file.path("..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "assess_contamination")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it does not exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp <- readRDS(file.path(phydir, "phy_sp.RDS"))
phy_sp_clr <- readRDS(file.path(phydir, "phy_sp_clr.RDS"))

# Common contaminant list
common_contams <- read.table(file.path(indir, "common_contaminants.csv"), sep = ",", header = TRUE)

# Taxa identified to be overabundant in controls in previous analysis
previous_contams <- read.csv(file.path(indir, "abundant_in_negs_taxa.csv"))

# Taxonomy table
tax <- read.table(file.path(tabledir, "taxonomy_table_CAT.tsv"), sep = "\t", header = TRUE)

# Pydamage table
pydamage_summary <- read.table(file.path(tabledir, "pydamage_summary_CAT.tsv"), sep = "\t", header = TRUE)

### Set up thresholds for filtering

# Set minimum sample content threshold
min_samp <- 10^4

# Set minimum relative abundance detection threshold
min_ab <- 10^-4

# Set minimum prevalence threshold (in at least one host species must pass it)
prev_thresh <- 0.2 

##############################
#### COLLECT INFO ON TAXA ####
##############################

#### RELATIVE ABUNDANCE IN SAMPLES VS CONTROLS ####
# Get relative abundances
phy_sp_r <- transform(phy_sp, "compositional")
phy_sp_m <- phy_sp_r %>% psmelt

# Taxa present in nonhuman mammals
taxa <- phy_sp_m %>% filter(Species != "Homo sapiens" & Abundance > 0) 
phy_sp_m <- phy_sp_m %>% filter(OTU %in% taxa)

##### PREVALENCE AND AVERAGE RELATIVE ABUNDANCE ####
# Get prevalence of OTU in samples, controls and blanks

abundance_df <-
  phy_sp_m %>%
  # Do a weighted average to account for different number of samples per species
  group_by(Species, OTU) %>% summarise(mean_abundance = mean(Abundance)) %>% ungroup %>%
  group_by(OTU) %>% summarise(mean_abundance = mean(mean_abundance)) %>% ungroup %>%
  arrange(mean_abundance) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU)))

# Order taxa by average relative abundance
otu_order <- abundance_df$OTU %>% as.character

# Plot
p_ma <- ggplot(abundance_df, aes(x = mean_abundance, fill = mean_abundance, y = OTU)) +
  geom_point(shape = 21) +
  scale_colour_viridis_c() +
  scale_x_log10() + xlab("mean abundance\nin samples") +
  # Add threshold line
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x.top = element_text())

# Identify low prevalence taxa (anything that doesn't exist in more than 20% of samples in least one species)
species <- phy_sp@sam_data$Species %>% levels
prevalence <- data.frame(row.names = species)

for (spe in species) {
  subset <- phy_sp %>% subset_samples(Species == spe)
  temp <- prevalence(subset, detection = min_ab, include.lowest = TRUE, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("OTU")

write.csv(prevalence, file.path(subdir, "taxon_prevalence_per_species.csv"), quote = FALSE, row.names = FALSE)

prevalence_summ <- prevalence %>%
  pivot_longer(cols = -OTU, names_to = "Species", values_to = "prevalence") %>%
  # Keep only DC samples
  filter(!grepl("control", Species) & !grepl("blank", Species)) %>%
  mutate(OTU = factor(OTU, levels = otu_order)) %>%
  group_by(OTU) %>% summarise(mean = mean(prevalence, na.rm = TRUE),
            q1 = quantile(prevalence, 0.25, na.rm = TRUE),
            q3 = quantile(prevalence, 0.75, na.rm = TRUE),
            min = min(prevalence, na.rm = TRUE),
            max = max(prevalence, na.rm = TRUE))

p_p <- ggplot(prevalence_summ, aes(x = mean, fill = mean, y = OTU)) +
  geom_errorbar(aes(xmin = q1, xmax = q3), linewidth = 0.5, colour = "grey") +
  geom_point(aes(colour = mean, alpha = 0.5)) +
  scale_colour_viridis_c(option = "magma") + xlab("\nprevalence per\nhost species") +
  geom_vline(xintercept = prev_thresh, linetype = "dashed") +
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), , axis.title.x.top = element_text())

#### DAMAGE PATTERNS ####

# Match lineage to species names
name_lineage_match <- tax %>% select(lineage, species)

damage_df <- full_join(name_lineage_match, pydamage_summary, by = "lineage") %>%
          rename(OTU = species) %>%
          mutate(OTU = factor(OTU, levels = otu_order)) %>%
          filter(!is.na(OTU)) %>% 
          # Where more than one lineage is represented as one species, get median, min and max across lineages
          group_by(OTU) %>% summarise(median = median(median, na.rm = TRUE),
                                      min = min(min),
                                      max = max(max), .groups = "drop")

# Plot
p_d <- ggplot(damage_df, aes(x = median + 0.01, y = OTU)) +
  geom_point(aes(colour = median + 0.01), alpha = 0.8, size = 0.5,
            position = position_jitter(width = 0.05)) +
  scale_x_continuous(trans = "log10") +
  scale_colour_viridis_c(option = "turbo", trans = "log10") + xlab("\ndamage patterns\n(damage_model_pmax)") +
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x.top = element_text())

#### Previous information: is that taxon in common contaminant lists
#### Or identified in our previous study as contaminant?

contam_list <- data.frame(OTU = otu_order,
                          genus = str_remove(otu_order, " .*"),
                          previous_study = (otu_order %in% previous_contams$x))

contam_list <- contam_list %>% left_join(common_contams, by = c("genus" = "Contaminant_genera")) %>%
        mutate(Weyrich2019 = (Source %in% c("Weyrich et al 2019, Table 1", "Both studies")),
               Salter2014 = (Source %in% c("Salter et al 2014, Table 1", "Both studies")),
               Source = NULL) %>%
        mutate(OTU = factor(OTU, levels = otu_order))

p_c <- ggplot(pivot_longer(contam_list, cols = c(Weyrich2019, Salter2014, previous_study),
                          names_to = "Source", values_to = "Present"),
            aes(x = Source, y = OTU)) +
     geom_tile(aes(fill = Present)) +
     scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "white"), name = "Present in contaminant list") +
     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
           axis.title.y = element_blank(), axis.title.x = element_blank(),
           legend.position = "top")

#### Combine all tables ####

# Combine sample/negative abundance ratio, prevalence, mean abundance and habitat info
assess_taxa <- abundance_df %>%
               left_join(rename(select(damage_df, c(OTU, median)), p_damage_max_median = median)) %>%
               left_join(rename(select(prevalence_summ, c(OTU, max)), prevalence_max = max)) %>%
               left_join(contam_list, by = "OTU")

# Save table
write.table(assess_taxa, file=file.path(subdir, "assess_taxa.csv"), sep=",", row.names=FALSE, quote=FALSE)

#### Plot ####
p <- plot_grid(p_ma + theme(legend.position="none"),
               p_p + theme(legend.position="none"),
               p_d + theme(legend.position="none"),
               p_c + theme(legend.position="none"),
               nrow = 1, align = "h", axis = "tb", rel_widths = c(0.75, 0.75, 1))

ggsave(file=file.path(subdir, "assess_taxa.png"), p, width=10, height=12)

###################################
#### STRUCTURE OF CONTAMINANTS ####
###################################

# Subset to contaminants
contams <- previous_contams$x
phy_contams <- phy_sp %>% prune_taxa(contams, .)

# Remove empty samples
phy_contams <- prune_samples(sample_sums(phy_contams) > 0, phy_contams) %>%
  subset_samples(!grepl("blank", Species))

# Run PCoA on jaccard distance
pcoa_samples <- ordinate(phy_contams, method = "PCoA", distance = "jaccard")

p_contams <- plot_ordination(phy_contams, pcoa_samples, color = "Species") +
               geom_point(size = 3, alpha = 0.7) +
               scale_color_manual(values = species_palette, name = "Host Species") +
               theme(legend.position = "right")

ggsave(file=file.path(subdir, "contaminant_taxa_pcoa.png"), p_contams, width=8, height=6)

#### Do some samples have more contaminants identified than others?
contam_prop <- phy_sp_m %>%
            mutate(identified.contam = (OTU %in% contams)) %>%
            group_by(new_name, identified.contam, Species) %>%
            # Calculate abundance and number of taxa that are identified as contaminants, per sample
            summarise(prop_abundance = sum(Abundance),
                      n_taxa = n_distinct(OTU[Abundance > 0])) %>%
            # Calculate proportion of taxa identified as contaminants per sample
            ungroup %>% group_by(new_name) %>%
            mutate(prop_taxa = n_taxa/sum(n_taxa))

write.csv(contam_prop, file.path(subdir, "contaminant_proportions_per_sample.csv"), quote = FALSE, row.names = FALSE)

# Plot proportion of abundance and number of taxa identified as contaminants per sample

p1 <- ggplot(contam_prop, aes(y = new_name, x = prop_abundance, fill = identified.contam)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("TRUE" = "#FF5733", "FALSE" = "grey"), name = "Contaminant") +
      facet_wrap(~Species, ncol = 4, scales = "free_y") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
      labs(title = "A) Relative abundance of contaminant vs non-contaminant taxa per sample")

p2 <- ggplot(contam_prop, aes(y = new_name, x = n_taxa, fill = identified.contam)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("TRUE" = "#FF5733", "FALSE" = "grey"), name = "Contaminant") +
      facet_wrap(~Species, ncol = 4, scales = "free_y") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
      labs(title = "B) Number of contaminant vs non-contaminant taxa per sample")

p <- plot_grid(p1 + theme(legend.position="none"),
               p2 + theme(legend.position="none"),
               nrow = 2, align = "v", axis = "lr")

ggsave(file=file.path(subdir, "contaminant_proportions_per_sample.png"), p, width=10, height=20)

#################################
#### COLLECT INFO ON SAMPLES ####
#################################

# Have not run this yet
if (FALSE) {
 
# Use different sources of information to assess if a sample is contaminated or not

#### DECOM ####
# Get decom results
decom_tbl <- data.frame(phy_sp@sam_data) %>%
  select(new_name, Common.name, Species, starts_with("p_")) %>%
  # turn all NAs to 0 (for plotting)
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  # Calculate oral to soil+skin ratio
  mutate(oral_contam_ratio = (p_aOral + p_mOral)/(p_Sediment.Soil + p_Skin))

decom_tbl_long <- decom_tbl %>%
  # Pivot longer
  pivot_longer(starts_with("p_"), names_to = "Source", values_to = "Proportion")

# Turn source into a factor
decom_tbl_long$Source <- factor(decom_tbl_long$Source,
                                levels = c("p_aOral", "p_mOral", "p_Sediment.Soil", "p_Skin", "p_Unknown"))

# Order samples by oral proportion
sample_levels <- decom_tbl %>%
  arrange(Order_grouped, Species, desc(p_aOral + p_mOral)) %>% pull(new_name)

decom_tbl_long$new_name <- factor(decom_tbl_long$new_name, levels = sample_levels)
decom_tbl$new_name <- factor(decom_tbl$new_name, levels = sample_levels)

# Rearrange
decom_tbl_long <- decom_tbl_long %>% arrange(Order_grouped, new_name, Source)

# Palette for decom barplot
source_palette <- list("p_aOral"="#E54457", "p_mOral"="#AB0A1D", "p_Sediment.Soil"="#B2AD0B", "p_Skin"="#FFCC73", "p_Unknown"="#AAAAAA")

# Plot decom results
p_d <-
  ggplot(data = decom_tbl_long, aes(y = new_name, x = Proportion, fill = Source, group = Common.name)) +
  geom_bar(stat = "identity", colour = NA) +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_fill_manual(values = source_palette, name = "",
                    labels = c("ancient oral", "modern oral", "sediment/soil", "skin", "unknown")) +
  scale_x_continuous(expand = c(0,0)) + 
  theme(legend.position = "top",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        plot.margin = margin(t=1, r=10, b=1, l=1)) + 
  guides(fill = guide_legend(ncol = 2))

# Plot oral to contam ratio
p_r <- ggplot(data = decom_tbl, aes(x = oral_contam_ratio, fill = oral_contam_ratio, y = new_name)) +
  geom_point(shape = 21) +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_colour_viridis_c() +
  theme(legend.position="none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  # Add line at ratio = 1
  geom_vline(xintercept = 1, linetype = "dashed") +
  xlab("oral to\nsoil+skin\nratio") +
  scale_x_log10(breaks = c(0.1, 1, 10))

#### Total sample content (taxa sums) ####

# Taxa sums per sample
content <- data.frame(total_abundance = sample_sums(phy_sp),
                      new_name = sample_names(phy_sp_clr),
                      Order_grouped = phy_sp@sam_data$Order_grouped)

p_c <- ggplot(content, aes(x = total_abundance, y = new_name)) +
  geom_bar(stat = "identity") +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  theme(legend.position="none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  # Add line at threshold
  geom_vline(xintercept = min_samp, linetype = "dashed") +
  xlab("classified\nreads") +
  scale_x_log10(breaks = c(10^2, 10^4, 10^6))

# Show species as a bar
species_bar <- decom_tbl %>% select(new_name, Species, Common.name, Order_grouped) %>% unique %>%
               arrange(new_name) %>%
               # Choose one label per species (in the middle of the species group)
               group_by(Common.name, Species) %>% mutate(order = row_number()) %>%
               mutate(label = ifelse(order == floor(mean(order)), Common.name, "")) %>%
               mutate(label = ifelse(is.na(label), as.character(Species), label))

p_bar <-
  ggplot(data = species_bar, aes(y = new_name, x=1, fill = Species, group = Common.name)) +
  geom_tile() + scale_fill_manual(values = species_palette, name = "") +
  facet_grid(rows = vars(Order_grouped), scales = "free", space = "free", switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank()) +
        scale_y_discrete(position = "right",
                        label = setNames(species_bar$label, species_bar$new_name)) +
    xlab("") + ylab("")

#### Run cuperdec ####
## Run cuperdec to identify good and bad samples

# Get the OTU table
otu_table <- phy_sp@otu_table %>% data.frame %>% rownames_to_column(var="Taxon")
cpdc_taxa <- load_taxa_table(otu_table) 

# Get isolation sources using the habitat table
# If it is reported as oral more frequently than environmental (soil + marine) it is considered oral
isol <- habitats %>% filter(habitat %in% c("oral", "marine", "soil")) %>%
  pivot_wider(names_from = "habitat", values_from = "occurences", id_cols = "taxon", values_fill = 0) %>%
  group_by(taxon) %>% filter(oral > (soil+marine)) %>%
  select("taxon") %>% unique %>% rename("Taxon" = "taxon") %>%
  mutate(Isolation_Source="oral")

cpdc_db <- load_database(isol, target = "oral")

# Metadata map
metadata_map <- load_map(data.frame(phy_sp@sam_data),
                         sample_col = "new_name",
                         source_col = "Order")

# Run cuperdec
curves <- calculate_curve(cpdc_taxa, database = cpdc_db)

filter_result <- simple_filter(curves, percent_threshold = 50)

cpdc_plot <- plot_cuperdec(curves, metadata_map, filter_result) + theme(plot.background=element_rect(fill="white"))

ggsave(file=file.path(subdir, "cuperdec_plot.png"), cpdc_plot, width=16, height=8)

# Add cuperdec info to main plot
cpdc_label <- data.frame(new_name = sample_levels,
                         Passed = filter_result$Passed[match(sample_levels, filter_result$Sample)]) %>%
              mutate(label = ifelse(Passed, "*", ""), oral_contam_ratio = 100) %>%
              mutate(Order_grouped = phy_sp@sam_data$Order_grouped[match(new_name, phy_sp@sam_data$new_name)])

p_r_cpdc <- p_r + geom_text(data = cpdc_label, aes(label = label), size = 3, colour = "#AB0A1D")

# Combine plots
p <- plot_grid(p_d, p_r_cpdc, p_c, p_bar, ncol = 4, align = "h", axis = "tb", rel_widths = c(1, 0.3, 0.3, 0.5))

ggsave(file=file.path(subdir, "assess_samples.png"), p, width=10, height=15)

### Combine tables
assess_samples <- decom_tbl %>% left_join(select(content, c(new_name, total_abundance)), by = "new_name")
assess_samples$passed_cuperdec <- filter_result$Passed[match(assess_samples$new_name, filter_result$Sample)]
assess_samples$passed_min_reads <- assess_samples$total_abundance > min_samp
assess_samples$passed_oral_contam_ratio <- assess_samples$oral_contam_ratio > 1
# Sometimes there is no decOM info, in that case consult only read count
assess_samples$passed <- ifelse(!is.na(assess_samples$passed_oral_contam_ratio),
                                       assess_samples$passed_min_reads & assess_samples$passed_oral_contam_ratio,
                                       assess_samples$passed_min_reads)

# Save table
write.table(assess_samples, file=file.path(subdir, "assess_samples.csv"), sep=",", row.names=FALSE, quote=FALSE)
}
###################
#### FILTERING ####
###################

#### Filter samples ####
# Remove samples with an oral/contam ratio < 1 samples with less reads than the threshold

#contaminated_samples <- assess_samples %>% filter(!passed_oral_contam_ratio & !is.neg) %>% pull(new_name)
#shallow_samples <- assess_samples %>% filter(!passed_min_reads & !is.neg) %>% pull(new_name)

#write.csv(contaminated_samples, file = file.path(subdir, "contaminated_samples.csv"), quote = FALSE, row.names = FALSE)
#write.csv(shallow_samples, file = file.path(subdir, "shallow_samples.csv"), quote = FALSE, row.names = FALSE)

shallow_samples <- which(sample_sums(phy_sp) < min_samp) %>% names

phy_sp_f <- prune_samples(!(sample_names(phy_sp) %in% shallow_samples), phy_sp)
#phy_sp_f <- prune_samples(!(sample_names(phy_sp_f) %in% contaminated_samples), phy_sp_f)

#### Filter taxa ####

# Relative abundance filtering: Relative abundance under the threshold turned to 0
otu_rel_ab <- transform_sample_counts(phy_sp_f, function(x) x/sum(x)) %>% otu_table
phy_sp_f@otu_table[otu_rel_ab < min_ab] <- 0

# Remove empty taxa
phy_sp_f <- prune_taxa(taxa_sums(phy_sp_f) > 0, phy_sp_f)

# Remove taxa with less than 20% prevalence in a single species
low_prevalence_taxa <- prevalence_summ %>%
  filter(max < prev_thresh) %>%
  pull(OTU)

write.csv(low_prevalence_taxa, file.path(subdir, "low_prevalence_taxa.csv"), quote = FALSE, row.names = FALSE)

phy_sp_f <- phy_sp_f %>% subset_taxa(!(taxa_names(phy_sp_f) %in% previous_contams$x))
phy_sp_f <- subset_taxa(phy_sp_f, !(taxa_names(phy_sp_f) %in% low_prevalence_taxa))

#### Collect info on decontaminated dataset ####
# Get number of taxa
phy_sp_f@sam_data$taxa_filt <- estimate_richness(phy_sp_f, measures="Observed")$Observed

# Remove empty taxa
phy_sp_f <- prune_taxa(taxa_sums(phy_sp_f) > 0, phy_sp_f)

# Save filtered phyloseq object
saveRDS(phy_sp_f, file.path(phydir, "phy_sp_f.RDS"))

#### Run cuperdec again ####
## Run cuperdec to identify good and bad samples

# Get the OTU table
otu_table <- phy_sp_f@otu_table %>% data.frame %>% rownames_to_column(var="Taxon")
cpdc_taxa <- load_taxa_table(otu_table) 

# Get isolation sources using the habitat table
# If it is reported as oral more frequently than environmental (soil + marine) it is considered oral
isol <- habitats %>% filter(habitat %in% c("oral", "marine", "soil")) %>%
  pivot_wider(names_from = "habitat", values_from = "occurences", id_cols = "taxon", values_fill = 0) %>%
  group_by(taxon) %>% filter(oral > (soil+marine)) %>%
  select("taxon") %>% unique %>% rename("Taxon" = "taxon") %>%
  mutate(Isolation_Source="oral")

cpdc_db <- load_database(isol, target = "oral")

# Metadata map
metadata_map <- load_map(data.frame(phy_sp_f@sam_data),
                         sample_col = "new_name",
                         source_col = "Order")

# Run cuperdec
curves <- calculate_curve(cpdc_taxa, database = cpdc_db)

filter_result <- simple_filter(curves, percent_threshold = 50)

cpdc_plot <- plot_cuperdec(curves, metadata_map, filter_result) + theme(plot.background=element_rect(fill="white"))

ggsave(file=file.path(subdir, "cuperdec_plot_after_filtering.png"), cpdc_plot, width=16, height=8)
