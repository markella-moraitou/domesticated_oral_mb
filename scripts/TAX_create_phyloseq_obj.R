##### GENERATE PHYLOSEQ #####

#### Use CAT taxonomy, abundance based on contig mapping, and sample metadata
#### to generate a phyloseq object for further analysis

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(renv)
library(phyloseq)
library(stringr)
library(tibble)
library(microbiome)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with input data
outdir <- normalizePath(file.path("..", "output")) # Directory with output data
tabledir <- normalizePath(file.path(outdir, "Final_tables")) # Directory with taxon abundance table
subdir <- normalizePath(file.path("..", "output", "community_analysis")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(subdir, "phyloseq_objects")) # subdirectory for the phyloseq objects

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
dir.create(phydir, showWarnings = FALSE)

#######################
#####  LOAD INPUT #####
#######################

#### Define file paths
# OTU table
otu_table_path <- file.path(tabledir, "abundance_table_CAT.tsv")
tax_table_path = file.path(tabledir, "taxonomy_table_CAT.tsv")

# Sample metadata
metadata_path <- file.path(indir, "sample_metadata.csv")

rc_path <- file.path(outdir, "read_count.csv")
rl_path <- file.path(outdir, "read_length.csv")
decom_path <- file.path(outdir, "3_mapped_seqs", "decOM_output", "decOM_output.csv")

#### Load files
# Load OTU table
otu_table <- read.table(otu_table_path, sep="\t", comment.char="", header=TRUE)
tax_table <- read.table(tax_table_path, sep="\t", comment.char="", header=TRUE)

# Load metadata
metadata <- read.csv(metadata_path) # Sample metadata

rc <- read.csv(rc_path) %>% rename_with( ~ paste0(., "_count")) # read count per step

rl <- read.csv(rl_path) %>% rename_with( ~ paste0(., "_avlength")) # average read length per step

decom <- read.csv(decom_path, row.names = NULL) %>% # decOM output
  # remove entries where no kmers have been counted
  filter(rowSums(!is.na(select(., starts_with("p_")))) > 0)

###############################################
#### COMBINE AND TIDY UP SAMPLE METADATA  #####
###############################################

# Combine all sample and host species metadata in one big table
meta <- metadata

# Add suffix .A to sample IDs to match OTU table
decom <- decom %>% mutate(Sink = case_when(!grepl("_", Sink) ~ paste0(Sink, ".A"),
                        TRUE ~ Sink))

meta <-
  meta %>%
  left_join(rc, by=c("Sample.ID"="sample_count")) %>%
  left_join(rl, by=c("Sample.ID"="sample_avlength")) %>% 
  left_join(decom, by=c("Sample.ID"="Sink")) %>%
  as.data.frame %>%
  # Get genus 
  mutate(Genus = str_remove(Species, " .*"))

meta$Genus = factor(meta$Genus, levels = c("Equus", "Ovis", "Sus", "Homo"))
meta$Domestication = factor(meta$Domestication, levels = c("wild", "feral", "domestic", "human"))

spe_levels <- meta %>% arrange(Genus, Domestication) %>% pull(Species) %>% unique
meta$Species <- factor(meta$Species, levels=spe_levels)

# Distinguish between fully domestic and feral sheep
meta$Group <- ifelse(meta$Species == "Ovis aries", 
                    paste(meta$Domestication, "sheep", sep = " "),
                    str_to_lower(meta$Common.name))

group_levels <- meta %>% arrange(Genus, Domestication, Group) %>% pull(Group) %>% unique

meta$Group <- factor(meta$Group, levels=group_levels)

#Add column indicating samples and controls
meta$is.neg <- grepl("blank|control", meta$Species)

## Get better samples names
rename <- meta %>%
  select(Sample.ID, Species, Common.name, is.neg) %>% rename(old_name = Sample.ID) %>%
  # new names will consist of the first letter of the genus, the first three of the species epithet and a number
  group_by(Common.name) %>% mutate(num=row_number() %>% str_pad(width = 2, pad = "0")) %>%
  separate(col=Common.name, into=c("part1", "part2"), fill="left", sep=" ") %>%
  # Use first 4 letters of last adjective (part3) and the entire last word (part4)
  mutate(new_name=case_when(!is.neg ~ paste(str_to_lower(part2), num, sep="_"),
                          is.neg ~ paste(str_sub(part1, 1, 3), str_to_lower(part2), num, sep="_"))) %>%
  select(old_name, new_name)

# Temporary bit for this unknown library
rename$new_name[rename$old_name=="maybe_DM_017"] <- paste0(rename$new_name[rename$old_name == "DM_017"], ".maybe")

meta$new_name <- rename$new_name[match(meta$Sample.ID, rename$old_name)]

# Keep a version with all lab metadata (not just for samples in OTU table)
meta_all <- meta

################################
#### PREPARE TABLES FOR PS #####
################################

#### OTU TABLE ####

# May remove this later: decimals in abundance are throwing errors when calculating richness
# therefore I will round everything down to the next lower integer
# this means anything with abundance less than 1 will appear as absent (fine with that!)

tbl <- otu_table %>% column_to_rownames("lineage") %>%
  mutate(across(everything(), floor))

# Get new names
tbl <- tbl %>% 
  rename_with(~rename$new_name[match(., rename$old_name)], everything())

OTU = otu_table(tbl, taxa_are_rows = TRUE)

#### SAMPLE DATA ####
# Keep only samples that exist in tbl
meta <- meta %>% filter(new_name %in% colnames(tbl))

# Get rownames
rownames(meta) <- meta$new_name

SAM = sample_data(meta)

#### TAX TABLE ####

taxonomy <- tax_table %>% select(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "lineage")) %>%
          filter(lineage %in% taxa_names(OTU)) %>%
          mutate(phylum = str_remove(phylum, "_[A-Z]+")) %>% # to match palettes
          as.matrix()

# Get taxa names as row names
row.names(taxonomy) <- taxonomy[, "lineage"]

TAX = tax_table(taxonomy) %>% unique

###################################
####  CREATE PHYLOSEQ TABLES  #####
###################################

# Create phyloseq object and filter out empty samples
cat("Creating phyloseq table\n")
phy <- phyloseq(OTU, SAM, TAX)

# To reduce the size of this table and make it more manageable, remove taxa with less than 500 abundance
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- subset_samples(phy, sample_sums(phy) > 0)

# Agglomerate to species level and keep only prokaryotes and archaea
cat("Agglomerating to species level\n")
phy_sp <- phy %>% tax_glom(taxrank="species") %>% subset_taxa(species != "no support") %>%
 subset_taxa(superkingdom %in% c("Bacteria", "Archaea"))

# Extra lineage column causes issues with some scripts
phy_sp@tax_table <- phy_sp@tax_table[,c("superkingdom", "phylum", "class", "order", "family", "genus", "species")]

phy_sp <- subset_taxa(phy_sp, taxa_sums(phy_sp) > 0)
phy_sp <- subset_samples(phy_sp, sample_sums(phy_sp) > 0)

# Get taxa names instead of ids
taxa_names(phy_sp) <- make.unique(as.vector(phy_sp@tax_table[,"species"]))

#### Collect some metadata

# Collect number of OTUs per sample
phy_sp@sam_data$taxa_raw <- estimate_richness(phy_sp, measures="Observed")$Observed

# Calculate oral to soil ratio according to DecOM results
phy_sp@sam_data <- phy_sp@sam_data %>% data.frame %>% mutate(oral_to_soil_ratio=(p_mOral + p_aOral)/p_Sediment.Soil) %>% sample_data

# CLR-normalisation
phy_sp_clr <- phy_sp %>% transform('clr')

#####################
#### SAVE OUTPUT ####
#####################

# Phyloseq objects
saveRDS(phy_sp, file.path(phydir, "phy_sp.RDS"))
saveRDS(phy_sp_clr, file.path(phydir, "phy_sp_clr.RDS"))

# Constituent parts of phyloseq object
write.table(otu_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_OTU.tsv"), sep = "\t", row.names=TRUE, quote=FALSE)
write.table(tax_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_TAX.tsv"), sep = "\t", row.names=TRUE, quote=FALSE)
write.table(data.frame(sample_data(phy_sp)), file.path(subdir, "phyloseq_objects", "phy_sp_SAM.tsv"), sep = "\t", row.names=TRUE, quote=TRUE)
write.table(otu_table(phy_sp_clr), file.path(subdir, "phyloseq_objects", "phy_sp_clr_OTU.tsv"), sep = "\t", row.names=TRUE, quote=TRUE)

# Entire taxonomy table
write.table(taxonomy, file.path(subdir, "taxonomy_all.tsv"), sep = "\t", row.names=FALSE, quote=TRUE)

# Entire metadata table
# Get a big metadata table with all lab samples
meta_all <- meta_all %>% left_join(data.frame(phy_sp@sam_data)) %>%
  mutate(is.neg=grepl("blank|control", Species))

write.table(meta_all, file.path(subdir, "metadata_all.tsv"), sep = "\t", row.names=FALSE, quote=TRUE)
