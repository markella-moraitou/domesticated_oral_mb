##### Topological MAG triplets analysis #####

#### Analysed MAGs to compare if the topologies are more consistent with ((wild, domestic),human)
#### or with ((domestic, human), wild)

#### LOAD PACKAGES ####
library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)
library(cowplot)
library(ggtree)
library(ggtreeExtra)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output")) # Directory with output
magdir <- normalizePath(file.path(outdir, "mags")) # Directory with MAG output
subdir <- normalizePath(file.path(magdir, "mag_triplets")) # subdirectory for the output of this script

if(!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

source(file.path("modules", "plot_setup.R"))
plot_setup(file.path("..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(magdir, "bac_meta_drep.tsv"), sep="\t", header=TRUE, quote = "", comment = "")
ar_meta <- read.table(file.path(magdir, "ar_meta_drep.tsv"), sep="\t", header=TRUE, quote = "", comment = "")

# MAG trees
bac_tree <- read.tree(file = file.path(magdir, "bac_tree_drep.tree"))
ar_tree <- read.tree(file = file.path(magdir, "ar_tree_drep.tree"))

presabs <- read.csv(file.path(magdir, "mag_mapping_stats", "hq_mag_presence_per_host.csv"))

############################
#### SEARCH TREE SPLITS ####
############################

# For this analysis we consider all sheep domestic
species_dom <- bac_meta %>% select(host_species, Domestication) %>%
                filter(!is.na(host_species) & Domestication != "feral") %>% unique

# Combine MAG taxonomy with presence/absence data
data <- rbind(bac_meta, ar_meta) %>% select(label, bin, domain, phylum, order, family, genus) %>%
    right_join(presabs, by = c("label" = "label", "bin" = "bin")) %>%
    mutate(Domestication = species_dom$Domestication[match(host_species, species_dom$host_species)]) %>%
    mutate(host = case_when(Domestication == "human" ~ "human",
                            TRUE ~ paste(Domestication, host_genus)))

# Write function that counts the number of species on each side of a split
split_balance <- function(tree, node, data) {
    # Get the two clades defined by the node and the number of tips they lead to
    clades <- c(tree$tip.label, tree$node.label)[Descendants(tree, node, "children")]
    # Get tip names
    clade1 <- as.character(clades[1])
    clade2 <- as.character(clades[2])
    # Get all tips in a clade
    tips1 <- c(tree$tip.label, tree$node.label)[Descendants(tree, clade1, "tips")[[1]]]
    tips2 <- c(tree$tip.label, tree$node.label)[Descendants(tree, clade2, "tips")[[1]]]
    # Get the number of host species in each clade
    hosts1 <- data %>% filter(label %in% tips1) %>%
                        pull(host) %>% unique %>% sort
    hosts2 <- data %>% filter(label %in% tips2) %>%
                        pull(host) %>% unique %>% sort
    # Collect into a data frame
    df <- data.frame(
        clade = c(clade1, clade2),
        n_hosts = c(length(hosts1), length(hosts2)),
        hosts = c(paste(hosts1, collapse=","), paste(hosts2, collapse=","))
    ) %>%
    # Arrange so that clade with most hosts is first
    arrange(desc(n_hosts)) %>% mutate(split = c(1, 2))
    return(df)
}

# Iterate over all internal nodes 
# and calculate the number of host species on each side of the split

## Bacteria
splits_df_bac <- data.frame()

node_list <- bac_tree$node.label

for (node in node_list) {
    df <- split_balance(bac_tree, node, data)
    df_wide <- df %>% 
        mutate(node = node) %>%
        pivot_wider(names_from = split, values_from = c(clade, n_hosts, hosts))
    splits_df_bac <- rbind(splits_df_bac, df_wide)
}

splits_df_bac$domain <- "Bacteria"

## Archaea
splits_df_ar <- data.frame()

node_list <- ar_tree$node.label

for (node in node_list) {
    df <- split_balance(ar_tree, node, data)
    df_wide <- df %>% 
        mutate(node = node) %>%
        pivot_wider(names_from = split, values_from = c(clade, n_hosts, hosts))
    splits_df_ar <- rbind(splits_df_ar, df_wide)
}

splits_df_ar$domain <- "Archaea"

splits_df <- rbind(splits_df_bac, splits_df_ar)

write.csv(splits_df, file = file.path(subdir, "mag_tree_splits.csv"), row.names=FALSE)

## Keep only splits where one side has 1 host species and the other two host species

splits_df_filt <- splits_df %>%
    filter((n_hosts_1 == 2 & n_hosts_2 == 1))

# Display the node type like ((A,B),C)
label_splits <- function(splits) {
    # Domestic-wild patterns to match
    dw_patterns <- c("Sus.*Sus",
                     "Ovis.*Ovis",
                     "Equus.*Equus") %>% paste(collapse="|")
    splits <- splits %>%
        mutate(all_hosts = paste0(hosts_1, ",", hosts_2)) %>%
        mutate(node_type = case_when(
            # Domestic and wild together of the same genus
            hosts_2 == "human" & grepl(dw_patterns, hosts_1) ~ "((D, W), H)",
            # Domestic and human together
            grepl("wild", hosts_2) & grepl("domestic.*human", hosts_1) & grepl(dw_patterns, all_hosts) ~ "((D, H), W)", 
            # Wild and human together as a control
            grepl("domestic", hosts_2) & grepl("human.*wild", hosts_1) & grepl(dw_patterns, all_hosts) ~ "((H, W), D)",
            TRUE ~ "other"
        ))
}

splits_df_filt <- label_splits(splits_df_filt)

write.csv(splits_df_filt, file = file.path(subdir, "mag_tree_splits_2_1.csv"), row.names=FALSE)

splits_summary <- splits_df_filt %>%
        group_by(node_type) %>%
        summarise(count = n(),
                  proportion = n() / nrow(splits_df_filt))

#######################
#### PERMUTATIONS #####
#######################

# Function to permute host labels and re-calculate split types
permute_splits <- function(bac_tree, ar_tree, data, n_permutations=1000) {
    permuted_results <- data.frame()
    host_perm_map <- data.frame(original = unique(data$host))
    for (i in 1:n_permutations) {
        cat("Permutation", i, "\n")
        # Permute host labels
        host_perm_map <- host_perm_map %>%
            mutate(permuted = sample(original))
        permuted_data <- data %>%
                        mutate(host = host_perm_map$permuted[match(host, host_perm_map$original)])
        # Recalculate splits
        splits_df_perm <- data.frame()
        
        # Bacteria
        node_list <- bac_tree$node.label
        
        for (node in node_list) {
            df <- split_balance(bac_tree, node, permuted_data)
            df_wide <- df %>% 
                mutate(node = node) %>%
                pivot_wider(names_from = split, values_from = c(clade, n_hosts, hosts))
            splits_df_perm <- rbind(splits_df_perm, df_wide)
        }
        
        # Archaea
        node_list <- ar_tree$node.label
        
        for (node in node_list) {
            df <- split_balance(ar_tree, node, permuted_data)
            df_wide <- df %>% 
                mutate(node = node) %>%
                pivot_wider(names_from = split, values_from = c(clade, n_hosts, hosts))
            splits_df_perm <- rbind(splits_df_perm, df_wide)
        }
        
        splits_df_perm <- label_splits(splits_df_perm)
        splits_summary <- splits_df_perm %>%
            group_by(node_type) %>%
            summarise(count = n()) %>%
            mutate(permutation = i)
        
        permuted_results <- rbind(permuted_results, splits_summary)
    }
    
    return(permuted_results)
}

# Run permutations for bacteria and archaea
set.seed(42)

permuted_res <- permute_splits(bac_tree, ar_tree, data, n_permutations=100)

# Calculate proportion
permuted_res <- permuted_res %>% group_by(permutation) %>%
                mutate(prop = count / sum(count))

write.csv(permuted_res, file = file.path(subdir, "mag_tree_splits_permuted.csv"), row.names=FALSE)

permuted_summary<- permuted_res %>% group_by(node_type) %>%
        summarise(median_count = median(count),
                  min_count = min(count),
                  q1_count = quantile(count, 0.25),
                  q3_count = quantile(count, 0.75),
                  max_count = max(count),
                  median_prop = median(prop),
                  min_prop = min(prop),
                  q1_prop = quantile(prop, 0.25),
                  q3_prop = quantile(prop, 0.75),
                  max_prop = max(prop))

# Combine observed and permuted results
summary_combine <-
            full_join(splits_summary, permuted_summary, by="node_type") %>% 
            # Change numeric NAs to 0
            mutate(count = ifelse(is.na(count), 0, count),
                   proportion = ifelse(is.na(proportion), 0, proportion)) %>%
            arrange(desc(count))

write.csv(summary_combine, file = file.path(subdir, "mag_tree_splits_summary.csv"), row.names=FALSE)

# Plot observed vs permuted
p1 <- ggplot(aes(x = node_type, y = count), data = permuted_res) +
        geom_boxplot(fill="white") +
        geom_point(data = summary_combine, aes(x = node_type, y = count), color="orange", size=3) +
        ylab("Number of triplets") +
        theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust=1))

p2 <- ggplot(aes(x = node_type, y = prop * 100), data = permuted_res) +
        geom_boxplot(fill="white") +
        geom_point(data = summary_combine, aes(x = node_type, y = proportion * 100), color="orange", size=3) +
        ylab("Percentage of triplets") +
        theme(axis.title.x = element_blank(), axis.text.x = element_text(hjust=1))

p <- plot_grid(p1, p2, ncol=2, align="h")

ggsave(p, file = file.path(subdir, "mag_tree_split_comparison.png"), width=6, height=4)

#######################
#### PLOT SUBTREES ####
#######################

# Function to plot subtrees for each family alongside a presence absence heatmap
plot_substree <- function(big_tree, node, mag_data, mag_presence) {
    # Get only that node
    sub_tree <- extract.clade(big_tree, node = node)
    sub_data <- mag_data %>% filter(label %in% sub_tree$tip.label)
    
    # Scale branch lengths for easier plotting
    sub_tree$edge.length <- sub_tree$edge.length / max(sub_tree$edge.length) * 0.1
    
    p_tree <- ggtree(sub_tree) %<+%
              select(sub_data, c(label, phylum, host_species)) +
              geom_tippoint(aes(colour = host_species), size=1) +
              geom_tiplab(size=2, aes(colour = host_species)) +
              scale_colour_manual(values = species_palette, name = "Assembled in")
    
    # Add presence/absence heatmap
    pres_data <- mag_presence %>% mutate(presence=TRUE) %>%
                    filter(label %in% sub_tree$tip.label) %>%
                    select(label, host_species, presence) %>%
                    pivot_wider(names_from = host_species, values_from = presence, values_fill = 0) %>%
                    column_to_rownames("label") %>% as.matrix
    
    p <- gheatmap(p_tree, pres_data, offset = 0.08, width = 1.5,
        colnames=TRUE, font.size = 2, legend_title="Presence") +
        scale_fill_manual(values = c("grey96", "orange"), labels = c("absent", "present"), name = "") +
        theme(legend.position = "none", plot.title = element_text(size = 8)) +
        ggtitle(paste0("Subtree at node ", node))
    return(p) 
}

# Plot all subtrees with interesting topologies
# For bacteria
splits_df_bac <- splits_df_filt %>% filter(node_type != "other" & domain == "Bacteria")

for (i in 1:nrow(splits_df_bac)) {
    node <- splits_df_bac$node[i]
    p <- plot_substree(bac_tree, node, data, presabs)
    ggsave(p, file = file.path(subdir, paste0("subtree_", node, ".png")), width=5, height=1)
}

# For archaea
splits_df_ar <- splits_df_filt %>% filter(node_type != "other" & domain == "Archaea")
for (i in 1:nrow(splits_df_ar)) {
    node <- splits_df_ar$node[i]
    p <- plot_substree(ar_tree, node, data, presabs)
    ggsave(p, file = file.path(subdir, paste0("subtree_", node, ".png")), width=5, height=1)
}

