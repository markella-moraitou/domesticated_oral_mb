##### Distances between nearest domestic-wild and domestic-human MAGs analysis #####

#### Analysed MAGs to compare if a MAG from a domestic host
#### is more closely related to a MAG from a wild host or from humans

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(phytools)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "output")) # Directory with output
magdir <- normalizePath(file.path(outdir, "mags")) # Directory with MAG output
subdir <- normalizePath(file.path(magdir, "mag_distances")) # subdirectory for the output of this script

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

###########################################
##### ADD MAGS IDENTIFIED VIA MAPPING #####
###########################################

# MAGs that were found to be in more than one host species should be represented
# as radiations in the tree with 0 distance

# Get metadata
data <- rbind(bac_meta, ar_meta) %>% filter(is.tip == TRUE)

tax_data <- data %>% select(label, bin, domain, phylum, order, family, genus)

host_data <- data %>% select(label, host_species, Common.name, host_genus, Domestication)

# Unique host metadata
host_species_data <- select(host_data, c(host_species, Common.name, host_genus, Domestication)) %>%
        # Consider all sheep as domestic
        filter(Domestication != "feral") %>%
        unique

# Identify MAGs to add
mags_to_add <- presabs %>% select(label, host_species, assembly_species) %>%
                # Keep only MAGs found in more than one host species
                group_by(label) %>% filter(n() > 1) %>%
                ungroup() %>% arrange(label) %>%
                left_join(host_species_data) %>%
                # Keep cases where the presence was identified by mapping and not assemblying
                filter(!assembly_species) %>%
                # Create a new label for each additional host species
                mutate(new_label = paste(Domestication,
                                         label,
                                         str_replace_all(host_species, " ", "_"),
                                         "mapped",
                                         sep = "-")) %>%
                # Add domain info
                left_join(tax_data, by = "label")

write.csv(mags_to_add, file = file.path(subdir, "mags_to_add.csv"), row.names = FALSE)

# Give existing tip names to match
rename_tips <- data.frame(label = c(bac_tree$tip.label, ar_tree$tip.label)) %>%
               mutate(host_species = host_data$host_species[match(label, host_data$label)]) %>%
               mutate(Domestication = host_species_data$Domestication[match(host_species, host_species_data$host_species)]) %>%
               mutate(new_label = paste(Domestication,
                                        label,
                                        str_replace_all(host_species, " ", "_"),
                                        "assembled",
                                        sep = "-"))

write.csv(mags_to_add, file = file.path(subdir, "rename_tips.csv"), row.names = FALSE)

#### Modify bacteria tree ####
bac_tree_edit <- bac_tree

bac_tree_edit$tip.label <- rename_tips$new_label[match(bac_tree$tip.label, rename_tips$label)]

bacs_to_add <- mags_to_add %>% filter(domain == "Bacteria")

# Add new tips with 0 distance
for(i in 1:nrow(bacs_to_add)) {
    mag <- rename_tips$new_label[rename_tips$label == bacs_to_add$label[i]]
    # Exit with error if MAG has length more than one
    if (length(mag) > 1) {
        stop(paste0("Error: MAG ", bacs_to_add$label[i], " has multiple tip labels in the tree"))
    }
    new_mag <- bacs_to_add$new_label[i]
    # Find the node number for the original MAG
    node_num <- which(bac_tree_edit$tip.label == mag)
    if (length(node_num) > 1) {
        warning(paste0("Warning: MAG ", mag, " has multiple tip labels in the tree, skipping addition of new tip"))
        next
    }
    # Add new tip with 0 distance
    bac_tree_edit <- bind.tip(bac_tree_edit, new_mag, where = node_num, edge.length = 0)
}

#### Modify archaea tree ####
ar_tree_edit <- ar_tree

ar_tree_edit$tip.label <- rename_tips$new_label[match(ar_tree$tip.label, rename_tips$label)]

ars_to_add <- mags_to_add %>% filter(domain == "Archaea")

# Add new tips with 0 distance
for(i in 1:nrow(ars_to_add)) {
    mag <- rename_tips$new_label[rename_tips$label == ars_to_add$label[i]]
    # Exit with error if MAG has length more than one
    if (length(mag) > 1) {
        stop(paste0("Error: MAG ", ars_to_add$label[i], " has multiple tip labels in the tree"))
    }
    new_mag <- ars_to_add$new_label[i]
    # Find the node number for the original MAG
    node_num <- which(ar_tree_edit$tip.label == mag)
    if (length(node_num) > 1) {
        warning(paste0("Warning: MAG ", mag, " has multiple tip labels in the tree, skipping addition of new tip"))
        next
    }
    # Add new tip with 0 distance
    ar_tree_edit <- bind.tip(ar_tree_edit, new_mag, where = node_num, edge.length = 0)
}

#### Calculated phylogenetic distances ####
bac_dist <- cophenetic.phylo(bac_tree_edit) %>% data.frame() %>% rownames_to_column(var = "MAG1") %>%
    pivot_longer(-MAG1, names_to = "MAG2", values_to = "Distance") %>%
    # Fix names
    mutate(MAG2 = gsub(MAG2, pattern = ".", replacement = "-", fixed = TRUE)) %>%
    # Order MAG names so that domestic are preferably in MAG1 %>%
    mutate(
        MAG1_orig = MAG1,
        MAG2_orig = MAG2,
        MAG1 = pmin(MAG1_orig, MAG2_orig),
        MAG2 = pmax(MAG1_orig, MAG2_orig)
    ) %>%
    ungroup() %>%
    select(-MAG1_orig, -MAG2_orig) %>%
    # Separate columns
    separate(MAG1, into = c("domestication_1", "MAG1", "host_species_1", "method_1"), sep = "-", extra = "warn") %>%
    separate(MAG2, into = c("domestication_2", "MAG2", "host_species_2", "method_2"), sep = "-", extra = "warn")

ar_dist <- cophenetic.phylo(ar_tree_edit)  %>% data.frame() %>% rownames_to_column(var = "MAG1") %>%
    pivot_longer(-MAG1, names_to = "MAG2", values_to = "Distance") %>%
    # Fix names
    mutate(MAG2 = gsub(MAG2, pattern = ".", replacement = "-", fixed = TRUE)) %>%
    # Order MAG names so that domestic are preferably in MAG1 %>%
    mutate(
        MAG1_orig = MAG1,
        MAG2_orig = MAG2,
        MAG1 = pmin(MAG1_orig, MAG2_orig),
        MAG2 = pmax(MAG1_orig, MAG2_orig)
    ) %>%
    ungroup() %>%
    select(-MAG1_orig, -MAG2_orig) %>%
    # Separate columns
    separate(MAG1, into = c("domestication_1", "MAG1", "host_species_1", "method_1"), sep = "-", extra = "warn") %>%
    separate(MAG2, into = c("domestication_2", "MAG2", "host_species_2", "method_2"), sep = "-", extra = "warn")

#### Combine ####

phy_dist <- rbind(bac_dist, ar_dist) %>%
    mutate(host_species_1 = str_replace(host_species_1, "_", " "),
           host_species_2 = str_replace(host_species_2, "_", " "))

write.csv(phy_dist, file = file.path(subdir, "mag_distances.csv"), row.names = FALSE)

# Keep only distances were one MAG is a domestic and the other is not
phy_dist_filt <- phy_dist %>%
                 filter((domestication_1 == "domestic" & domestication_2 != "domestic") |
                        (domestication_2 == "domestic" & domestication_1 != "domestic"))

##########################################################
#### CALCULATE DISTANCES TO NEAREST WILD OR HUMAN MAG ####
##########################################################

# Write function that for each domestic MAG, get the closest MAG from a wild host and from humans
get_distances <- function(mag, species, distance_table) {
    dom_wild <- list("Equus caballus" = "Equus quagga",
                     "Ovis aries" = "Ovis ammon",
                     "Sus domesticus" = "Sus scrofa")
    # Filter distances for the given MAG
    dist_filt <- distance_table %>%
        # Keep only distances involving the given MAG and host species (all domestic should be in columns MAG1 and host_species_1)
        filter((MAG1 == mag & host_species_1 == species)) %>%
        # Are we comparing with humans or wild counterparts?
        mutate(comparing_with = domestication_2) %>%
        # Was the MAG idenitified via assembly or mapping?
        mutate(presence = case_when(method_1 == "mapped" | method_2 == "mapped" ~ "mapped",
                                    TRUE ~ "assembled"))
    
    # Get the smallest distance for each comparison
    dist_min <- dist_filt %>%
        group_by(comparing_with) %>%
        slice_min(Distance) %>%
        select(comparing_with, Distance, presence) %>%
        # Add relevant data
        mutate(domesticate_MAG = mag,
               domesticate_host = species) %>% unique() %>%
        ungroup() %>% mutate(presence = case_when(any(presence == "mapped") ~ "mapped",
                                        TRUE ~ "assembled"))
    return(dist_min)    
}

## Do this for all MAGs found in domestic hosts
dom_mags <- phy_dist_filt %>%
            select(MAG1, host_species_1) %>%
            rename(MAG = MAG1,
                   host_species = host_species_1) %>%
            mutate(host_genus = str_remove(host_species, " .*")) %>% unique

dw_dh_distances <- data.frame()

for(i in 1:nrow(dom_mags)) {
    mag = dom_mags$MAG[i]
    species = dom_mags$host_species[i]
    dist_res <- get_distances(mag, species, phy_dist_filt)
    dw_dh_distances <- rbind(dw_dh_distances, dist_res)
    if(i %% 10 == 0) {
        cat(paste0("Processed ", i, " / ", nrow(dom_mags), " domestic MAGs\n"))
    }
}

dw_dh_distances <- unique(dw_dh_distances)

write.csv(dw_dh_distances, file = file.path(subdir, "domestic_mag_dw_dh_distances.csv"), row.names = FALSE)

# Make wider
dw_dh_wide <- dw_dh_distances %>%
    pivot_wider(names_from = comparing_with, values_from = Distance, id_cols = c(domesticate_MAG, domesticate_host, presence))

# Plot
p <- ggplot(data = dw_dh_wide, aes(x = wild, y = human, colour = domesticate_host, shape = presence)) +
     geom_jitter(alpha = 0.5, size = 2, height = 0.01, width = 0.01) +
     scale_colour_manual(values = species_palette, name = "Domesticate host") +
     scale_shape_manual(values = c("mapped" = 4, "assembled" = 16), name = "MAG presence via",
                        labels = c("assembly", "mapping")) +
     facet_grid(domesticate_host ~ .) +
     labs(x = "Distance to nearest wild counterpart MAG",
          y = "Distance to nearest human MAG") +
     geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
     theme(legend.direction = "vertical", legend.position = "right")

# Label mags that are closer to human MAGs
p <- p + geom_text(data = subset(dw_dh_wide, human < wild),
                   aes(label = domesticate_MAG),
                   vjust = 1, hjust = 0.5, size = 1.7, color = "black") +
    xlim(c(0,1.8))

ggsave(p, filename = file.path(subdir, "dw_dh_distances.png"), width = 6, height = 5)
