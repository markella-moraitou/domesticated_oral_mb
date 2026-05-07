##### GET PER SAMPLE STATS ON MEGAHIT CONTIGS #####

#### Takes as input the contig info file generated at the M1 step
#### and summarises it per sample

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(tidyr)
library(dplyr)

#### VARIABLES AND WORKING DIRECTORY ####

# Load arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  cat("Argument missing. Please provide:.\n
  Arg 1: Path of contig metadata table to be summarised\n")
} else {
  # Check arguments
  input_path <- args[1] # location input table (read count or read length) to process
  if (!file.exists(input_path)) {
    stop("Arg 1: file does not exist\n")}
}

#### Load input ####
data <- read.table(file = input_path, sep="\t", header = TRUE)

#### Process table ####
out <- data %>% group_by(sample) %>%
    # Get contig count and some summary statistics for multi (coverage) and length
    mutate(
        contig_count = n_distinct(contig),
        contig_count_lt1000 = n_distinct(contig[len > 1000]),
        multi_median = median(multi),
        len_min = min(len),
        len_median = median(len),
        len_max = max(len),
        
    ) %>%
    # Then get breakdown of flags per sample
    group_by(sample, flag, contig_count, contig_count_lt1000, multi_median, len_min, len_median, len_max) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = flag, values_from = count, names_prefix = "flag_")

#### Print output ####
output_path <- file.path(dirname(input_path), "contig_info_per_sample.txt")
write.table(out, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)