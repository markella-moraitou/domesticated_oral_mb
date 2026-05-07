import argparse
import pathlib
import pandas as pd
from prep_CAT_abundance_table import parse_depth_files, parse_cat_files

# This script take the raw DRAM output (annotations.tsv),
# the named contig annotation outputs from CAT (*.contig2classification.txt)
# and the mapping stats from KMA
# and combines them to a table that shows the abundance of different genes stratified by host taxon

def summarise_annotations(raw):
    # This function takes the annotations.tsv output from DRAM and keeps only KEGG, Merops and CAZY genes
    # (in this order if one feature has been annotated with more than one database)
    # Summarise raw data
    print("Simplifing annotations table")
    raw_filtered = raw.loc[:, ["feature", "scaffold", "fasta", "ko_id", "kegg_hit", "peptidase_id", "peptidase_hit", "peptidase_identity", "cazy_ids", "cazy_hits"]]
    raw_filtered = raw_filtered[raw_filtered.feature.notna()]
    # Get features with KEGG annotations
    kegg_genes = raw_filtered[raw_filtered["ko_id"].notnull() &
                              raw_filtered["ko_id"].str.match(r"^K\d+$")].loc[:, ["feature", "scaffold", "fasta", "ko_id", "kegg_hit"]].rename(columns = {"ko_id":"gene_id", "kegg_hit":"gene_description"})
    kegg_genes["database"] = "KEGG"
    # Get features with MEROPS ANNOTATIONS (but no KEGG)
    merops_genes = raw_filtered[(raw_filtered["peptidase_id"].notnull() &
                                 raw_filtered["ko_id"].isnull()) &
                                 raw_filtered["peptidase_identity"] > 0.8].loc[:, ["feature", "scaffold", "fasta", "peptidase_id", "peptidase_hit"]].rename(columns = {"peptidase_id":"gene_id", "peptidase_hit":"gene_description"})
    merops_genes["gene_description"] = merops_genes["gene_description"].str.replace("M.* - ", "", regex=True)
    merops_genes["database"] = "MEROPS"
    # CAZY genes
    cazy_genes = raw_filtered[(raw_filtered["cazy_ids"].notnull() &
                               raw_filtered["ko_id"].isnull()) &
                               raw_filtered["peptidase_id"].isnull()].loc[:, ["feature", "scaffold", "fasta", "cazy_ids", "cazy_hits"]].rename(columns = {"cazy_ids":"gene_id", "cazy_hits":"gene_description"})
    cazy_genes["database"] = "CAZY"
    # Combine in one table
    # Create a "contig" column in annotation, that should match the contig_mags dataframe
    genes_df = pd.concat([kegg_genes, merops_genes, cazy_genes])
    genes_df.rename(columns={"scaffold": "contigName"}, inplace=True)
    genes_df = genes_df.loc[:, ["contigName", "gene_id", "gene_description", "database"]]
    return(genes_df)

def main():
    parser = argparse.ArgumentParser(description='Combine functional and taxonomic annotations of contigs with contig abundance')
    parser.add_argument('--dram_annotations', type=pathlib.Path, help='Unbinned DRAM output (raw)')
    parser.add_argument('--cat_annotations', type=pathlib.Path, help='Directory with CAT output files')
    parser.add_argument('--tax_ranks', type=str, default = "superkingdom,phylum,class,order,family,genus,species", help='Taxonomic ranks to be used for the table (options are superkingdom, phylum, class, order, family, genus, species, and can be provided as a comma-separated list e.g. "species,genus")')
    parser.add_argument('--mapping_dir', type=pathlib.Path, help='Path containing read depth files')
    parser.add_argument('--mapping_suffix', default = "_mapping_stats.txt", type=str, help='Suffix of read depth files')
    parser.add_argument('--mapping_minid', default = 98, type=int, help='Minimum percent identity to filter by (default: 98)')
    parser.add_argument('--mapping_mincov', default = 10, type=int, help='Minimum breadth of coverage to filter by (default: 10)')
    
    args = parser.parse_args()
    
    # Taxonomic ranks to use
    ranks_ls = args.tax_ranks.split(",")
    ranks_ls = [ x.strip() for x in ranks_ls ]
    
    #### Process depths
    depths = parse_depth_files(path = args.mapping_dir, suf = args.mapping_suffix, min_id = args.mapping_minid, min_breadth = args.mapping_mincov).astype(str)
    #### Process CAT output
    print("Loading CAT annotations from " + str(args.cat_annotations) + " and filtering for " + str(args.tax_ranks))   
    tax = parse_cat_files(args.cat_annotations, tax_ranks = ranks_ls).astype(str)
    # Keep only contigs present in depth table
    tax = tax[tax["contigName"].isin(depths["contigName"])].copy()
    
    # Across all column remove suffix and support value to keep unique names
    all_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    tax_filt = tax[["contigName", "lineage"] + all_ranks].copy()
    tax_filt[all_ranks] = tax_filt[all_ranks].apply(lambda col: col.str.replace(r"[a-z]__", "", regex = True).str.replace(": ?\d\.\d+$", "", regex = True))

    print("After keeping only contigs present in depth table: " + str(tax_filt.shape[0]) + " contigs\n")
    print("Unique lineages: " + str(tax_filt["lineage"].nunique()) + "\n")
    
    #### Process DRAM output
    print("Loading annotations...")
    func = pd.read_table(args.dram_annotations, dtype = {3:str, 4:str, 5:str, 6:str, 14:str, 15:str, 16:str, 22:str})
    func.rename(columns={"Unnamed: 0": "feature"}, inplace=True)
    func_simple = summarise_annotations(func)
    # Keep only contigs present in depth table
    func_simple = func_simple[func_simple["contigName"].isin(depths["contigName"])].copy()
    print("After keeping only contigs present in depth table:")
    print("Total number of annotated contigs:" + str((func_simple['contigName']).nunique()))
    print("Total number of unique genes:" + str(func_simple["gene_id"].nunique()) + "\n")
    
    #### Combine
    # Merge func with tax and depth
    print("Combining tables...\n")
    combined = pd.merge(tax_filt, func_simple, how="outer", on="contigName")
    combined = pd.merge(combined, depths, how = "outer", on = "contigName")
    combined["mapped_reads"] = pd.to_numeric(combined["mapped_reads"])
    
    # Set feature as index and save as TSV
    print("Saving combined table...\n")
    combined.to_csv("combined_contig_annotations.tsv", index = False, sep = "\t")
    
    #combined = pd.read_csv("combined_contig_annotations.tsv", sep = "\t")
    # Average depth per gene, stratified by taxon
    print("Getting gene abundance stratified by taxon...\n")
    combined_counts = combined.dropna(subset = "gene_id").groupby(
        ["Sample", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "gene_id", "gene_description", "database"
         ]).agg(mapped_reads = ("mapped_reads", "sum"))
    
    print("Saving counts...\n")
    combined_counts = combined_counts.reset_index()
    combined_counts.to_csv("gene_abundance_stratified.tsv", index = False, sep = "\t")

if __name__ == "__main__":
    main()