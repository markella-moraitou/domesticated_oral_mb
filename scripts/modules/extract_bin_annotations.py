import sys
import argparse
import pandas as pd
#import numpy as np
#import os

# This script take the raw DRAM output and a file mapping contigs to MAGs with columns "contig" and "bin"
# and outputs a table similar to the DRAM output that groups annotations by bin, not unbinned fasta file

def main():
    parser = argparse.ArgumentParser(description='Extract MAG annotations from unbinned DRAM output, using a file mapping contigs to MAGs')
    parser.add_argument('--contig_mag_mapping', type=str, help='File mapping contigs to MAGs with columns "contig" and "bin" (CSV)')
    parser.add_argument('--dram_annotations', type=str, help='Unbinned DRAM output (raw) (TSV)')
    args = parser.parse_args()
    
    # Load tables
    print("Loading mapping file...")
    contig_mags = pd.read_csv(args.contig_mag_mapping, index_col=False)
    
    # Select only bin and contig columns
    contig_mags = contig_mags[["contig", "bin"]]
    print("\n")
    print("Number of bins:" + str(contig_mags['bin'].nunique()))
    print("Number of contigs in bins:" + str(contig_mags['contig'].nunique()))
    print("\n")
    
    print("Loading annotations...")
    annotations = pd.read_table(args.dram_annotations, dtype = {3:str, 4:str, 5:str, 6:str, 14:str, 15:str, 16:str, 22:str})
    annotations.rename(columns={"Unnamed: 0": ""}, inplace=True)
    print("Total number of annotated contigs:" + str((annotations['fasta'] + annotations['scaffold']).nunique()))
    print("Total number of features:" + str(annotations[""].nunique()))
    print("\n")
    
    # Create a "contig" column in annotation, that should match the contig_mags dataframe
    annotations['contig'] = annotations['fasta'].str.replace(r'_final_contigs', '') + "_" + annotations['scaffold']
    
    # Get annotation for contigs that ended up binned in MAGs
    print("Creating new tables with features per bin...")
    mag_annotations = contig_mags.merge(annotations, on='contig', how='inner')
    print("Number of bins with annotations:" + str(mag_annotations["bin"].nunique()))
    print("Number of contigs with annotations:" + str(mag_annotations["contig"].nunique()))
    print("Features in bins:" + str(mag_annotations[""].nunique()))
    
    # Set bin as the new fasta column and drop bin and contig
    mag_annotations["fasta"] = mag_annotations["bin"]
    mag_annotations = mag_annotations.drop(["bin", "contig"], axis=1)
    
    # Set feature as index and save as TSV
    mag_annotations = mag_annotations.set_index('')
    mag_annotations.to_csv("MAG_annotations.tsv", index = True, sep = "\t") 

if __name__ == "__main__":
    main()