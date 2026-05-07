import argparse
import pathlib
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client
import logging

## Use .res file to calculate abundance of bins from KMA mapping

logging.getLogger("distributed.shuffle._scheduler_plugin").setLevel(logging.ERROR)

def calculate_mapped_reads(frag_file, res_file):
    frag_file = frag_file.set_index('template')
    # Count occurrences of each index value
    out = frag_file.index.value_counts()
    out = out.reset_index()
    out.columns = ['template', 'mapped_reads']
    out['mapped_reads'] = out['mapped_reads'].astype('int64')
    # Repartition for faster processing
    out = out.repartition(npartitions=200).simplify()
    res_file = res_file.repartition(npartitions=200).simplify()
    # Add identity and coverage information
    merged = dd.merge(out, res_file, on='template', how='outer')
    # Separate contig name and length
    merged['template'] = merged['template'].fillna("").astype(str)
    merged['contigName'] = merged['template'].str.split(' ', n=1).str[0]
    merged['contigLength'] = merged['template'].str.split("len=", n=2).str[1].astype('int64')
    merged = merged[['contigName', 'contigLength', 'mapped_reads', 'identity', 'coverage']]
    return merged

def mapping_stats_per_bin(reads_per_contig, contig_bin_map):
    # Summarise mapping statistics per bin
    # Mapped reads: sum of mapped reads for all contigs in bin
    # Coverage: weighted average based on contig length
    # Identity: weighted average based on contig length
    # Ensure correct dtypes
    reads_per_contig['contigName'] = reads_per_contig['contigName'].astype(str)
    contig_bin_map['contigName'] = contig_bin_map['contigName'].astype(str)
    # Merge data frames
    merged = dd.merge(reads_per_contig, contig_bin_map, on=['contigName', 'contigLength'], how='outer')
    # Repartition and simplify to make computations faster
    merged = merged.repartition(npartitions=200)
    merged = merged.simplify().persist()
    print("Unique bins in merged data: " + str(merged['bin'].nunique().compute()))
    print("Unique contigs in merged data: " + str(merged['contigName'].nunique().compute()))
    print("Contigs with no mapped reads: " + str(merged['mapped_reads'].isna().sum().compute()))
    # Fill NaN values for contigs with no mapped reads
    merged['mapped_reads'] = merged['mapped_reads'].fillna(0)
    merged['identity'] = merged['identity'].fillna(0)
    merged['coverage'] = merged['coverage'].fillna(0)
    # Prep calculations for weighted average before turning to pandas
    merged['identityxlength'] = merged['identity'] * merged['contigLength']
    merged['coveragexlength'] = merged['coverage'] * merged['contigLength']
    # Group by bin and sum mapped reads
    bin_abundance = merged.groupby('bin').agg({'mapped_reads': 'sum',
                                               'identityxlength': 'sum',
                                               'coveragexlength': 'sum',
                                               'contigLength': 'sum'}).reset_index()
    bin_abundance = bin_abundance.compute()
    # Compute weighted averages
    bin_abundance['identity'] = bin_abundance['identityxlength'] / bin_abundance['contigLength']
    bin_abundance['coverage'] = bin_abundance['coveragexlength'] / bin_abundance['contigLength']
    bin_abundance = bin_abundance[['bin', 'mapped_reads', 'identity', 'coverage', 'contigLength']]
    bin_abundance = bin_abundance.rename(columns={'contigLength': 'total_Length'})
    print("Summary stats")
    print(bin_abundance.describe())
    # Remove bins with identity less than 50% and coverage less than 10%
    bin_abundance = bin_abundance[bin_abundance['identity'] >= 50]
    bin_abundance = bin_abundance[bin_abundance['coverage'] >= 10]
    return bin_abundance

def main():
    argparser = argparse.ArgumentParser(description='Calculate contig abundance from KMA .res file')
    argparser.add_argument('--frag_file', type=pathlib.Path, help='.frag file from KMA mapping')
    argparser.add_argument('--res_file', type=pathlib.Path, help='.res file from KMA mapping')
    argparser.add_argument('--out_prefix', type=str, help='Output prefix for abundance table')
    argparser.add_argument('--processes', type=int, help='Number of processes to use', default=1)
    argparser.add_argument('--contig_bin_map', required = False, type=pathlib.Path, help='CSV file mapping contig names to bin names')
    args = argparser.parse_args()
    frag = args.frag_file
    res = args.res_file
    proc = args.processes
    map = args.contig_bin_map
    # Set up multi-core processing
    client = Client(n_workers=proc, threads_per_worker=1)  # Adjust as needed
    # Load files
    print("Loading KMA .frag file: " + str(frag) + "...")
    df_frag = dd.read_csv(frag, sep = "\t", header=None, blocksize="100MB")
    df_frag = df_frag.loc[:,0:5] # Keep only relevant columns
    df_frag.columns = ["read", "n_templates", "score", "start", "end", "template"]
    # Get correct filetypes
    df_frag['template'] = df_frag['template'].astype(str)
    print("Loading KMA .res file: " + str(res) + "...")
    df_res = dd.read_csv(res, sep = "\t", blocksize="100MB")
    df_res = df_res[['#Template', 'Query_Identity', 'Template_Coverage']]
    df_res.columns = ['template', 'identity', 'coverage']
    df_res['template'] = df_res['template'].astype(str)
    # Calculate mapped reads per contig
    print("Calculating number of mapped reads per template...")
    mapped_reads = calculate_mapped_reads(df_frag, df_res)
    if map:
        print("Loading contig-bin map: " + str(map) + "...")
        df_map = dd.read_csv(map, blocksize="100MB")[["contig", "bin", "length"]].rename(columns={"contig": "contigName", "length": "contigLength"})
        # Turn length column to int
        df_map['contigLength'] = df_map['contigLength'].astype(int)
        # Keep only contigs already in df_frag
        contig_names = mapped_reads['contigName'].compute()
        df_map = df_map[df_map['contigName'].isin(contig_names)]
        print("Calculating mapping statistics per bin...")
        # Calculate mapped reads per bin
        mapped_reads_per_bin = mapping_stats_per_bin(mapped_reads, df_map)
        mapped_reads = mapped_reads_per_bin
    else:
        mapped_reads = mapped_reads.compute()
    # Save table
    print("Saving mapping statistics to: " + args.out_prefix + "_mapping_stats.txt")
    mapped_reads.to_csv(args.out_prefix + "_mapping_stats.txt", sep = "\t", index = False)

if __name__ == "__main__":
    main()    
