
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import matplotlib.gridspec as gridspec

# This script plots a heatmap based on the ANI values

#### Heatmap Data ####
# Load data
def fastani_heatmap(ani_df, dereplicated_bins = None):
    print("Preparing FastANI heatmap...")
    if dereplicated_bins is None:
        dereplicated_bins = []  
    # Simplify
    ani_df = ani_df.iloc[:,0:3]
    ani_df = ani_df.replace(to_replace=r'MEGAHIT-.*-', value='', regex=True)
    ani_df = ani_df.replace(to_replace=r'\.fa.*', value='', regex=True).drop_duplicates()
    ani_df.columns = ["bin1", "bin2", "ANI"]
    # Get wide format
    pivot_ani_df = ani_df.pivot(index='bin1', columns='bin2', values='ANI')
    pivot_ani_df = pivot_ani_df.fillna(0)
    # Get order of rows/columns
    linkage = sch.linkage(pivot_ani_df, method='centroid')
    dendro = sch.dendrogram(linkage, no_plot=True)
    order = list(pivot_ani_df.iloc[1,dendro['leaves']].index.values)
    # Reorder
    pivot_ani_df = pivot_ani_df.loc[order, order]
    # Get back NA values
    pivot_ani_df.replace(0, np.nan, inplace=True)
    # Create a mask for highlighting specific bins
    highlight_mask = pd.DataFrame(0, index=pivot_ani_df.index, columns=pivot_ani_df.columns)
    for bin_name in dereplicated_bins:
        if bin_name in pivot_ani_df.index and bin_name in pivot_ani_df.columns:
            highlight_mask.loc[bin_name, bin_name] = 1
    print("Highlight mask:")
    print(highlight_mask)
    # Plot heatmap
    def plot_heatmap(ax):
        sns.heatmap(pivot_ani_df,
                    cmap = 'viridis', ax = ax,
                    cbar_ax = None, cbar = False,
                    annot=highlight_mask.replace(0, np.nan) # Overlay dots where mask is 1
                    )
        ax.set_title('FastANI Heatmap')
    return plot_heatmap, order

#### PyDamage ####
def pydamage_scatter(dam_df, order):
    print("Preparing PyDamage scatterplot...")
    dam_df = dam_df[dam_df['damage_model_p'].notna()]
    # Simplify bin names
    dam_df = dam_df.replace(to_replace=r'MEGAHIT-.*-', value='', regex=True)
    dam_df = dam_df.replace(to_replace=r'\.fa.*', value='', regex=True).drop_duplicates()
    # Make bin categorical and order according to the FastANI heatmap
    dam_df.bin = pd.Categorical(dam_df.bin, categories = order, ordered = True)
    dam_df = dam_df.sort_values('bin')
    # Log-transform p-value with a pseudocount of 10^-4
    dam_df["neg_log10_p"] = -np.log10(dam_df.damage_model_p + 10**(-4))
    # Plot
    def plot_scatter(ax):
        sns.boxplot(x='bin', y='neg_log10_p', data=dam_df, ax=ax, showfliers=False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title('pyDamage Scatter Plot')
    return plot_scatter

# Extract taxonomic assignments at a specific rank from CAT output
def extract_tax(value, rank):
    if ';' + rank + '__' in value:
        return value.split(';' + rank + '__')[1].split(';')[0]
    else:
        return np.nan

#### Taxonomic assignment ####
def taxonomy(bins_df, order):
    print("Preparing taxonomy plot...")
    # Extract order info form classification column
    bins_df.loc[:, 'phylum'] = bins_df['classification'].apply(lambda x: extract_tax(x, "p"))
    # Simplify bin names
    bins_df = bins_df.replace(to_replace=r'MEGAHIT-.*-', value='', regex=True)
    bins_df = bins_df.replace(to_replace=r'\.fa.*', value='', regex=True).drop_duplicates()
    # Reorder according to heatmap order
    bins_df.bin = pd.Categorical(bins_df.bin, categories = order, ordered = True)
    bins_df = bins_df.sort_values('bin')
    # Get colour for plotting
    unique_values = list(bins_df['phylum'].unique())
    palette = sns.color_palette("husl", len(unique_values))
    color_dict = dict(zip(unique_values, palette))
    bins_df['color'] = bins_df['phylum'].map(color_dict)
    # Get a dataframe for plotting
    plot_df = pd.DataFrame({'phylum' : bins_df['color'].apply(lambda x: palette.index(x))})
    plot_df.index = list(bins_df['bin'].values)
    def plot_taxonomic_bar(ax):
        sns.heatmap(plot_df.transpose(), cmap=palette, ax = ax, cbar = False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title('Bin phylum classification')
    # Plot a legend for the taxonomic bar
    def plot_legend(ax):
        for i, value in enumerate(unique_values):
            ax.bar(i, 1, color = color_dict[value], label = value)
        ax.set_xticks(range(len(unique_values)))
        # Place axus x labels inside the bars and rotate them vertically
        ax.set_xticklabels(unique_values)
        ax.set_yticks([])
        ax.set_title('Phylum legend')
    return plot_taxonomic_bar, plot_legend

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Summary plot of the assembled MAGs, showing ANI, taxonomic classification and pydamage results')
    parser.add_argument('--fastani', type=str, help='Path to FastANI output TSV')
    parser.add_argument('--pydamage', type=str, help='Path to pyDamage output; Must contain "contig" and "bin" column, as well as the pydamage output columns')
    parser.add_argument('--bin_metadata', type=str, help='Path to bin metadata, containing taxonomic info')
    parser.add_argument('--dereplicated_bins', required = False, type=str, help='List of dereplicated bins (optional)')
    parser.add_argument('--output', type=str, nargs='+', help='Name for the output plot (PNG)')
    args = parser.parse_args()
    output = args.output[0]
    print("Loading input...")
    fastani = pd.read_table(args.fastani, sep = " ", header = None)
    print(fastani.head())
    pydamage = pd.read_csv(args.pydamage, header = 0, index_col = 0)
    # Only keep pydamage output from bins in fastani table
    pydamage_filt = pydamage[pydamage['bin'].isin(fastani[0])]
    print(pydamage_filt.head())
    bin_meta = pd.read_csv(args.bin_metadata, header = 0).loc[:,["bin", "classification"]]
    bin_meta_filt = bin_meta[bin_meta['bin'].isin(fastani[0])]
    print(bin_meta_filt.head())
    # If dereplicated bins are provided, provide this as argument to fastani_heatmap function
    if args.dereplicated_bins:
        derep_bins = pd.read_csv(args.dereplicated_bins, header = None).squeeze().tolist()
        # Get FastANI heatmap, and also extract bin order
        result = fastani_heatmap(ani_df = fastani, dereplicated_bins = derep_bins)
    heatmap = result[0]
    order = result[1]
    # Create scatterplot from python data and plot the same way as in the heatmap
    scatter=pydamage_scatter(pydamage_filt, order)
    # Get taxonomic classification sideplot
    tax_bar = taxonomy(bin_meta_filt, order)
    tax_plot = tax_bar[0]
    tax_legend = tax_bar[1]
    # Save combined plot
    gs = gridspec.GridSpec(4, 1, height_ratios=[1, 5, 0.5, 0.5])
    fig = plt.figure(figsize=(30, 40))
    # Plot 1
    ax1 = fig.add_subplot(gs[0])
    scatter(ax1)
    ax1.set_xticklabels([], visible = False)
    # Plot 2
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    heatmap(ax2)
    ax2.set_xticklabels([], visible = False)
    # Plot 3
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    tax_plot(ax3)
    # Legend for plot 3
    ax4 = fig.add_subplot(gs[3])
    tax_legend(ax4)
    plt.savefig(output, dpi=300, bbox_inches='tight')  # Save as PNG with high resolution

if __name__ == "__main__":
    main()
