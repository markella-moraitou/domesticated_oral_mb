## Contents

**scripts directory**

├──```modules``` **directory**: scripts used as modules in main scripts

*Scripts to be run in preparation for the pipeline:*
         
├──```0_demultiplex_inhouse.sh```: Demultiplexing according to barcode combination, for libraries from this study.

├──```0_download_ncbi_genomes.sh```: Downloading reference genomes from NCBI.

├──```0_download_prev_published_data.sh```: Download previously published metagenomes using the ```previously_published_data.csv``` file.

├──```0_get_ameta_dbs.sh```: Download data necessary for running the aMeta pipeline (https://github-com.translate.goog/NBISweden/aMeta).

*Scripts for the preprocessing of metagenomic reads and assembly-free taxonomic classification*
         
├──```1_seq_cleanup.sh```: Adapter, barcode and low quality trimming, merging, and mean quality and length filtering using fastp (https://github.com/OpenGene/fastp).

├──```2_deduplication.sh```: Remove PCR duplicates using BBTools dedupe module (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/dedupe-guide/).

├──```3a_prepare_mapping_refs.sh```: Prepare references for mapping metagenome against the human, host and PhiX genomes.

├──```3b_mapping.sh```: Map against references using bwa aln (https://github.com/lh3/bwa) and separate the unmapped reads for downstream analysis.

├──```4_tax_classification.sh```: Taxonomic classification using KrakenUniq (https://github.com/fbreitwieser/krakenuniq).

├──```4_tax_classification_ameta.sh```: Taxonomic classificationg using aMeta.

*Scripts for assembly and analysis of MAGs (Metagenome-Assembled Genomes)*

├──```M1a_MAG_assemblies.sh```: Assembling metagenomic contigs using MEGAHIT (https://github.com/voutcn/megahit).

├──```M1b_contig_dereplication.sh```: Clustering contigs using MMseqs2 (https://github.com/soedinglab/MMseqs2).

├──```M3a_select_bins.sh```: Select bins for downstream analysis and collect metadata.

├──```M3b_checkM.sh```: Run checkM to get completeness and contamination per MAG (https://github.com/Ecogenomics/CheckM).

├──```M3c_bin_taxonomy.sh```: Run GTDBtk to get taxonomic assignments (https://github.com/Ecogenomics/GTDBTk).

├──```M4a_dereplicate_bins.sh```: Dereplicate MAGs using dRep (https://github.com/MrOlm/drep).

├──```M4b_ANI.sh```: Calculate and plot average nucleotide differences using FastANI (https://github.com/ParBLiSS/FastANI).

├──```M4c_bin_abundances.sh```: Maps metagenomics reads to a combined MAG references using KMA and get abundances (https://github.com/genomicepidemiology/kma).

*Scripts for the functional profiling using contigs*

├──```F1_contig_annotations.sh```: Annotate contigs using DRAM (https://github.com/WrightonLabCSU/DRAM).

├──```F2_extract_MAG_annotations.sh```: Extract the functional annotations for each bin from the unbinned DRAM output.

*Scripts for the taxonomic classification using contigs*

├──```T1_contig_taxonomy.sh```: Assign taxonomy to contigs using CAT (https://github.com/AnyiHu/CAT).

├──```T2_contig_mapping.sh```: Map metagenomic reads to cluster contig reference and calculate abundances (https://github.com/genomicepidemiology/kma).

*Get final tables for R analysis*

└──```Final_tables.sh```: Get tables for analysis of contigs combining CAT taxonomies, depths (from nfcore/mag) and DRAM annotations
