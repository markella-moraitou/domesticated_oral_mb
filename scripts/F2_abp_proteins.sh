#!/bin/bash -l

# This script looks for abp proteins in the contig ORFs

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR} # Output subdirectory for the output of this specific script
genes=${outdir}F1_contig_annotations/raw/working_dir/clustered_contigs_rep_seq/genes.annotated.faa # Input faa sequences
contigmapping=${outdir}/T2_contig_mapping
contigtax=${outdir}/T1_contig_taxonomy
subdir=${outdir}/F2_abp_proteins/ # Output subdirectory

mkdir -p $subdir # Create subdirectory if not there

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${contigdir}"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

cd $subdir

#####################
####  RUN HMMER  ####
#####################

abpa_fasta=${indir}/abpA.fasta

# Align
mafft --auto $abpa_fasta > abpA_alignment.fa

# Build HMM profile
hmmbuild abpA_profile.txt abpA_alignment.fa

# Search
hmmsearch -o abpA_hmmer_results.out -A abpA_hmmer_alignment.txt --tblout abpA_hmmer_table.txt -E 1 abpA_profile.txt $genes 

################################
#### GET CONTIG ABUNDANCES  ####
################################

# Get contig names where abpA was
grep -v "#"  abpA_hmmer_table.txt | awk '{print $1}' | sed 's/clustered_contigs_rep_seq_//' | sed 's/_[0-9]*$//' > abpA_contigs.tmp

# Iterate through mapping stats from all samples to the combined contig database
ls $contigmapping/*_mapping_stats.txt| while read file
do  
    # If first file, get header
    if [ "$file" == "$(ls $contigmapping/*_mapping_stats.txt | head -1)" ]; then
        echo -e "sample\t$(head -1 $file)" > abpA_mapping_stats.txt
    fi
    sample=$(basename $file _mapping_stats.txt)
    # Subset to the contigs containing abp and add the mapped sample as a new column
    grep -w -f abpA_contigs.tmp $file | awk -v sample=$sample '{print sample "\t" $0}' >> abpA_mapping_stats.txt
done

# Also get taxonomy of these contigs

head -1 ${contigtax}/all.contig2classification_named.txt > abpA_tax.txt
grep -f abpA_contigs.tmp ${contigtax}/all.contig2classification_named.txt | grep -v "no taxid assigned" >> abpA_tax.txt

rm abpA_contigs.tmp
