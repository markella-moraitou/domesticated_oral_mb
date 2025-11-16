#!/bin/bash -l

#SBATCH -p main
#SBATCH -n 100
#SBATCH --mem 200G
#SBATCH -t 10:00:00
#SBATCH -J megahit_derep
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script uses MEGAHIT to perfrom de novo assembly from the metagenomic reads
# Also concatenates and dereplicates to use a single reference for taxonomic and functional profiling and mapping

## Change to the directory containing the script
[ ! -z $SLURM_SUBMIT_DIR ] && cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR} # Output subdirectory for the output of this specific script

subdir=${outdir}M1_contig_assemblies/ # Output subdirectory for MAG assembly

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Load modules and activate environment
conda activate oral_mb_evol
module load bioinfo-tools megahit R/4.2

## Print some info
echo "Taking input from: ${outdir}/2_mapped_seqs/"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z $SLURM_NTASKS ]; then
  processes=1
else
  processes=$SLURM_NTASKS
fi
echo "Number of processes: ${processes}"

cd $subdir

#####################################
#### CONCATENATE AND DEREPLICATE ####
#####################################

# Check if file exists and delete if so
[ -f all_contigs_1000bp.fa ] && rm all_contigs_1000bp.fa
[ -f all_contigs_1000bp.fa.gz ] && rm all_contigs_1000bp.fa.gz

# Concatenate all contigs into a single file adding sample name to 
ls *_final_contigs_1000bp.fa | while read contigs
do
  sample=${contigs%_final_contigs_1000bp.fa}
  cat $contigs | sed "s/>/>${sample}_/g" >> all_contigs_1000bp.fa
done

pigz -p $processes all_contigs_1000bp.fa

# Dereplicate using MMseqs

mmseqs easy-cluster all_contigs_1000bp.fa.gz clustered_contigs tmp --min-seq-id 0.98 -c 0.8 --cov-mode 5 --threads 2 --split-memory-limit 90G
