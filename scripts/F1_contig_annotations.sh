#!/bin/bash -l

#SBATCH -p long
#SBATCH -n 40
#SBATCH -t 5-00:00:00
#SBATCH -J func_contigs
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err
#SBATCH --profile=task

# This script annotates contigs using DRAM

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
contigs=${outdir}M1_contig_assemblies/clustered_contigs_rep_seq.fasta # Input contigs
subdir=${outdir}/F1_contig_annotations/ # Output subdirectory for DRAM output
dram_refdir=${REFDIR}/DRAM_data # DRAM database

mkdir -p $subdir # Create subdirectory if not there

## Load modules and activate environment
conda activate DRAM

## Print some info
echo "Taking input from: ${contigdir}"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z $SLURM_NTASKS ]
then
    processes=1
else
    processes=$SLURM_NTASKS
fi

echo "Number of processes: ${processes}"

cd $subdir

####################
####  RUN DRAM  ####
####################

# Set up databases
DRAM-setup.py import_config --config_loc config.txt

# Annotate
DRAM.py annotate -i $contigs -o raw --threads 20 --min_contig_size 1000

# Distill annotations
DRAM.py distill -i raw/annotations.tsv -o distillation