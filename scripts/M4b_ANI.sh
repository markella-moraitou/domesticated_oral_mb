#!/bin/bash -l

#SBATCH -p shared
#SBATCH -n 15
#SBATCH -t 12:00:00
#SBATCH --mem=60GB
#SBATCH -J fastani
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script runs FastANI

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

subdir=${outdir}M4_dereplicated_bins/ # Output subdirectory for MAG assembly
magdir=${outdir}M3_MAG_metadata/MQ_HQ_BINS/ # Input directory with 'good' MAGs
maginf=${outdir}M3_MAG_metadata/mq_hq_bin_metadata.csv # Input file with metadata for the MAGs

## Determine samples to be run
sample_list=${indir}/file_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z "$SLURM_NTASKS" ]; then
  processes=1
else
  processes=$SLURM_NTASKS
fi

echo "Number of tasks: ${SLURM_NTASKS}"

cd $subdir

##################
#### FASTANI #####
##################

# Get list of HQ genomes
ls ${magdir}/* > mq_hq_genome_list.tmp

# Run FastANI using the same list as reference and query
fastANI --rl mq_hq_genome_list.tmp --ql mq_hq_genome_list.tmp --fragLen 1000 -t $processes -o fastANI_results.tsv

# Simplify file names
awk '{gsub(".*//", "", $1); gsub(".*//", "", $2); print}' fastANI_results.tsv > temp
mv temp fastANI_results.tsv

# Plot heatmap
python $scriptdir/modules/ANI_heatmap.py --fastani fastANI_results.tsv \
      --bin_metadata ../M3_MAG_metadata/mq_hq_bin_metadata.csv --pydamage ../M3_MAG_metadata/contig_metadata.csv --dereplicated_bins dereplicated_bins_list.txt \
      --output ANI_heatmap.png 
