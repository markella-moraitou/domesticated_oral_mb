#!/bin/bash -l

#SBATCH -p shared
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH --mem=60GB
#SBATCH -J derep
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script using dRep to dereplicate MAGs from the same species

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

subdir=${outdir}M4_dereplicated_bins/ # Output directory for this script
magdir=${outdir}M3_MAG_metadata/MQ_HQ_BINS/ # Input directory with 'good' MAGs
maginf=${outdir}M3_MAG_metadata/mq_hq_bin_metadata.csv # Input file with metadata for the MAGs

mkdir -p $subdir # Create subdirectory if not there

## Load modules
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${magdir}"
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

echo "Number of tasks: ${processes}"

##################
#### PROCESS  ####
##################

cd $subdir

echo "Dereplicating $magdir"
dRep dereplicate $subdir -p $processes -g ${magdir}/*fa.gz --ignoreGenomeQuality --genomeInfo $maginf -sa 0.98

ls dereplicated_genomes/*fa.gz > dereplicated_bins_list.txt