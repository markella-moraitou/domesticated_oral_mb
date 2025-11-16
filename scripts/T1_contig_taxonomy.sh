#!/bin/bash -l

#SBATCH -p main
#SBATCH -n 256
#SBATCH -t 1-00:00:00
#SBATCH -J cat
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Assign taxonomy to contigs using CAT

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

subdir=${outdir}T1_contig_taxonomy/ # Output subdirectory
contigs=${outdir}M1_contig_assemblies/clustered_contigs_rep_seq.fasta # Input contigs
catdb=$REFDIR/CAT_DB # Where the CAT db will be created

mkdir -p $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${contigdir}"
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
echo "Number of tasks: ${processes}"

export PATH=$PROJDIR/software/CAT_pack/CAT_pack:$PATH

#################
#### RUN CAT ####
#################

cd $subdir 

# Annotate
CAT_pack contigs -c $contigs --nproc $processes -d $catdb/cat_db/db --out all --taxonomy_folder $catdb/cat_db/tax

# Add names to classifications and summarise
CAT_pack add_names -i all.contig2classification.txt -o all.contig2classification_named.txt -t $catdb/cat_db/tax --only_official
