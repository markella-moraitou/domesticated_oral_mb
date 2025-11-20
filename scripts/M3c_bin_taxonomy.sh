#!/bin/bash -l

#SBATCH -p main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH -t 1-00:00:00
#SBATCH --mem=300GB
#SBATCH -J bin_taxonomy
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script runs GTDB-Tk on high quality drafts

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

magdir=${outdir}M3_MAG_metadata/MQ_HQ_BINS/ # Input directory with 'good' MAGs
subdir=${outdir}M3_MAG_metadata/GTDB_tk # Output subdirectory for GTDB taxonomy

mkdir -p $subdir

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate gtdbtk # Had to create a separate environment due to conflicts

## Print some info
echo "Taking input from: ${magdir}"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z "$SLURM_NTASKS" ] && [ -z "$SLURM_CPUS_PER_TASK" ]; then
  processes=2
else
  processes=$(echo "$SLURM_NTASKS * $SLURM_CPUS_PER_TASK" | bc)
fi

echo "Number of process: $processes"
cd $subdir

##################
#### CLASSIFY ####
##################

export GTDBTK_DATA_PATH=$REFDIR/gtdbtk_r220_data/release220/

# Run classify_wf
mkdir -p "classify_wf"

echo "RUNNING classify_wf"
gtdbtk classify_wf --cpu $processes --pplacer_cpus 1 --extension gz \
    --genome_dir $magdir --mash_db $GTDBTK_DATA_PATH \
    --out_dir "classify_wf"

# Convert trees to iTOL format
find classify_wf -name "gtdbtk.*.tree" | while read tree; do
    echo "Converting $tree to iTOL format"
    gtdbtk convert_to_itol --input $tree --output ${tree%.tree}_itol.tree
done

# Run de_novo_wf for bacteria
mkdir -p "de_novo_wf"

echo "RUNNING de_novo_wf FOR BACTERIA"
gtdbtk de_novo_wf --cpu $processes --extension gz --bacteria --outgroup_taxon p__Chloroflexota \
    --genome_dir $magdir --out_dir "de_novo_wf" \
    --gtdbtk_classification_file "classify_wf/gtdbtk.bac120.summary.tsv"

# Convert trees to iTOL format
find de_novo_wf -name "gtdbtk.*.tree" | while read tree; do
    echo "Converting $tree to iTOL format"
    gtdbtk convert_to_itol --input $tree --output ${tree%.tree}_itol.tree
done

# Run de_novo_wf for archaea
mkdir -p "de_novo_wf_archaea"

echo "RUNNING de_novo_wf FOR ARCHAEA"
gtdbtk de_novo_wf --cpu $processes --extension gz --archaea --outgroup_taxon p__Altiarchaeota \
    --genome_dir $magdir --out_dir "de_novo_wf_archaea" \
    --gtdbtk_classification_file "classify_wf/gtdbtk.ar53.summary.tsv"

# Convert trees to iTOL format
find de_novo_wf_archaea -name "gtdbtk.*.tree" | while read tree; do
    echo "Converting $tree to iTOL format"
    gtdbtk convert_to_itol --input $tree --output ${tree%.tree}_itol.tree
done

##########################
#### COLLECT METADATA ####
##########################

conda activate oral_mb_evol

# Add this info to hq_bin_metadata.csv
# First make sure bins match and turn to csv
sed 's/.fa\t/.fa.gz\t/g' classify_wf/gtdbtk.bac120.summary.tsv | sed 's/,/;/g' | sed 's/\t/,/g' > gtdbtk.bac120.summary.csv
sed 's/.fa\t/.fa.gz\t/g' classify_wf/gtdbtk.ar53.summary.tsv | sed 's/,/;/g' | sed 's/\t/,/g' > gtdbtk.ar53.summary.csv

# Combine archaea and bacteria summaries
python $scriptdir/modules/combine_tables.py --sep "," --axis "rows" gtdbtk.bac120.summary.csv gtdbtk.ar53.summary.csv > gtdbtk.combined.summary.csv

# Add to HQ bin metadata
cd ..
python $scriptdir/modules/combine_tables.py --sep "," --axis "columns" mq_hq_bin_metadata.csv ${subdir}/gtdbtk.combined.summary.csv > mq_hq_bin_metadata.tmp

[ $(wc -l < mq_hq_bin_metadata.tmp) -eq $(wc -l < mq_hq_bin_metadata.csv ) ] && mv mq_hq_bin_metadata.tmp mq_hq_bin_metadata.csv
