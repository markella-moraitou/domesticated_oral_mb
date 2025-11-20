#!/bin/bash -l

#SBATCH -p memory
#SBATCH -n 100
#SBATCH -t 10:00:00
#SBATCH --mem=800GB
#SBATCH -J checkM
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script runs checkM and identifies high quality drafts

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

subdir=${outdir}M3_MAG_metadata/ # Output subdirectory for MAG assembly

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
  processes=2
else
  processes=$(echo "$SLURM_NTASKS * 2" | bc)
fi

echo "Number of tasks: ${SLURM_NTASKS}"

cd $subdir

####################
#### Run CheckM ####
####################

[ -d "CheckM" ] || mkdir "CheckM"

checkm lineage_wf ./BINS ./CheckM -t $processes -x ".fa.gz" --tab_table -f ./CheckM/checkm_output.tsv 

#### Collect metadata

echo "Bin Id,Marker lineage,# genomes,# markers,# marker sets,0,1,2,3,4,5+,Completeness,Contamination,Strain heterogeneity" > bins_checkM.csv

awk -F "\t" 'BEGIN {OFS=","} NR > 1 {
    printf "%s.fa.gz", $1 
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}' CheckM/checkm_output.tsv >> bins_checkM.csv

# Combine with existing bin metadata
cp bin_metadata.csv bin_metadata.tmp
python $scriptdir/modules/combine_tables.py --sep "," --axis "columns" bin_metadata.tmp bins_checkM.csv > bin_metadata.csv

[ -f bin_metadata.csv ] && rm bin_metadata.tmp bins_checkM.csv

###############################
#### IDENTIFY MQ & HQ MAGS ####
###############################

min_completeness=50
max_contamination=10

#### Select MAGs that meet the minimum completeness and maximum contamination standards
python $scriptdir/modules/filter_tables.py --input bin_metadata.csv --sep $',' --expression "Completeness > ${min_completeness} and Contamination < ${max_contamination}" > mq_hq_bin_metadata.csv

# Create a new folder with the HQ bins
mkdir -p "MQ_HQ_BINS"

awk -F "," 'NR > 1 {print $1}' mq_hq_bin_metadata.csv | while read bin ; do
  ln -s $subdir/BINS/$bin ./MQ_HQ_BINS/$bin
done
