#!/bin/bash -l

#SBATCH -n 5
#SBATCH -t 3:00:00
#SBATCH -J find_input
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script finds the input files for the project from the file list,
# creates links to the demultiplexed fastq files and runs FASTQC

################
#### SET UP ####
################

## Change to the directory containing the script
[ -z $SLURM_SUBMIT_DIR ] || cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Demultiplexed dir: ${DEMUXDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR}
subdir=${OUTDIR}/0_all_raw_reads/ # Directory were all reads will be collected

refdir=${REFDIR} # Directory where the reference files are located
pubdir=${REFDIR}/published_metagenomes # Directory where the published metagenomes are located

mkdir -p $outdir
mkdir -p $subdir
mkdir -p ${scriptdir}/logs

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info

echo "Saving output to: ${subdir}"

echo "Sample list:"
echo ${sample_list}
cat $sample_list

echo "Global modules:"
module list
echo "Conda modules:"
conda list

if [ -z "$SLURM_NTASKS" ]; then
  n_proc=1
else
  n_proc=$SLURM_NTASKS
fi

echo "Number of tasks: ${n_proc}"

## Run in 3 threads per CPU
n_threads=$( echo "${n_proc} * 3" | bc )

cd $subdir/

#################
#### PROCESS ####
#################

#### Create symlinks to the demultiplexed fastq files
# Process samples in provided samples list - First filter for the newly generated data
tail -n+2 ${sample_list} | while IFS=, read -r sample species study; do 
    if [[ $study == "This Study" ]]
    then
        fwd_file=$DEMUXDIR/${sample}_F.fastq.gz
        rev_file=$DEMUXDIR/${sample}_R.fastq.gz
    else
        fwd_file=$pubdir/${sample}_F.fastq.gz
        rev_file=$pubdir/${sample}_R.fastq.gz
    fi
    # expand variables in fwd_file with envsubst
    if [ ! -f $fwd_file ] || [ ! -f $rev_file ]; then
        echo "Couldn't find input for $sample. Skip"
        continue
    fi
    echo "Creating symlinks for $sample"
    ln -s $fwd_file ${subdir}/${fwd_name}
    ln -s $rev_file ${subdir}/${rev_name}
done

####################
####  METADATA  ####
####################

# Run FastQC
echo "RUNNING FASTQC"
fastqc *.gz -f fastq -o . -t $n_threads

# Aggregate reports in MultiQC
echo "RUNNING MULTIQC"
multiqc -f --interactive .

# Remove intermediate FastQC reports
rm *fastqc*

# Extract read count and length in a sepate file
# Read count table
echo "COLLECTING METADATA"
echo "sample,demult_F,demult_R" > $OUTDIR/read_count.csv
 
# Average length table
echo "sample,demult_F,demult_R" > $OUTDIR/read_length.csv

tail -n+2 $sample_list | while IFS=, read sample _ _ 
do
    
    fwd=${sample}_F
    rev=${sample}_R
    
    # Get number of reads from multiqc_data
    rc_F=$(grep "^$fwd" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
    rc_R=$(grep "^$rev" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
    
    # Add to file
    echo "$sample,$rc_F,$rc_R" >> $OUTDIR/read_count.csv
    
    # Get read length
    rl_F=$(grep "^$fwd" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $10}')
    rl_R=$(grep "^$rev" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $10}')

    echo "$sample,$rl_F,$rl_R" >> $OUTDIR/read_length.csv
done
