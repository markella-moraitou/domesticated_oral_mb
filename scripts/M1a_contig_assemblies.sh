#!/bin/bash -l

#SBATCH -n 100
#SBATCH -t 10:00:00
#SBATCH -J megahit
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
mapdir=${outdir}3_mapped_seqs/ # Input directory with mapped sequences

mkdir -p $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${outdir}/3_mapped_seqs/"
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

############################
#### BINNING & ASSEMBLY ####
############################

cd $subdir 

tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _
do
  # Find short reads
  input="${mapdir}/${sample}_unmapped.fastq.gz"
  if [ ! -f $input ]
  then
    echo "Input file for ${sample} doesn't exist"
    continue
  fi
  
  if [ -f $subdir/${sample}_final_contigs.fa ]
  then
    echo "Output for ${sample} already exists. Skipping..."
    continue
  fi
  
  megahit -r $input -t $processes -o $subdir/${sample}

  # If any contigs were generated, move them to main output directory
  output=$subdir/${sample}/final.contigs.fa
  
  if [ -s $output ]
  then
    mv $output $subdir/${sample}_final_contigs.fa
  fi
done

echo "ASSEMBLY DONE"

################################
#### KEEP ONLY LONG CONTIGS ####
################################

# Keep only contigs longer than 1000
ls *_final_contigs.fa | while read contigs
do
  sample=${contigs%_final_contigs.fa}
  awk '/^>/ {header=$0} !/^>/ {if (length($0) > 1000) print header ORS $0}' ${contigs} > ${sample}_final_contigs_1000bp.fa
done

##########################
#### COLLECT METADATA ####
##########################

# Collect metadata about all assembled contigs, based on info on FASTA headers
echo -e "sample\tcontig\tflag\tmulti\tlen" > contig_info.txt

ls *_final_contigs.fa | while read contigs
do
  sample=${contigs%_final_contigs.fa}
  awk -v sample=$sample '/^>/ {printf "%s\t", sample
            printf "%s\t", substr($1, 2);
            for (i=2; i<=NF; i++) 
              { if ($i ~ /^flag=/) { printf "%s\t", substr($i, 6) } 
                else if ($i ~ /^multi=/) { printf "%s\t", substr($i, 7) } 
                else if ($i ~ /^len=/) { printf "%s\n", substr($i, 5) } }}' $contigs >> contig_info.txt
done

# Summarise contig info per sample
Rscript $scriptdir/modules/contig_stats_per_sample.R contig_info.txt
