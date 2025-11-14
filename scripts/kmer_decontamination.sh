#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J kmer_decontam
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script uses k-mer similarity to decontaminate metagenomic samples using aKmerBroom

## Change to the directory containing the script
cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

export indir=${PROJDIR}/input/ # Directory where the provided input is located
export scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

export outdir=${OUTDIR} # Output subdirectory for the output of this specific script
export subdir=${outdir}/kmer_decontam/ # Output subdirectory for aKmerBroom output
export mapdir=${outdir}/2_mapped_seqs/ 

export script_path=$PROJDIR/software/aKmerBroom/

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Determine samples to be run
export sample_list=${indir}/file_list_test.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol
## Manually change PYTHONPATH
export PYTHONPATH=/path/to/desired/python/packages:$PYTHONPATH

## Print some info
echo "Taking input from: ${indir}"
echo "Saving mapped fastq to: ${subdir}"

echo "Global modules:"
module list

echo "Conda modules:"
conda list

echo "Number of tasks: ${SLURM_NTASKS}"

################
#### SET UP ####
################

cd $subdir

# Download ancient_kmers.bloom if not there
[ -f "ancient_kmers.bloom" ] || wget https://zenodo.org/record/7587160/files/ancient_kmers.bloom -O ancient_kmers.bloom

#################
#### PROCESS ####
#################

#### Write function ####
# To parallise this, I will create subdirectories for each job, run aKmerBroom, run, then copy output to parent directory

aKmerBroom() {
  sample=$1
  input="${mapdir}/${sample}_unmapped.fastq.gz"
  if [ ! -f $input ]
  then
    echo "Input file for ${sample} doesn't exist"
  elif [ -f ${sample}_annotated_with_anchor_kmers.fastq.gz ]
  then
    echo "Output for ${sample} already exists, skipping"
  else
  
    # Make new subdirectory for this sample
    [ -d $sample ] && rm -rf $sample # Remove subdir if already there
    mkdir $sample
    cd $sample
    
    [ -d "data" ] || mkdir "data"

    # Get symbolic link to the scripts and data we need
    ln -s $script_path/akmerbroom.py
    ln -s $script_path/scripts
    ln -s $subdir/ancient_kmers.bloom data/ancient_kmers.bloom
    
    # Unzip input file and save as 'unknown_reads.fastq'
    echo "Unzipping $sample file"
    zcat $input > data/unknown_reads.fastq
    
    # Run aKmerBroom
    echo "Running aKmerBroom on $sample"
    python akmerbroom.py --ancient_bloom
    
    # Rename output
    mv output/annotated_reads_with_anchor_kmers.fastq $subdir/${sample}_annotated_with_anchor_kmers.fastq
    mv output/annotated_reads.fastq $subdir/${sample}_annotated.fastq
    
    cd $subdir

    pigz ${sample}_annotated_with_anchor_kmers.fastq
    pigz ${sample}_annotated.fastq
  fi
}

export -f aKmerBroom

#### Run aKmerBroom ####

# First get sample list
samples=($(awk -F "," 'NR>1 {print $1}' $sample_list))

# Run in parallel
parallel -j ${SLURM_NTASKS} aKmerBroom ::: ${samples[*]} | tee parallel_log.txt

##################################
#### RUN QC AND SAVE METADATA ####
##################################

module load bioinfo-tools python/3.8.7 FastQC MultiQC 

# Run FastQC
fastqc *_annotated_with_anchor_kmers.fastq.gz -f fastq -o . -t $SLURM_NTASKS 

# Generate MultiQC report
multiqc *.zip -f --interactive --title decontaminated

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_annotated_with_anchor_kmers//g" decontaminated_multiqc_report_data/multiqc_fastqc.txt > decontaminated_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py decontaminated_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "decontaminated"
python $scriptdir/modules/multiqc_to_csv.py decontaminated_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "decontaminated"

rm *tmp

## Source tracking with decOM ##

# Activate environment
conda activate decOM

bash ${scriptdir}/modules/source_tracking_decOM.sh "_annotated_with_anchor_kmers.fastq.gz" $SLURM_NTASKS
