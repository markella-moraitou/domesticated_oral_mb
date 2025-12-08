#!/bin/bash -l

#SBATCH -p main
#SBATCH -n 200
#SBATCH -t 15:00:00
#SBATCH -J mapping
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script takes deduplicated reads and
# maps against the human and host genomes to maintain only unmapped reads
# Then runs decOM to assess unmapped reads

################
#### SET UP ####
################

## Change to the directory containing the script
[ -z $SLURM_SUBMIT_DIR ] || cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR} # Output subdirectory for the output of this specific script
subdir=${outdir}/3_mapped_seqs/ # Output subdirectory for the mapping reads

conc_refdir=${REFDIR}/concatenated_references/ # Directory to store concatenated references (host-human)

mkdir -p $subdir # Create subdirectory if not there

## Instructions to install decOM - More in https://github.com/CamilaDuitama/decOM.git
## Installing decOM
#mkdir -p ${PROJDIR}/software

#cd ${PROJDIR}/software
#git clone https://github.com/CamilaDuitama/decOM.git
#cd decOM
#conda env create -n decOM --file environment.yml
#conda deactivate

# Downloading sources
#wget https://zenodo.org/record/6513520/files/decOM_sources.tar.gz
#tar -xf decOM_sources.tar.gz

export PATH=${PROJDIR}/software/decOM:${PATH}

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
#module load bioinfo-tools bwa samtools BEDTools FastQC MultiQC/1.12 miniconda3
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${outdir}/2_deduplicated_seqs"
echo "Saving mapped fastq to: ${subdir}"

echo "Sample list:"
echo ${sample_list}

echo "Global modules:"
module list
echo "Conda modules:"
conda list

echo "Number of tasks: ${SLURM_NTASKS}"
cd $subdir

#################
#### MAPPING ####
#################

# This requires the host_genome_per_species.csv file
# Check if it exists, if not exit
if [ ! -f ${outdir}/concatenated_genome_per_species.csv ]
then
  echo "concatenated_genome_per_species.csv not found! Exiting..."
  exit
fi

cd $subdir

echo "Aligning based on ${outdir}/concatenated_genome_per_species.csv"
cat ${outdir}/concatenated_genome_per_species.csv

tail -n+2 $sample_list | while IFS=, read -r sample species _ _
do
  species=${species// /"_"} # Replace whitespace with underscore
  
  # Find input file in the fastq output
  input=${outdir}/2_deduplicated_seqs/${sample}_deduped.fastq.gz
   
  # Check if output file exists already, if yes skip
  if [[ ! -f ${sample}_unmapped.fastq  ]] && [[ -f ${sample}_unmapped.fastq.gz ]]
  then
    echo "Output for ${sample} exists, skipping"
    continue
  fi
  
  # Check if input files exists if not, print a message and continue
  if [[ ! -f $input ]]
  then
    echo "${input} does not exist"
    continue
  fi
  
  # Get reference from ${outdir}/concatenated_genome_per_species.csv 
  ref=$(awk -F "," -v species=$species '$1==species {print $2}' ${outdir}/concatenated_genome_per_species.csv)
  echo "Aligning $sample to $(basename ${ref%.gz} )"
  # Align with bwa aln -- using the parameters suggested in the nfcore-eager pipeline
  # Look for unzipped files, as the reference genomes were indexed before they were compressed
  bwa aln -n 0.04 -k 2 -l 1024 -o 2 -t $SLURM_NTASKS ${ref%.gz} $input > ${sample}_alignment.sai
  # Convert .sai to .sam and then to a sorted .bam
  bwa samse -r "@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina" ${ref%.gz} ${sample}_alignment.sai ${input} | samtools sort - -@ $SLURM_NTASKS -O bam > ${sample}.bam
  # Generate .csi index
  samtools index -c ${sample}.bam
  # Export unmapped reads as a fastq file
  samtools view -Sb -f 4 -@ 10 ${sample}.bam > ${sample}_unmapped.bam
  samtools fastq ${sample}_unmapped.bam > ${sample}_unmapped.fastq
  # Do the same for mapped reads
  samtools view -Sb -F 4 -@ 10 ${sample}.bam > ${sample}_mapped.bam
  samtools fastq ${sample}_mapped.bam > ${sample}_mapped.fastq
  # Compressed output
  pigz -p $SLURM_NTASKS ${sample}_unmapped.fastq
  pigz -p $SLURM_NTASKS ${sample}_mapped.fastq
done

##################################
#### RUN QC AND SAVE METADATA ####
##################################

## For unmapped reads ##
# Run FastQC
fastqc *_unmapped.fastq.gz -f fastq -o . -t $SLURM_NTASKS 

# Generate MultiQC report
multiqc *.zip -f --interactive --title unmapped

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_unmapped//g" unmapped_multiqc_report_data/multiqc_fastqc.txt > unmapped_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py unmapped_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "unmapped"
python $scriptdir/modules/multiqc_to_csv.py unmapped_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "unmapped"

## For mapped reads ##
# Run FastQC
fastqc *_mapped.fastq.gz -f fastq -o . -t $SLURM_NTASKS 

# Generate MultiQC report
multiqc *.zip -f --interactive --title mapped

# Remove FastQC files
rm -f *fastqc*

# Get correct sample names
sed "s/_mapped//g" mapped_multiqc_report_data/multiqc_fastqc.txt > mapped_general_stats.tmp

# Run python script
python $scriptdir/modules/multiqc_to_csv.py mapped_general_stats.tmp ${outdir}/read_count.csv "Total Sequences" "mapped"
python $scriptdir/modules/multiqc_to_csv.py mapped_general_stats.tmp ${outdir}/read_length.csv "avg_sequence_length" "mapped"

rm *tmp

## Get stats on host-mapped, human-mapped and PhiX-mapped

# First get stats per chromosome
ls ${subdir}*_mapped.bam | while read bam
do
  sample=$(basename ${bam%.bam})
  samtools idxstats $bam > ${sample}_idxstats.txt
done

# Summarise mapping stats per genome using in-house script (generates genome_mapping_stats.tsv output in the same folder)
python ${scriptdir}/modules/genome_mapping_stats.py "_idxstats.txt" ${outdir}/reference_genome_contig_list.csv

## Source tracking with decOM ##

# Activate environment
conda activate decOM

bash ${scriptdir}/modules/source_tracking_decOM.sh "_unmapped.fastq.gz" 5
