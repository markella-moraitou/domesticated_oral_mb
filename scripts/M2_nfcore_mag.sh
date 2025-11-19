#!/bin/bash -l

# This script uses the nfcore-mag pipeline to construct and validate metagenome assembled genomes

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR} # Output subdirectory for the output of this specific script

subdir=${outdir}M2_nfcore_mag/ # Output subdirectory for MAG assembly
seqdir=${outdir}3_mapped_seqs/ # Input directory with mapped sequences
contigdir=${outdir}M1_contig_assemblies/ # Input directory with assembled contigs

mkdir -p $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment

#conda create --name env_nf nextflow
conda activate env_nf

# Fix SSL certificate path issues - use environment variables to determine conda environment location
# Get the conda environment path dynamically
CONDA_ENV_PATH=$(conda info --base)/../.conda/envs/env_nf
if [ ! -d "$CONDA_ENV_PATH" ]; then
    # Fallback: try to find the actual conda environment path
    CONDA_ENV_PATH=$(conda info --envs | grep "env_nf" | awk '{print $NF}' | head -1)
fi

# Set SSL certificate environment variables using the detected conda environment path
export SSL_CERT_FILE="${CONDA_ENV_PATH}/ssl/cacert.pem"
export REQUESTS_CA_BUNDLE="${CONDA_ENV_PATH}/ssl/cacert.pem"
export CURL_CA_BUNDLE="${CONDA_ENV_PATH}/ssl/cacert.pem"

# Verify the SSL certificate file exists
if [ ! -f "$SSL_CERT_FILE" ]; then
    echo "Warning: SSL certificate file not found at $SSL_CERT_FILE"
    echo "Falling back to system certificates"
    unset SSL_CERT_FILE REQUESTS_CA_BUNDLE CURL_CA_BUNDLE
fi

# Activate the modules, you can also choose to use a specific version with e.g. `Nextflow/21.10`.
module load python/3.9.5 nextflow nf-core nf-core-pipelines
module load PDC singularity

## Print some info
echo "Taking input from: ${seqdir}"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

###############
#### SETUP ####
###############

cd $subdir

## Create sample sheet
if [ ! -f ./samplesheet.csv ] || [ ! -f ./assemblies.csv ]
then
  echo "sample,group,short_reads_1,short_reads_2,long_reads" > ./samplesheet.csv
  echo "id,group,assembler,fasta" > ./assemblies.csv
  
  tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _ _
  do
    seq_input=${seqdir}/${sample}_unmapped.fastq.gz
    contig_input=${contigdir}/${sample}_final_contigs.fa
    if [[ ! -f $seq_input ]] || [[ ! -f $contig_input ]]
    then
      echo "Input missing for ${sample}. Skipping..."
      continue
    fi
    # Calculate longest contig. If shorter than 500kb, skip.
    short=1
    while read i
    do
      if [[ ${i} -gt 500 ]]
      then
        short=0
        echo "Adding ${sample}"
        echo -e "$sample,0,$seq_input,," >> ./samplesheet.csv
        echo -e "$sample,0,MEGAHIT,$contig_input" >> ./assemblies.csv
        break
      fi
    done < <(grep -v '>' $contig_input  2> /dev/null | awk '{print length($0)}' 2> /dev/null)
    if [[ $short -eq 1 ]]
    then
      echo "$sample contigs too short. Skipping..."
    fi
  done
fi

######################
#### RUN PIPELINE ####
######################

nextflow run nf-core/mag -r 3.2.1 -profile pdc_kth --project $PROJ \
                          --input ./samplesheet.csv --assembly_input ./assemblies.csv --single_end --outdir $subdir\
                          --skip_gtdbtk \
                          --skip_prodigal --skip_prokka --skip_metaeuk --skip_concoct \
                          --busco_db ${REFDIR}/bacteria_odb10 --busco_clean \
                          --binning_map_mode own \
                          --refine_bins_dastool --postbinning_input "both" \
                          --ancient_dna  --pydamage_accuracy 0.2 \
                          -resume

##############
#### GZIP ####
##############

find GenomeBinning/ -type f -name "*.fa" -exec pigz -p 2 {} \;
