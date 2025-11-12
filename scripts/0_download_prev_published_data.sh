#!/bin/bash -l

#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J prev_published_data
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Download previously published data
# Samples were selected and download links acquired via AncientMetagenomeDir

################
#### SET UP ####
################

## Change to the directory containing the script
[ -z $SLURM_SUBMIT_DIR ] || cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
published=${REFDIR}/published_metagenomes/ # Directory where the published metagenomes will be downloaded
outdir=${OUTDIR}

# Sample list
sample_list=$indir/human_dental_calculus.csv

[ -d $REFDIR ] || mkdir $REFDIR # Create directory for all references
[ -d $published ] || mkdir $published # Create directory for concatenated references

## Load modules and activate environment
conda activate oral_mb_evol

###########################
##### DOWNLOAD GENOMES ####
###########################

cd $published

# Get download links for input file
tail -n+2 ${sample_list} | while IFS="," read sample_name links md5s _
do
    # Get the links for forward and reverse into an array, same for corresponding md5sums
    IFS=';' read -r -a links_array <<< "$links"
    echo $links
    IFS=';' read -r -a md5_array <<< "$md5s"
    echo $md5s
    # Check that we have the same number of elements in both arrays
    if [[ ${#links_array[@]} -ne ${#md5_array[@]} ]]; then
        echo "Download links and md5sums don't match for $sample_name. Skip"
        continue
    fi
    # Check if the file is already downloaded
    if ([ -f ${sample_name}_F.fastq.gz ] && [ -f ${sample_name}_R.fastq.gz ]) || ([ -f ${sample_name}_1.fq.gz ] && [ -f ${sample_name}_2.fq.gz ])
    then
        echo "File ${sample_name} already downloaded. Skip"
        continue
    else
        echo "Downloading $sample_name"
    fi

    # In both cases, download and store md5sums in a file
    wget -O ${sample_name}_F.fastq.gz ${links_array[0]}
    echo -e "${md5_array[0]}\t${sample_name}_F.fastq.gz" >> md5sums.tsv

    wget -O ${sample_name}_R.fastq.gz ${links_array[1]}
    echo -e "${md5_array[1]}\t${sample_name}_R.fastq.gz"  >> md5sums.tsv
done

# Check md5sums
md5sum -c md5sums.tsv >> md5sums_output.txt

####################
####  METADATA  ####
####################

# Run FastQC
echo "RUNNING FASTQC"
fastqc *.gz -f fastq -o . -t $SLURM_NTASKS

# Aggregate reports in MultiQC
echo "RUNNING MULTIQC"
multiqc -f --interactive .

# Remove intermediate FastQC reports
rm *fastqc*

# Extract read count and length in a sepate file
# Read count table
echo "COLLECTING METADATA"
echo "sample,demult_F,demult_R" > ./read_count.csv
 
# Average length table
echo "sample,demult_F,demult_R" > ./read_length.csv
 
tail -n+2 $sample_list | while IFS=, read sample _ fwd_file rev_file study
do
    if [[ $study == "this_study" ]]
    then
        continue
    fi
    
    fwd=$(basename ${fwd_file%.f*.gz})
    rev=$(basename ${rev_file%.f*.gz})
    
    # Get number of reads from multiqc_data
    rc_F=$(grep "^$fwd" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
    rc_R=$(grep "^$rev" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $5}')
    
    # Add to file
    echo "$sample,$rc_F,$rc_R" >> ./read_count.csv
    
    # Get read length
    rl_F=$(grep "^$fwd" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $10}')
    rl_R=$(grep "^$rev" multiqc_data/multiqc_fastqc.txt | awk -v FS="\t" '{print $10}')

    echo "$sample,rl_F,rl_R" >> ./read_length.csv
done
