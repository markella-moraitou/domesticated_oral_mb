#!/bin/bash -l

#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J cleanup
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Cleaning up and merging reads (both this study and previously available)
# Trims poly-G and poly-X tails (resulting from NovaSeq two-colour chemistry), low-quality bases and Ns (sliding windows of 3)
# Removes adapters and barcodes
# Also, filters out sequences less than 30 bp

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

outdir=${OUTDIR}
demultdir=${OUTDIR}/0_all_raw_reads/ # Directory with demultiplexed sequences
subdir=${OUTDIR}/1_cleaned_seqs/ # Output subdirectory for the output of this specific script
repdir=${OUTDIR}/1_cleaned_seqs/fastq_reports/ # Output subdirectory for the fastp reports

[ -d $subdir ] || mkdir $subdir # Create subdir if not there
[ -d ${scriptdir}/logs ] || mkdir ${scriptdir}/logs # Create directory for script logs if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info

echo "Saving fastq output to: ${subdir}"
echo "Saving reports to: ${repdir}"

echo "Sample list:"
echo ${sample_list}
cat $sample_list

echo "Conda modules:"
conda list

echo "Number of tasks: ${SLURM_NTASKS}"

## Run in 3 threads per CPU
n_proc=$( echo "${SLURM_NTASKS} * 3" | bc )

cd $subdir/

#################
#### PROCESS ####
#################

#### Adapter removal, quality trimming and merging with fastp ####
# For samples from this study, also remove barcodes by trimming first 7 nucleotides from both forward and reverse reads

tail -n+2 ${sample_list} | while IFS=, read -r sample species study
do
  fwd_file=${sample}_F.fastq.gz
  rev_file=${sample}_R.fastq.gz
  if [ -f ${sample}_merged.fastq.gz ]; then
    echo "Output files for $sample already exist, skipping"
    continue
  fi
  # Exclude samples not from this study - they will be processed a bit differently
  if [[ $study == "This Study" ]]
  then
    # Process with fastq
    # --trim_poly_g: removes low quality reads (both sides, mean quality 30, window size 3)and poly-G tails (due to NovaSeq 2-colour chemistry)
    # -3 --cut_tail_window_size 3 --cut_tail_mean_quality 30: trim seqs with qual < 30 in sliding windows from 3'
    # --adapter_fasta: specify adapter sequences (Ns not allowed, so I need to specify one for each index)
    # --trim_front1 7 and --trim_front2 7: remove barcodes by trimming the first 7 nucleotides from both forward and reverse
    # --merge, --overlap_len_require 11 and --overlap_diff_limit 3: merge F and R if there is at least 11 of overlap with at most 3 mismatches (probability of this happening at random <1%)
    # --length_required 30: filters out reads shorter than 30
    # Saving merged reads in one file
    # and unpaired F and R in separate files (regardless if they could not be merged or if the other mate failed the filters)
    echo "Processing $sample with barcode trimming"
    fastp \
      --merge \
      --in1 $demultdir/$fwd_file \
      --in2 $demultdir/$rev_file \
      --merged_out ${subdir}/${sample}_merged.fastq.gz \
      --unpaired1 ${subdir}/${sample}_unpaired_F.fastq.gz \
      --out1 ${subdir}/${sample}_unmerged_F.fastq.gz \
      --unpaired2 ${subdir}/${sample}_unpaired_R.fastq.gz \
      --out2 ${subdir}/${sample}_unmerged_R.fastq.gz \
      --trim_poly_g \
      --trim_poly_x \
      -3 --cut_tail_window_size 3 --cut_tail_mean_quality 30 \
      --adapter_fasta ${indir}/adapter_seqs.fasta \
      --trim_front1 7 \
      --trim_front2 7 \
      --overlap_len_require 11 \
      --overlap_diff_limit 3 \
      --length_required 30 \
      --dedup \
      --json ${repdir}/${sample}.fastp.json  \
      --html ${repdir}/${sample}.fastp.html  \
      --thread $n_proc
  else
    echo "Processing $sample without barcode trimming"
    fastp \
      --merge \
      --in1 $demultdir/$fwd_file \
      --in2 $demultdir/$rev_file \
      --merged_out ${subdir}/${sample}_merged.fastq.gz \
      --unpaired1 ${subdir}/${sample}_unpaired_F.fastq.gz \
      --out1 ${subdir}/${sample}_unmerged_F.fastq.gz \
      --unpaired2 ${subdir}/${sample}_unpaired_R.fastq.gz \
      --out2 ${subdir}/${sample}_unmerged_R.fastq.gz \
      --trim_poly_g \
      --trim_poly_x \
      -3 --cut_tail_window_size 3 --cut_tail_mean_quality 30 \
      --adapter_fasta ${indir}/adapter_seqs.fasta \
      --overlap_len_require 11 \
      --overlap_diff_limit 3 \
      --length_required 30 \
      --dedup \
      --json ${repdir}/${sample}.fastp.json  \
      --html ${repdir}/${sample}.fastp.html  \
      --thread $n_proc
  fi
done

##################################
#### RUN QC AND SAVE METADATA ####
##################################

## The following bit generates a MultiQC report and then extracts info from it to add to read_count.csv and read_length.csv
## These to files will be updated after every step, to keep track of the pre-processing pipeline

# Since I have to run it twice, first to extract columns "FastQC_mqc-generalstats-fastqc-total_sequences" from a few different MultiQC output tables and add them to read_count.csv and
# then to extract columns "FastQC_mqc-generalstats-fastqc-avg_sequence_length" (again from a few tables) and add them to read_length.tsv,
# I will define the following actions in a function that I can call

cd $subdir

## Get QC reports for each type of output

output_categories=(merged unmerged unpaired)

for category in "${output_categories[@]}"
do  
    # Run FastQC
    fastqc *_${category}*.fastq.gz -f fastq -o . -t $SLURM_NTASKS 
    # Generate MultiQC report
    multiqc *.zip -f --interactive --title ${category}
    # Remove FastQC files
    rm -f *fastqc*
done

# Define function
get_metadata() {
  local out_table="$1"
  local data_col="$2"

  # Do the entire thing seperately for the three output categories (merged, unmerged, unpaired)
  # This means that a few columns will be added in the out_table
  output_categories=(merged unmerged unpaired)
  
  for category in "${output_categories[@]}"
  do
    ## Extract data from MultiQC general stats table 
    
    # If these are unmerged or unpaired files, split the output in forward and reverse
    if [[ $category == "unmerged" || $category == "unpaired" ]]
    then 
      # Get info on forward reads, save on a temporary file, while normalising sample names (just sample name without suffix)
      awk 'NR==1 || $1 ~ /_F/ {print $0}' ${category}_multiqc_report_data/multiqc_fastqc.txt | \
       sed "s/_${category}_.//g" > ${category}_F_general_stats.tmp
      # Set data_table and out_col arguments for multiqc_to_csv.py custom script
      data_table=${category}_F_general_stats.tmp # The table from which data will be extracted
      out_col=${category}_F # The name of the column that will be created as an output
      # Run python script
      # Arg 1: data_table (file path; the tsv file from which data will be extracted)
      # Arg 2: out_table (file path; the csv file where a column will be added)
      # Arg 3: data_col (string; the name of the column in data_table that will be extracted)
      # Arg 4: out_col (string; the name of the new column in out_table)
      python $scriptdir/modules/multiqc_to_csv.py "${data_table}" "${out_table}" "${data_col}" "${out_col}"
      # Same for reverse reads
      awk 'NR==1 || $1 ~ /_R/ {print $0}' ${category}_multiqc_report_data/multiqc_fastqc.txt | \
       sed "s/_${category}_.//g" > ${category}_R_general_stats.tmp
      data_table=${category}_R_general_stats.tmp # The table from which data will be extracted
      out_col=${category}_R # The name of the column that will be created as an output
      # Run python script
      python $scriptdir/modules/multiqc_to_csv.py "${data_table}" "${out_table}" "${data_col}" "${out_col}"
    # Else, if it is a merged file the process will be done only one
    else
      # Still need to get a table with the correct sample names
      sed "s/_${category}//g" ${category}_multiqc_report_data/multiqc_fastqc.txt > ${category}_general_stats.tmp
      # Set data_table and out_col arguments
      data_table=${category}_general_stats.tmp # The table from which data will be extracted
      out_col=${category} # The name of the column that will be created as an output
      # Run python script
      python $scriptdir/modules/multiqc_to_csv.py "${data_table}" "${out_table}" "${data_col}" "${out_col}"
    fi
  rm -f *.tmp
  done
}

# Run for read count
get_metadata ${outdir}/read_count.csv "Total Sequences"

# Run for read length
get_metadata ${outdir}/read_length.csv "avg_sequence_length"
