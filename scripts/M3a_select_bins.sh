#!/bin/bash -l

#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J select_bins
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script selects bins for downstream analysis and collects metadata

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

# Input files and dirs
bindir=${outdir}M2_nfcore_mag/GenomeBinning/ 
bin_metadata=${outdir}M2_nfcore_mag/GenomeBinning/bin_summary.tsv
pydamagedir=${outdir}M2_nfcore_mag/Ancient_DNA/pydamage/
sample_metadata=${indir}/sample_metadata.csv

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

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

#############################
#### IDENTIFY FINAL MAGS ####
#############################

# For each sample, if bin refinement ran, select refined bins,
# if bin refinement didn't run (because there was no MaxBin2 output, select MetaBat2 bins

# Write a short function that takes a file list (output of 'ls') and the bin type
# and prints a table with four columns: sample, bintype, bin file name, and bin file path
print_bin_table () {
  readarray -t file_array <<< $1
  for file in ${file_array[*]}; do
    filename=${file##*/}
    echo "$sample,$2,$filename,$file" 
  done
}

# Store final bin paths in final_bins.csv
echo "sample,bin_type,bin,bin_path" > bins.csv

refined=0
only_maxbin=0
only_metabat=0
odd=0
missing=0

echo " ======= SELECTING BINS TO USE DOWNSTREAM ======= "
while read sample; do
  # Search for refined bins, if they exist add them to the final bins CSV
  refined_files=$(ls $bindir/DASTool/bins/MEGAHIT-**Refined-${sample}.*fa.gz 2> /dev/null | sort)
  if [[ -n $refined_files ]]; then
    echo "Found refined bins for $sample."
    print_bin_table "${refined_files}" "refined" >> bins.csv
    refined=$((refined+1))
  # If they are not there check if MaxBin2 or MetaBAT2 output is there (if both exist, refinement should have run)
  else
    maxbin_files=$(ls $bindir/MaxBin2/bins/MEGAHIT-MaxBin2-${sample}.*fa.gz 2> /dev/null | sort)
    metabat_files=$(ls $bindir/MetaBAT2/bins/MEGAHIT-MetaBAT2-${sample}.*fa.gz 2> /dev/null | sort )
    if [[ -n $maxbin_files && -n $metabat_files ]]; then
      echo "No refined bins for $sample, but both MetaBAT2 and MaxBin2 output exist."
      odd=$((odd+1)); fi
    if [[ -n $maxbin_files && ! -n $metabat_files ]]; then
      echo "No refined bins for $sample, because only MaxBin2 output exists."
      print_bin_table "${maxbin_files}" "MaxBin2" >> bins.csv
      only_maxbin=$((only_maxbin+1)); fi
    if [[ ! -n $maxbin_files && -n $metabat_files ]]; then
      echo "No refined bins for $sample, because only MetaBAT2 output exists."
      print_bin_table "${metabat_files}" "MetaBAT2" >> bins.csv
      only_metabat=$((only_metabat+1)); fi
    if [[ ! -n $maxbin_files && ! -n $metabat_files ]]; then
      echo "No bins for $sample."
      missing=$((missing+1)); fi
  fi
done < <(awk -F "," 'NR >1 {print $1}' $sample_list)

echo -e " $refined samples with refined bins,\n $only_maxbin samples with only MaxBin2 bins,\n $only_metabat sample with only MetaBAT2 bins,\n $odd samples without refined bins despite both binners running,\n $missing with no bins"

# Create symbolic links for all final bins
[ -d "BINS" ] || mkdir "BINS"

awk -F "," 'NR > 1 {print $0}' bins.csv | while IFS=, read _ _ filename path ; do
   cp $path ./BINS/$filename 
done

##################################
#### COLLECT CONTIG METADATA  ####
##################################

# Get metadata for contigs that ended up 

## Get which contigs ended up in which output bin (out of those selected in the previous step)
echo "contig,length,bin,sample" > contigs_per_bin.csv

total=$(wc -l bins.csv | cut -d " " -f1)
total=$((total - 1))
elapsed=0

echo " ======= IDENTIFY WHICH CONTIGS ENDED UP IN WHICH BINS  ======= "
while IFS="," read sample _ bin _; do
  zcat BINS/$bin | grep ">" | cut -d' ' -f1,4 | while IFS=" " read contig length; do
    echo "${sample}_${contig#>},${length/len=},$bin,${sample}" >> contigs_per_bin.csv
  done
  elapsed=$((elapsed+1))
  echo "${elapsed}/${total}: Extracting contigs from bin ${bin}"
done < <(tail -n +2  bins.csv)

echo "Total contigs found: $(wc -l contigs_per_bin.csv)"

echo " ======= COMBINE PYDAMAGE OUTPUT ======= "

## Combine and summarise all pydamage output
# Get headers
head -1 $(find ${pydamagedir}/analyze/**/pydamage_results/ -name "pydamage_results.csv" | head -1) > pydamage_results_combined.csv

total=$(find ${pydamagedir}/analyze/**/pydamage_results/ -name "pydamage_results.csv" | wc -l |  cut -d " " -f1)
elapsed=0

# First add sample name before contig label, to get unique identifiers matching contigs_per_bin.csv
find ${pydamagedir}/analyze/**/pydamage_results/ -name "pydamage_results.csv" | while read path; do
  tmp=${path#**/analyze/MEGAHIT-}
  sample=${tmp%/pydamage_results/pydamage_results.csv}
  elapsed=$((elapsed + 1))
  echo "${elapsed}/${total}: Processing pydamage file for $sample"
  # Filter contigs_per_bin.csv 
  awk -F, -v sample=$sample '$4==sample {print $1 ","}' contigs_per_bin.csv > contigs_${sample}.tmp
  # Append sample name and check if contig was used for the final bins from that sample
  awk -v sample="$sample" 'NR >1 {print sample "_" $0}' "$path" | grep -f contigs_${sample}.tmp >> pydamage_results_combined.csv
done

rm contigs_*.tmp # remove temporary files

# Combine the two tables
python $scriptdir/modules/combine_tables.py --sep "," --axis "columns" contigs_per_bin.csv pydamage_results_combined.csv > contig_metadata.csv

[ -f contig_metadata.csv ] && rm contigs_per_bin.csv pydamage_results_combined.csv # if output exists remove intermediate files

##############################
#### COLLECT BIN METADATA ####
##############################

echo " ======= COLLECT BIN METADATA ======= "

## Filter bin_summary.tsv from nfcore/mag to keep only selected bins, make csv instead of tsv
# Remove all the depth columns

awk -F '\t' '
  NR==1 {
    # Go through headers and identify those containing the word "Depth"
    for (i = 1; i <= NF; i++) {
      if ($i !~ /Depth/ && i > 1) {
        keep_indices[i] = 1
        keep_headers[i] = $i
      }
    }
    # Print the "keep" headers and "Depth"
    for (i in keep_indices) {
      printf "%s,", keep_headers[i]
    }
    printf "\n"
  }
  NR > 1 {
    for (i in keep_indices) {
      printf "%s,", $i
    }
    printf "\n"
  }' $bin_metadata > bin_summary.tmp

## Filter bin_summary.tmp from nfcore/mag to keep only selected bins, make csv instead of tsv
head -1 bin_summary.tmp | sed 's/\t/,/g'  > bin_summary.csv
grep -f <(awk -F "," 'NR > 1{print $3}' bins.csv | sed 's/\.gz//' ) bin_summary.tmp | sed 's/\t/,/g' |sed 's/.fa,/.fa.gz,/g'  >> bin_summary.csv

## Get info about the samples they originated from
echo "bin,Sample,Species,Common.name,Study,Domestication,Location,Country,Period" > sample_metadata.csv

while IFS="," read sample _ bin _; do
  # Find sample in metadata, select relevant column
  echo ${bin},$( grep -a "^${sample}," $sample_metadata | csvcut -c 1,2,3,4,5,6,7,8 | uniq) >> sample_metadata.csv
done < <(tail -n +2 bins.csv )

## Summarise python info per bin

echo " ======= SUMMARISE PYDAMAGE OUTPUT PER BIN ======= "

python $scriptdir/modules/summarise_pydamage_res.py \
            --damage_file contig_metadata.csv \
            --groupvar "bin" # Summarise by bin

# Make csv
awk -F '\t' 'BEGIN {OFS=","} {print $1, $2, $3, $4, $5, $6, $7}' pydamage_summary.tsv > pydamage_summary.csv
rm pydamage_summary.tsv

# Combine the tables
python $scriptdir/modules/combine_tables.py --sep "," --axis "columns" sample_metadata.csv bin_summary.csv pydamage_summary.csv <(awk -F "," 'BEGIN {OFS=","} {print $3, $2, $4}' bins.csv) > bin_metadata.csv

[ -f bin_metadata.csv ] && rm sample_metadata.csv bin_summary.csv bin_summary.tmp pydamage_summary.csv bins.csv # if output exists remove intermediate files
