#!/bin/bash -l

#SBATCH -p long
#SBATCH -n 10
#SBATCH -t 4-00:00:00
#SBATCH -J get_abundances
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# This script maps all samples to a concatenated file of dereplicated MAGs to get abundances

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

subdir=${outdir}M4_dereplicated_bins/mapping # Output subdirectory for MAG mapping
magdir=${outdir}M4_dereplicated_bins/dereplicated_genomes # Input directory with dereplicated 'good' MAGs
maginf=${outdir}M3_MAG_metadata # Input dir with metadata for the MAGs
readdir=${outdir}3_mapped_seqs/ # Directory where the unmapped reads are stored

mkdir -p $subdir

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
  processes=1
else
  processes=$SLURM_NTASKS
fi

echo "Number of tasks: ${processes}"

cd $subdir

#############################
#### MAP SAMPLES TO MAGS ####
#############################

# Check if file exists and delete if so
[ -f dereplicated_MAGs.fa ] && rm dereplicated_MAGs.fa
[ -f dereplicated_MAGs.fa.gz ] && rm dereplicated_MAGs.fa.gz

ls ${magdir}/*.fa.gz | while read mag 
do
  sample=${mag#${magdir}/MEGAHIT-*-}
  sample=${sample%.[0-9]*.fa.gz}
  zcat $mag | sed "s/>/>${sample}_/g" >> dereplicated_MAGs.fa
done

pigz -p ${processes} dereplicated_MAGs.fa

#### INDEX
echo "Indexing..."
kma index -i ${subdir}/dereplicated_MAGs.fa.gz -o database 2> index.log

#### MAP
tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _
do
    echo "Mapping sample ${sample} (${species})"
    # Find mapped reads files
    read_file=$(ls ${readdir}/${sample}_unmapped.fastq.gz 2> /dev/null | head -n 1)
    
    # Run KMA to map reads to contigs and generate SAM file
    echo "Mapping..."
    kma -i $read_file -o ${subdir}/${sample}_kma \
        -t_db database -t ${processes} -nc -1t1
    
    echo "Finished mapping sample ${sample}"
done

######################
#### Collect info ####
######################

# Filter contig metadata to only include length above 1000bp (shorter were not used for binning)
awk -F "," '$2 > 1000 {print $0}' ${maginf}/contig_metadata.csv > contig_metadata_lt1000.csv

# Log messages
echo "Calculating abundances from mapping results..." > abundances.log

#### Calculate abundances (mapped reads) from frag.gz files
ls *_kma.frag.gz | while read filename
do
    sample=${filename##*/}
    sample=${sample%%_kma.frag.gz}
    echo "Getting abundance for sample ${sample}"
    res=${sample}_kma.res
    pigz -d -p $processes $filename
    python $scriptdir/modules/contig_abundance.py \
        --frag_file ./${filename%.gz} --contig_bin_map contig_metadata_lt1000.csv --res_file $res \
        --processes $processes --out_prefix ${sample} >> abundances.log
    pigz -p $processes ${filename%.gz}
    [ ! -f ${sample}_mapping_stats.txt ] && (echo "ERROR: abundance calculation for sample ${sample} failed!" 2>&1 | tee -a abundances.log)
done

#### Get single table with all mapping stats for every samples
echo -e "sample\tbin\tmapped_reads\tidentity\tcoverage\ttotal_Length" > mapping_stats_bins.txt

ls *mapping_stats.txt | while read filename
do
    sample=${filename##*/}
    sample=${sample%%_mapping_stats.txt}
    echo "Adding mapping stats for sample ${sample}"
    awk -v s=$sample 'NR>1 {print s"\t"$0}' $filename >> mapping_stats_bins.txt
done

#### Calculate mapped reads per file
# Sum contig mapped reads per sample to get the total number of reads used to characterise the community
echo -e "Sample\tbin_reads" > bin_reads.txt

ls *mapping_stats.txt | while read filename
do
    sample=${filename##*/}
    sample=${sample%%_mapping_stats.txt}
    echo "Getting total mapped reads for sample ${sample}"
    rc=$(awk '{sum += $2} END {print sum}' $filename)
    echo -e "${sample}\t${rc}" >> bin_reads.txt
done

# Add to big read count file
python $scriptdir/modules/multiqc_to_csv.py bin_reads.txt ${outdir}/read_count.csv "bin_reads" "bin_reads"
