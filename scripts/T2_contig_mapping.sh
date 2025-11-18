#!/bin/bash -l

#SBATCH -p main
#SBATCH -n 100
#SBATCH -t 2-00:00:00
#SBATCH -J contig_abundances
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

# Map metagenomic reads to cluster contig reference and calculate abundances

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

subdir=${outdir}T2_contig_mapping/ # Output subdirectory
contigs=${outdir}M1_contig_assemblies/clustered_contigs_rep_seq.fasta # Input contigs
readdir=${outdir}3_mapped_seqs # Input directory with mapped sequences

[ -d $subdir ] || mkdir $subdir # Create subdirectory if not there

## Determine samples to be run
sample_list=${indir}/sample_list.csv # The file list determining the samples to be processed

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo "Taking input from: ${contigdir}"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

processes=$SLURM_NTASKS
[ -z $processes ] && processes=1
echo "Number of tasks: ${processes}"

#################
#### RUN KMA ####
#################

cd $subdir 

# Index
if [ -f database.comp.b ]; then
    echo "Index already exists, skipping indexing step"
else
    echo "Index does not exist, creating index"
    kma index -i $contigs -o database 2> index.log
fi

tail -n+2 ${sample_list} | while IFS=, read -r sample species _ _
do
    # Check if output already exists
    if [ -f ${subdir}/${sample}_kma.res ]; then
        echo "Output for sample ${sample} already exists, skipping..."
        continue
    fi
    echo "Mapping sample ${sample} (${species})"
    
    # Find mapped reads files
    read_file=$(ls ${readdir}/${sample}_unmapped.fastq.gz 2> /dev/null | head -n 1)
    
    # Run KMA to map reads to contigs and generate SAM file
    echo "Mapping..."
    kma -i $read_file -o ${subdir}/${sample}_kma \
        -t_db database -t ${processes} -nc -1t1 2>> ${sample}_mapping.log
        
    echo "Finished mapping sample ${sample}"
done

######################
#### Collect info ####
######################

# Log messages
echo "Calculating abundances from mapping results..." > abundances.log

#### Calculate abundances (mapped reads) from frag.gz files
ls *_kma.frag.gz | while read filename
do
    sample=${filename##*/}
    sample=${sample%%_kma.frag.gz}
    res=${sample}_kma.res
    pigz -d -p $processes $filename
    echo "Getting abundance for sample ${sample}"
    python $scriptdir/modules/contig_abundance.py \
        --frag_file ./${filename%.gz} --res_file $res \
        --out_prefix ${sample} --processes $processes >> abundances.log
    pigz -p $processes ${filename%.gz}
    [ ! -f ${sample}_mapping_stats.txt ] && (echo "ERROR: abundance calculation for sample ${sample} failed!" 2>&1 | tee -a abundances.log)
done

#### Calculate mapped reads per file
# Sum contig mapped reads per sample to get the total number of reads used to characterise the community
echo -e "Sample\tcontig_reads" > contig_reads.txt

ls *mapping_stats.txt | while read filename
do
    sample=${filename##*/}
    sample=${sample%%_mapping_stats.txt}
    echo "Getting total mapped reads for sample ${sample}"
    rc=$(awk '{sum += $3} END {print sum}' $filename)
    echo -e "${sample}\t${rc}" >> contig_reads.txt
done

# Add to big read count file
python $scriptdir/modules/multiqc_to_csv.py contig_reads.txt ${outdir}/read_count.csv "contig_reads" "contig_reads"
