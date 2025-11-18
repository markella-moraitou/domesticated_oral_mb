#!/bin/bash -l

# Get tables for analysis of contigs combining CAT taxonomies, depths (from nfcore/mag) and DRAM annotations

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

indir=${PROJDIR}/input/ # Directory where the provided input is located
scriptdir=${PROJDIR}/scripts/ # Directory where the scripts are saved

outdir=${OUTDIR} # Output subdirectory for the output of this specific script

subdir=${outdir}/Final_tables/ # Output subdirectory
catdir=${outdir}/T1_contig_taxonomy/ # Input CAT taxonomy
annotdir=${outdir}/F1_contig_annotations/raw # Input DRAM annotations
depthdir=${outdir}/T2_contig_mapping # Input read depths
damagedir=${outdir}/M2_nfcore_mag/Ancient_DNA/pydamage # Input pydamage analysis

mkdir -p $subdir # Create subdirectory if not there

## Load modules and activate environment
conda activate oral_mb_evol

## Print some info
echo -e "Taking input from:\n${catdir}\n${annotdir}\n${depthdir}/"
echo "Saving output to: ${subdir}"

echo "Global modules:"
module list
echo "Conda modules:"
conda list

####################
#### GET TABLES ####
####################

cd $subdir 

#### ABUNDANCE TABLE ####

python $scriptdir/modules/prep_CAT_abundance_table.py \
  --mapping_dir $depthdir \
  --mapping_dir $depthdir \
  --taxonomy_dir $catdir \
  --mapping_minid 98 \
  --mapping_mincov 0 \
  --tax_ranks "species,genus,family" > prep_CAT_abundance_table.log

#### TAXONOMY AND GENE TABLE ####

python $scriptdir/modules/combine_contig_func_tax.py \
      --dram_annotations $annotdir/annotations.tsv \
      --cat_annotations $catdir \
      --mapping_dir $depthdir \
      --mapping_minid 98 \
      --mapping_mincov 50 \
      --tax_ranks "species,genus,family" > combined_contig_annotations.log

# Compress
gzip combined_contig_annotations.tsv
gzip gene_abundance_stratified.tsv

#### PYDAMAGE SUMMARY ####

python $scriptdir/modules/summarise_pydamage_res.py \
            --damage_dir $damagedir/analyze \
            --taxonomy_dir $catdir \
            --groupvar "taxon" # Summarise only by taxon

mv pydamage_summary.tsv pydamage_summary_CAT.tsv 