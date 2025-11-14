#!/bin/bash -l

# This script downloads databases for MAG analysis

## Change to the directory containing the script
cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

echo "Raw data dir: ${RAWDIR}"
echo "Project dir: ${PROJDIR}"
echo "Output dir: ${OUTDIR}"
echo "Reference dir: ${REFDIR}"

conda activate gtdb2cat

if [ -z "$SLURM_NTASKS" ]; then
  processes=2
else
  processes=$(echo "$SLURM_NTASKS * 2" | bc)
fi

cd $REFDIR

############################
#### DOWNLOAD DATABASES ####
############################

## CheckM: download and unpack
mkdir checkm_data_2015_01_16
cd checkm_data_2015_01_16
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
pigz -d -p $processes checkm_data_2015_01_16.tar.gz
tar -xf checkm_data_2015_01_16.tar
rm checkm_data_2015_01_16.tar
cd ..

## GTDB
mkdir gtdbtk_r220_data
cd gtdbtk_r220_data
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
pigz -d -p $processes gtdbtk_r220_data.tar.gz
tar -xf gtdbtk_r220_data.tar

## BUSCO
wget https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz
gunzip bacteria_odb10.2024-01-08.tar.gz
tar -xf bacteria_odb10.2024-01-08.tar
rm bacteria_odb10.2024-01-08.tar

## Create CAT database using GTDB
gtdb=$REFDIR/gtdbtk_r220_data/release220/ # GTDB
catdb=$REFDIR/CAT_DB # Where the CAT db will be created

[ -d $catdb ] || mkdir $catdb 

cd $catdb

# Using GTDB database for this
cp $gtdb/taxonomy/* $catdb

# Get names.dmp and nodes.dmp using jupyter notebook adapted from Nikos Pappas and the MGXlab https://github.com/MGXlab/gtdb2cat
python $scriptdir/modules/gtdbtk2cat.py

CAT_pack prepare \
        --db_fasta $catdb/db/gtdb.nr.fa \
        --names $catdb/taxonomy/names.dmp \
        --nodes $catdb/taxonomy/nodes.dmp \
        --acc2tax $catdb/db/prot.accession2taxid.txt \
        --db_dir $catdb/cat_db \
        --verbose
