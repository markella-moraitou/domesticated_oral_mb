#!/bin/bash -l

#SBATCH -p long
#SBATCH -n 50
#SBATCH -t 3-00:00:00
#SBATCH -J mcmcglmm
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/job-%x.%j.out
#SBATCH --error=logs/job-%x.%j.err

set -x # Print each command before executing
set -e # Exit on any error

## Monitor memory usage
echo "Memory info at start:"
free -h
cat /proc/meminfo | grep MemTotal

## Change to the directory containing the script
[ ! -z $SLURM_SUBMIT_DIR ] && cd $SLURM_SUBMIT_DIR

## Define paths
source ../config.txt # Should contain the paths for the raw data, the intermediate output, the project directory, the references and the conda env

## Load modules and activate environment
module load bioinfo-tools
conda activate R_env
echo "Activated R_env"

cd /cfs/klemming/projects/snic/sllstore2017021/MARKELLA/mammal_om_evol_stats/scripts/community_analysis/
echo "Changed working dir"
echo $PWD
Rscript diff_abund_tests.R

#cd /cfs/klemming/projects/snic/sllstore2017021/MARKELLA/mammal_om_evol_stats/scripts/function/
#echo "Changed working dir"
#echo $PWD
#Rscript statistical_tests.R

echo "Memory info at end:"
free -h
