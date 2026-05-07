#!/bin/bash -l

#### Arguments ####
# $1: Suffix of arguments for which decOM should be run
# $2: Number of tasks
suffix=$1
n_tasks=$2

echo "Running decOM"
echo "Suffix: $suffix"
echo "N tasks: $n_tasks"

#### Set up ####

# Define path to sources, which sould be downloaded before using decOM
source=${REFDIR}/decOM_sources/p_sources
echo "Sources in $source"
decom_in=./decOM_input

[ ! -d $source ] && echo "Sources not found! Download them and provide path." && exit

[ -d $decom_in ] || ( mkdir $decom_in && echo "Creating directory $decom_in" ) # Create subdirectory for the input files

# Get sinks.txt file to use with -p_sinks as well as keys files to use with -p_keys, using the sample list

[ -f $decom_in/sinks.txt ] && rm -f $decom_in/sinks.txt && echo "Removing existing sinks.txt file"

ls *${suffix} | while read input
do
  s=${input%$suffix}
  s=${s%.A} # The .A suffix seems to be throwing decOM off
  if [[ ! -f $input ]]
  then
    echo "File $input doesn't exist. Skipping"
  else
    echo "Creating input files for $s"
    echo $s >> $decom_in/sinks.txt
    echo "${s} : ${input}" > $decom_in/${s}.fof
  fi
done

#### Run decOM ####

sinks=$decom_in/sinks.txt
keys=$decom_in/
decOM-MST -p_sources $source -m $source/../map.csv -p_sinks $sinks -p_keys $keys --plot False -mem 200GB -t $n_tasks -o ./decOM_output