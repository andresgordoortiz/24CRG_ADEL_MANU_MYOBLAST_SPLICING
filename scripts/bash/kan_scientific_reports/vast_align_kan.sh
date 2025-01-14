#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=240

# queue
#SBATCH --qos=short
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name vast-align
# job array directive
#SBATCH --array=4

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail


###############
# run command #
###############

#Define file list and select the file for the current array job
files=($PWD/data/raw/kan_scientific_reports/*_1.fastq.gz)
file1=${files[$SLURM_ARRAY_TASK_ID]}
file2=${file1/_1.fastq.gz/_2.fastq.gz}

basename=$(basename "$file1" _1.fastq.gz)
mkdir -p $PWD/data/processed/kan_scientific_reports/vast_out


singularity_image="docker://andresgordoortiz/vast-tools:latest"
VASTDB_PATH=$1
# Run vast-tools align using Singularity
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $PWD/data/processed/kan_scientific_reports/vast_out:/vast_out \
    $singularity_image vast-tools align \
    "$file1" "$file2" \
    -sp mm10 \
    -o /vast_out \
    --IR_version 2 \
    -c 8 \
    -n "$basename"



###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
