#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=10

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=2G

# job name
#SBATCH --job-name downloadfasta

# job array directive
#SBATCH --array=0-23

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
mkdir -p $PWD/data/raw/tao_molbiotech
cd $PWD/data/raw/tao_molbiotech

echo "Running command from file:"
sed "$((SLURM_ARRAY_TASK_ID + 1))q;d" fastq_files_tao.sh | bash

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
