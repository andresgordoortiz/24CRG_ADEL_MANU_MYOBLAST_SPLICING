#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=3

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=2G

# job name
#SBATCH --job-name downloadfasta

# job array directive
#SBATCH --array=0-11

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
mkdir -p $PWD/data/raw/zhang_nature
cd $PWD/data/raw/zhang_nature

# Check if there are less than 24 fastq.gz files
file_count=$(ls -1 *.fastq.gz 2>/dev/null | wc -l)
if [ "$file_count" -lt 12 ]; then
    sed "$((SLURM_ARRAY_TASK_ID + 1))q;d" "$PWD/$1" | bash
else
    echo "There are already 12 or more fastq.gz files in the directory."
fi

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
