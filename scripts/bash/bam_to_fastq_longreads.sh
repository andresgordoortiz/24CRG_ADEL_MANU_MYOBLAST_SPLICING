#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=30

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

# job name
#SBATCH --job-name bam_to_fastq
# job array directive
#SBATCH --array=0-3

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
module load SAMtools/1.20-GCC-13.2.0

# Define file list and select the file for the current array job
files=($PWD/data/processed/longreads_elisabeth/*.bam)
bam_file=${files[$SLURM_ARRAY_TASK_ID]}
fastq_file="$PWD/data/processed/longreads_elisabeth/$(basename "${bam_file}" .bam).fastq.gz"

# Convert BAM to FASTQ and compress the output
samtools bam2fq "${bam_file}" | gzip > "${fastq_file}"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

