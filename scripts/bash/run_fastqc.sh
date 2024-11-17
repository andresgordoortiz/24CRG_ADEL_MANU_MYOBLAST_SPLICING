#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=60

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=5G
#SBATCH --cpus-per-task=4

# job name
#SBATCH --job-name fastqc_multiqc

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


################
# run fastqc   #
################
singularity exec --bind $PWD/data/processed/longreads_elisabeth docker://biocontainers/fastqc:v0.11.9_cv8 fastqc -t 4 $PWD/data/processed/longreads_elisabeth/*.fastq.gz

################
# run multiqc  #
################
cd $PWD/data/processed/longreads_elisabeth
module load MultiQC/1.22.3-foss-2023b
multiqc .

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds