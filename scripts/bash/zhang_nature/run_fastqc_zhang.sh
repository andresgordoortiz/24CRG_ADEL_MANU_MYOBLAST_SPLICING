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
#SBATCH --cpus-per-task=8

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
mkdir -p $PWD/data/processed/zhang_nature/fastqc
# Run FastQC using Singularity
singularity exec --bind $PWD/data/raw/zhang_nature \
    docker://biocontainers/fastqc:v0.11.9_cv8 \
    fastqc -t 8 $PWD/data/raw/zhang_nature/*.fastq.gz \
    --outdir $PWD/data/processed/zhang_nature/fastqc

################
# run multiqc  #
################
singularity exec --bind $PWD/data/processed/zhang_nature:/zhang_nature \
    docker://multiqc/multiqc:latest \
    /bin/bash -c "cd /zhang_nature && multiqc . -n zhang_nature_molbiotech_multiqc_report.html"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds