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
#SBATCH --array=0-9

# job name
#SBATCH --job-name trimm_fastqc

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
mkdir -p $PWD/data/processed/trapnell2010/fastqc

files=($PWD/data/raw/trapnell2010/*_1.fastq.gz)
file1=${files[$SLURM_ARRAY_TASK_ID]}
file2=${file1/_1.fastq.gz/_2.fastq.gz}

singularity run --bind $PWD/data/raw/trapnell2010:/data/raw \
                --bind $PWD/data/processed/trapnell2010:/data/processed \
                https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0 \
                trim_galore --paired --fastqc -j 4 \
                ${file1} ${file2} \
                -o /data/processed \
                --fastqc_args "--outdir /data/processed/fastqc -t 8"



###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds