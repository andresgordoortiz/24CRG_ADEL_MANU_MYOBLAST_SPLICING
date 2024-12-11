#!/usr/bin/bash

##############################
# Full Processing Pipeline CRG Adel Lab for Trapnell 2010 dataset
# This script processes the data by performing the following steps:
# 1. Concatenate reads
# 3. Align reads
# 4. Generate a multiQC report
# 5. Run vast combine
##############################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

if [ -z "$1" ]; then
    echo "Error: No VASTDB_PATH provided."
    echo "Usage: $0 /path/to/vastdb"
    exit 1
fi

VASTDB_PATH=$1
# First job - downlaod files
echo "Submitting first job: download files..."
jid1=$(sbatch $PWD/scripts/bash/trapnell/download_fastq_trapnell.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid1"

# Second job - run fastqc and multiqc
echo "Submitting second job: fastqc..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/trapnell/run_fastqc_trapnell.sh | tr -cd '[:digit:].')
echo "...second job ID is $jid2"

# Third job - vast align
echo "Submitting third job: vast align..."
jid3=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/trapnell/vast_align_trapnell.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...third job ID is $jid3"

# Fourth job - vast combine
echo "Submitting fourth job: vast combine..."
jid4=$(sbatch --dependency=afterok:$jid3 $PWD/scripts/bash/trapnell/vast_combine_trapnell.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...third job ID is $jid4"

echo "All jobs submitted!"