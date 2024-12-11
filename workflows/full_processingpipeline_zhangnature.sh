#!/usr/bin/bash

##############################
# Full Processing Pipeline CRG Adel Lab for the zhang nature Dataset
# This script processes the data using the following steps:
# 1. Concatenate reads
# 2. Trim reads
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
# First job - Download FastQC
echo "Submitting first job: Downloading fastq..."
jid1=$(sbatch $PWD/scripts/bash/zhang_nature/download_fastq_zhang.sh data/raw/zhang_nature/zhang_fasta_files.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid1"

# Second job - fastqc
echo "Submitting first job: FastQC..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/zhang_nature/run_fastqc_zhang.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid2"

# Second job - align reads (dependent on first job)
echo "Submitting second job: Align reads..."
jid3=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/zhang_nature/vast_align_zhang.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...second job ID is $jid3"

# Third job - run vast combine (dependent on second job)
echo "Submitting third job: Run vast combine..."
jid4=$(sbatch --dependency=afterok:$jid3 $PWD/scripts/bash/zhang_nature/vast_combine_zhang.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...third job ID is $jid4"


echo "All jobs submitted!"