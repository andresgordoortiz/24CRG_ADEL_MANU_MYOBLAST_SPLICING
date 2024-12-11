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
jid1=$(sbatch $PWD/scripts/bash/tao_molbiotech/download_fastqc_tao.sh /data/raw/tao_molbiotech/fastq_files_tao.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid1"

# Second job - fastqc
echo "Submitting first job: FastQC..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/tao_molbiotech/run_fastqc_tao.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid2"

# Third job - align reads
echo "Submitting second job: Align reads..."
jid3=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/tao_molbiotech/vast_align_tao.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...second job ID is $jid3"

# Fourth job - run vast combine (dependent on second job)
echo "Submitting third job: Run vast combine..."
jid4=$(sbatch --dependency=afterok:$jid3 $PWD/scripts/bash/tao_molbiotech/vast_combine_tao.sh $VASTDB_PATH | tr -cd '[:digit:].')
echo "...third job ID is $jid4"


echo "All jobs submitted!"