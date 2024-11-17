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

mkdir -p $PWD/tmp2/fastqc

# First job - fastqc
echo "Submitting first job: FastQC..."
jid1=$(sbatch $PWD/scripts/bash/tao_molbiotech/run_fastqc_tao.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid1"

# Second job - align reads (dependent on first job)
echo "Submitting second job: Align reads..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/tao_molbiotech/vast_align_tao.sh | tr -cd '[:digit:].')
echo "...second job ID is $jid2"

# Third job - run vast combine (dependent on second job)
echo "Submitting third job: Run vast combine..."
jid3=$(sbatch --dependency=afterok:$jid2 $PWD/scripts/bash/tao_molbiotech/vast_combine_tao.sh | tr -cd '[:digit:].')
echo "...third job ID is $jid3"


echo "All jobs submitted!"