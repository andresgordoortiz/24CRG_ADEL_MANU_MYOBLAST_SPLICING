#!/usr/bin/bash

##############################
# Full Processing Pipeline CRG Adel Manu Lab for Pladienolide B dataset
# This script processes the Pladienolide B dataset (unpublished) by performing the following steps:
# 1. Concatenate reads
# 2. Trim reads
# 3. Align reads
# 4. Generate a multiQC report
# 5. Run vast combine
# 6. Run vast compare
##############################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# Third job - align reads
echo "Submitting third job: Align reads..."
jid1=$(sbatch $PWD/scripts/bash/vast_align_trapnell.sh | tr -cd '[:digit:].')
echo "...third job ID is $jid1"

# Fifth job - run vast combine (dependent on third job)
echo "Submitting fifth job: Run vast combine..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/vast_combine_trapnell.sh | tr -cd '[:digit:].')
echo "...fifth job ID is $jid2"

echo "All jobs submitted!"