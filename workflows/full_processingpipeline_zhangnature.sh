#!/usr/bin/bash

##############################
# Full Processing Pipeline CRG Adel Lab for the zhang nature Dataset
# This script processes the data using the following steps:
# 1. Concatenate reads
# 2. Trim reads
# 3. Align reads
# 4. Generate a multiQC report
# 5. Run vast combine
# 6. Run vast compare
##############################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# Third job - align reads
echo "Submitting third job: Trimming and FastQ..."
jid1=$(sbatch $PWD/scripts/bash/zhang_nature/processing_trim_reads_zhang.sh | tr -cd '[:digit:].')
echo "...first job ID is $jid1"

# Third job - align reads
echo "Submitting third job: Align reads..."
jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/zhang_nature/vast_align_zhang.sh | tr -cd '[:digit:].')
echo "...third job ID is $jid2"

# Fifth job - run vast combine (dependent on third job)
echo "Submitting fifth job: Run vast combine..."
jid3=$(sbatch --dependency=afterok:$jid2 $PWD/scripts/bash/zhang_nature/vast_combine_zhang.sh | tr -cd '[:digit:].')
echo "...fifth job ID is $jid3"

# Fifth job - run vast combine (dependent on third job)
echo "Submitting fifth job: MultiQC..."
jid4=$(sbatch --dependency=afterok:$jid3 $PWD/scripts/bash/multiqc.sh | tr -cd '[:digit:].')
echo "...fifth job ID is $jid3"

echo "All jobs submitted!"