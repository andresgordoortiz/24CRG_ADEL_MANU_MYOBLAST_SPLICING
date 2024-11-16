#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/tmp/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/tmp/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=1

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=500

# job name
#SBATCH --job-name vast-compare

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
group_a="SRR037946,SRR037945"
group_b="SRR037948,SRR037947,SRR037949"
group_c="SRR037950,SRR037951"
group_d="SRR037952,SRR037953,SRR037954"
name_a="24 hours"
name_b="60 hours"
name_c="120 hours"
name_d="168 hours"


# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate vasttools

cd $PWD/data/processed/pladienolideb/vast_out
/users/mirimia/projects/vast-tools/vast-tools compare INCLUSION_LEVELS_FULL-mm10-9.tab \
    -a $group_a \
    -b $group_b \
    --min_dPSI 0 \
    --min_range 0 \
    --GO --print_dPSI --print_sets ---print_AS_ev \
    -name_A $name_a  -name_B $name_b \
    -sp mm10 > summary_stats_control_vs_pladb10mM.txt

conda deactivate



###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
