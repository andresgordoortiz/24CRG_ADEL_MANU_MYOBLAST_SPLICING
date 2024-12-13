#!/usr/bin/bash


##################
# slurm setting #
##################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_MYOBLAST_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=120

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=15G
#SBATCH --cpus-per-task=2

# job name
#SBATCH --job-name run_RMarkdown

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

# Define URLs and filenames
URL1="https://vastdb.crg.eu/downloads/mm10/SPLICE_SITE_SCORES-mm10.tab.gz"
URL2="https://vastdb.crg.eu/downloads/mm10/EVENT_INFO-mm10.tab.gz"

FILE1="$PWD/notebooks/final/SPLICE_SITE_SCORES-mm10.tab.gz"
FILE2="$PWD/notebooks/final/EVENT_INFO-mm10.tab.gz"

# Create the directory if it doesn't exist
mkdir -p "$PWD/notebooks/final"

# Check and download the first file
if [ ! -f "$FILE1" ]; then
    echo "$FILE1 not found. Downloading..."
    wget "$URL1" -O "$FILE1"
else
    echo "$FILE1 already exists. Skipping download."
fi

# Check and download the second file
if [ ! -f "$FILE2" ]; then
    echo "$FILE2 not found. Downloading..."
    wget "$URL2" -O "$FILE2"
else
    echo "$FILE2 already exists. Skipping download."
fi
# Unzip the downloaded files
echo "Unzipping files..."
gunzip -c "$FILE1" > "${FILE1%.gz}"
gunzip -c "$FILE2" > "${FILE2%.gz}"
# Set the working directory inside the container to /workspace
singularity run --bind "$(pwd)/notebooks/final:/shared" \
  docker://andresgordoortiz/splicing_analysis_r_crg:v1.3 \
  bash -c "cd /; Rscript -e \"rmarkdown::render('/shared/myoblast_transcript_analysis.rmd')\""

mkdir -p $PWD/results/tables
mv $PWD/notebooks/final/*.csv $PWD/results/tables
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
