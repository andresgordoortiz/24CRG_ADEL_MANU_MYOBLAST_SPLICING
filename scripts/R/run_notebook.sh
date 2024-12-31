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
URL3="https://vastdb.crg.eu/downloads/mm10/PROT_IMPACT-mm10-v3.tab.gz"

FILE1="$PWD/notebooks/final/SPLICE_SITE_SCORES-mm10.tab.gz"
FILE2="$PWD/notebooks/final/EVENT_INFO-mm10.tab.gz"
FILE3="$PWD/notebooks/final/PROT_IMPACT-mm10-v2.3.tab.gz"

UNZIPPED_FILE1="${FILE1%.gz}"
UNZIPPED_FILE2="${FILE2%.gz}"
UNZIPPED_FILE3="${FILE3%.gz}"

# Create the directory if it doesn't exist
mkdir -p "$PWD/notebooks/final"

# Check and download the first file if the unzipped version doesn't exist
if [ ! -f "$UNZIPPED_FILE1" ]; then
    if [ ! -f "$FILE1" ]; then
        echo "$FILE1 not found. Downloading..."
        wget "$URL1" -O "$FILE1"
    else
        echo "$FILE1 already exists. Skipping download."
    fi
    echo "Unzipping $FILE1..."
    gunzip -c "$FILE1" > "$UNZIPPED_FILE1"
else
    echo "$UNZIPPED_FILE1 already exists. Skipping download and unzip."
fi


# Check and download the second file if the unzipped version doesn't exist
if [ ! -f "$UNZIPPED_FILE2" ]; then
    if [ ! -f "$FILE2" ]; then
        echo "$FILE2 not found. Downloading..."
        wget "$URL2" -O "$FILE2"
    else
        echo "$FILE2 already exists. Skipping download."
    fi
    echo "Unzipping $FILE2..."
    gunzip -c "$FILE2" > "$UNZIPPED_FILE2"
else
    echo "$UNZIPPED_FILE2 already exists. Skipping download and unzip."
fi

if [ ! -f "$UNZIPPED_FILE3" ]; then
    if [ ! -f "$FILE3" ]; then
        echo "$FILE3 not found. Downloading..."
        wget "$URL3" -O "$FILE3"
    else
        echo "$FILE3 already exists. Skipping download."
    fi
    echo "Unzipping $FILE3..."
    gunzip -c "$FILE3" > "$UNZIPPED_FILE3"
else
    echo "$UNZIPPED_FILE3 already exists. Skipping download and unzip."
fi
# Set the working directory inside the container to /workspace
singularity run --bind "$(pwd)/notebooks/final:/shared" \
  docker://andresgordoortiz/splicing_analysis_r_crg:v1.3 \
  bash -c "cd /; Rscript -e \"rmarkdown::render('/shared/myoblast_transcript_analysis.rmd')\""

mkdir -p $PWD/results/tables
mv $PWD/notebooks/final/results/tables $PWD/results
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
