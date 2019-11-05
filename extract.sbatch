#! /bin/bash -eu
#SBATCH --job-name=extract
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=%j-output.log
#SBATCH --error=%j-error.log

./extract_data.sh
