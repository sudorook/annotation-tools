#! /bin/bash -eu
#SBATCH --job-name=sqlite
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=%j-output.log
#SBATCH --error=%j-error.log

./import_sqlite.py
