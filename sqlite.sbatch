#! /bin/bash -eu
#SBATCH --job-name=sqlite
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --output=%j-output.log
#SBATCH --error=%j-error.log

./import_sqlite.py
