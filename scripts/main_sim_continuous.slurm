#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=logs/main_sim_continuous.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1000mb
#SBATCH --time=96:00:00
#SBATCH --qos=amolstad
#SBATCH --account=amolstad

export OMP_NUM_THREADS=1

ml R/4.1
Rscript scripts/main_sim_continuous.R

