#!/bin/bash
#SBATCH --array=1-80
#SBATCH --time=6-00:00:00
#SBATCH --mem=1000
#SBATCH --output=slurm_logs/slurm-%A_%a.out
#SBATCH --error=slurm_logs/slurm-%A_%a.err

# Load the R module
module load R/4.2.1-foss-2022a

# Set the R library path
export R_LIBS_USER=/home4/p300503/R/x86_64-pc-linux-gnu-library/4.4

# Set the replicate number based on the SLURM array index
replicate_number=$SLURM_ARRAY_TASK_ID

# RUN DAISIE on array
Rscript --vanilla DAISIE_models_constant.R ${replicate_number}