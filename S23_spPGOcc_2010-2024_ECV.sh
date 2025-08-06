#!/bin/bash --login
#SBATCH --job-name=S23_spPGOcc_2010-2024_ECV
#SBATCH --cpus-per-task=10
#SBATCH --mem=256G
#SBATCH --time=130:00:00
#SBATCH --qos=normal
#SBATCH --partition=general
#SBATCH --account=a_reside
#SBATCH --constraint=epyc3
#SBATCH --batch=epyc3
#SBATCH --array=0
#SBATCH --output=log_S23_spPGOcc_2010-2024_ECV_%A_%a.out
#SBATCH --error=log_S23_spPGOcc_2010-2024_ECV_%A_%a.err

----------------------------------------------------
# Load the right module (here R 4.4.0)
module load r/4.2.1-foss-2022a



---------------------------------------------------
# Execution
srun Rscript scripts/S23_spPGOcc_2010-2024_ECV.R $SLURM_ARRAY_TASK_ID

