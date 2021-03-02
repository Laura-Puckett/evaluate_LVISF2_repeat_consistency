#!/bin/bash

#SBATCH --job-name=LVIScomp
#SBATCH --time=00:40:00
#SBATCH --mem=2G
#SBATCH --chdir=/scratch/lkp58/LVIS/accuracy_assessment/
#SBATCH --output=./logs/%A_%a
#SBATCH --array=12-14 # 1-19

### ---- Inputs ----
inputs_path='./inputs.csv'

#start=$(date +%s)
#date_time_inner=`date +%Y%m%d_%H%M%S`
#echo "The starting date_time: "\$date_time_inner
#echo
#echo "SLURM_JOBID: " $SLURM_JOBID
#echo "SLURM ARRAY TASK ID: " $SLURM_ARRAY_TASK_ID
#echo


# Run the R script
module load R/3.6.2 proj/6.1.0 gdal/3.0.4
Rscript ./code/join_points_from_txt_filelist.R  $inputs_path $SLURM_ARRAY_TASK_ID

# - ENDING -
#echo "Ended at:"
#date