#!/bin/bash

#SBATCH --job-name=LVIScomp
#SBATCH --time=04:00:00
#SBATCH --mem=3G
#SBATCH --chdir=/scratch/lkp58/LVIS/accuracy_assessment/
#SBATCH --output=./logs/%A

### ---- Inputs ----
# "/projects/above_gedi/users/lpuckett/ShovelCreek/LVISF2_ABoVE_2019/1a_points/WGS84/"

distance_threshold=1 # meters
input_dir1="/projects/above_gedi/users/lpuckett/Rainbow2/LVIS_2017/1_points/WGS84/"
input_dir2="/projects/above_gedi/users/lpuckett/Rainbow2/LVISF2_ABoVE_2019/1_points/WGS84/"
filename1="LVIS2_ABoVE2017_0703_R1803_087324.shp"
filename2="LVISF2_ABoVE2019_0801_R2003_087625.shp"

#start=$(date +%s)
#date_time_inner=`date +%Y%m%d_%H%M%S`
#echo "The starting date_time: "\$date_time_inner
#echo
#echo "SLURM_JOBID: " $SLURM_JOBID
#echo "SLURM ARRAY TASK ID: " $SLURM_ARRAY_TASK_ID
#echo


# Run the R script
module load R/3.6.2 proj/6.1.0 gdal/3.0.4
Rscript ./code/join_points.R $distance_threshold $input_dir1 $input_dir2 $filename1 $filename2

# - ENDING -
#echo "Ended at:"
#date