#!/bin/bash
#SBATCH -p msismall
#SBATCH --job-name=beta-EMIT
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=750g
#SBATCH --tmp=390g
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guzman@umn.edu
cd /home/cavender/guzman/Imagery
module load gdal/3.7.0-gcc-8.2.0-tm2evb5
module load geos/3.11.2-gcc-8.2.0-km5gytt
module load proj/7.2.1-gcc-8.2.0-6kfkugi
module load R/4.2.2-gcc-8.2.0-vp7tyde
Rscript --vanilla /panfs/jay/groups/17/cavender/guzman/MMDB/R/Models/07_apply-model_EMIT.R


