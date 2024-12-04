#!/bin/bash
#SBATCH -p agsmall
#SBATCH --job-name=NMDS-functional
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=499g
#SBATCH --tmp=500g
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guzman@umn.edu
cd /home/cavender/guzman/MMDB
module load R/4.4.0-openblas-rocky8mod
Rscript --vanilla /home/cavender/guzman/MMDB/R/04_NMDS_functional.R

