#!/bin/bash
#SBATCH -p agsmall
#SBATCH --job-name=beta
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=499g
#SBATCH --tmp=390g
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guzman@umn.edu
cd /home/cavender/guzman/MMDB
module load R/4.2.2-gcc-8.2.0-vp7tyde
Rscript --vanilla /home/cavender/guzman/MMDB/R/Diversity/03-beta_phylogenetic.R
Rscript --vanilla /home/cavender/guzman/MMDB/R/Diversity/03-beta_functional.R
Rscript --vanilla /home/cavender/guzman/MMDB/R/Diversity/03-beta_taxonomic.R


