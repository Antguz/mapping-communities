#!/bin/bash
#SBATCH -p sapphire
#SBATCH --job-name=mantel_msr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=200g
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aguzman@fas.harvard.edu
cd /n/netscratch/jcavender_lab/Everyone/guzman/MMDB
module load R/4.4.1-fasrc01
export OMP_NUM_THREADS=100
Rscript --vanilla R/Diversity/02_mantel_msr.R
