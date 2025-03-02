#!/bin/bash
#SBATCH -p sapphire
#SBATCH --job-name=taxo
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem=300g
#SBATCH --tmp=50g
#SBATCH --time=55:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aguzman@fas.harvard.edu
cd /n/home11/guzman/
module load R/4.4.1-fasrc01
Rscript --vanilla /n/netscratch/jcavender_lab/Everyone/guzman/MMDB/R/Diversity/02_NMDS_taxonomic.R


