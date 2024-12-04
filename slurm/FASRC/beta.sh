#!/bin/bash
#SBATCH -p sapphire
#SBATCH --job-name=beta
#SBATCH --nodes=1
#SBATCH --ntasks=112
#SBATCH --mem=300g
#SBATCH --tmp=100g
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aguzman@fas.harvard.edu
cd /n/home11/guzman/
module load R/4.4.1-fasrc01
Rscript --vanilla /n/netscratch/jcavender_lab/Everyone/guzman/MMDB/R/Diversity/01-beta_phylogenetic.R
Rscript --vanilla /n/netscratch/jcavender_lab/Everyone/guzman/MMDB/R/Diversity/01-beta_functional.R
Rscript --vanilla /n/netscratch/jcavender_lab/Everyone/guzman/MMDB/R/Diversity/01-beta_taxonomic.R

