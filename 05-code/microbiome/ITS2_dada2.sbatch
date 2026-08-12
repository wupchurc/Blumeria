#!/bin/bash -l
# ITS2- identifying fungal species with dada2

#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu

module load R/4.4.0-openblas-rocky8
Rscript ITS_taxa.R
