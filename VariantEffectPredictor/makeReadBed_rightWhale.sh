#!/bin/bash
#SBATCH --job-name=makeReadBed_rightWhale.R
#SBATCH --output=makeReadBed_rightWhale.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript makeReadBed_rightWhale.R
