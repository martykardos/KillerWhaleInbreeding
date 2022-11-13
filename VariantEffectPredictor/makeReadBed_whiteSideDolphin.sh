#!/bin/bash
#SBATCH --job-name=makeReadBed_whiteSideDolphin
#SBATCH --output=makeReadBed_whiteSideDolphin.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript makeReadBed_whiteSideDolphin.R



