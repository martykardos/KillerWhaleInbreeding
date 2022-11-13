#!/bin/bash
#SBATCH --job-name=makeReadBed_longFinPilotWhale
#SBATCH --output=makeReadBed_longFinPilotWhale.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript makeReadBed_longFinPilotWhale.R
