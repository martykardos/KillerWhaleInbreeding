#!/bin/bash
#SBATCH --job-name=makeReadBed_ASM322739v1_HiC_pruned
#SBATCH --output=makeReadBed_ASM322739v1_HiC_pruned.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript makeReadBed_ASM322739v1_HiC_pruned.R
