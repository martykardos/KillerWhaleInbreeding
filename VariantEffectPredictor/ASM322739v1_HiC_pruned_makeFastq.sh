#!/bin/bash
#SBATCH --job-name=ASM322739v1_HiC_pruned_makeFastq
#SBATCH --output=ASM322739v1_HiC_pruned_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript ASM322739v1_HiC_pruned_makeFastq.R
