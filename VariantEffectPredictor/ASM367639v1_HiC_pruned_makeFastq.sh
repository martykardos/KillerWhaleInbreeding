#!/bin/bash
#SBATCH --job-name=ASM367639v1_HiC_pruned_makeFastq
#SBATCH --output=ASM367639v1_HiC_pruned_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript ASM367639v1_HiC_pruned_makeFastq.R
