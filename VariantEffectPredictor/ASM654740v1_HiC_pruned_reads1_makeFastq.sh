#!/bin/bash
#SBATCH --job-name=ASM654740v1_HiC_pruned_reads1_makeFastq
#SBATCH --output=ASM654740v1_HiC_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript ASM654740v1_HiC_pruned_reads1_makeFastq.R
