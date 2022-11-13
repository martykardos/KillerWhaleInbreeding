#!/bin/bash
#SBATCH --job-name=Eubalaena_glacialis_HiC_pruned_reads1_makeFastq
#SBATCH --output=Eubalaena_glacialis_HiC_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript Eubalaena_glacialis_HiC_pruned_reads1_makeFastq.R
