#!/bin/bash
#SBATCH --job-name=Ovir.te_1.0_HiC_pruned_reads1_makeFastq
#SBATCH --output=Ovir.te_1.0_HiC_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript Ovir.te_1.0_HiC_pruned_reads1_makeFastq.R
