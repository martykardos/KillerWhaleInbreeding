#!/bin/bash
#SBATCH --job-name=Mesoplodon_europaeus_HiC.softMask_pruned_reads1_makeFastq
#SBATCH --output=Mesoplodon_europaeus_HiC.softMask_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
R CMD BATCH Mesoplodon_europaeus_HiC.softMask_pruned_reads1_makeFastq.R
