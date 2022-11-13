#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_ASM322739v1_HiC
#SBATCH --output=rCode_trim_fasta_ASM322739v1_HiC.log
#SBATCH --nodelist=node24
#SBATCH --mem=90GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript rCode_trim_fasta_ASM322739v1_HiC.R
