#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_Eubalaena_glacialis_HiC.softMask
#SBATCH --output=rCode_trim_fasta_Eubalaena_glacialis_HiC.softMask.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus/EubalaenaTrim
module load R
Rscript rCode_trim_fasta_Eubalaena_glacialis_HiC.softMask.R
