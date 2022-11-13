#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_Bison_UMD1.0_HiC
#SBATCH --output=rCode_trim_fasta_Bison_UMD1.0_HiC.log
#SBATCH --nodelist=node25
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus/bisonTrim
module load R 
Rscript rCode_trim_fasta_Bison_UMD1.0_HiC.R
