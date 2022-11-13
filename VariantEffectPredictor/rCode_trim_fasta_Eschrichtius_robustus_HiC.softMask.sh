#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_Eschrichtius_robustus_HiC.softMask
#SBATCH --output=rCode_trim_fasta_Eschrichtius_robustus_HiC.softMask.log
#SBATCH --nodelist=node16
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus/EschrichtiusTrim
module load R 
Rscript rCode_trim_fasta_Eschrichtius_robustus_HiC.softMask.R
