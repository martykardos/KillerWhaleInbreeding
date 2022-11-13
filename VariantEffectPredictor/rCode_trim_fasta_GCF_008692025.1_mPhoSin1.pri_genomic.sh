#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_GCF_008692025.1_mPhoSin1.pri_genomic
#SBATCH --output=rCode_trim_fasta_GCF_008692025.1_mPhoSin1.pri_genomic.log
#SBATCH --nodelist=node15
#SBATCH --mem=90GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R 
Rscript rCode_trim_fasta_GCF_008692025.1_mPhoSin1.pri_genomic.R
