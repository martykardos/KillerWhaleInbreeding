#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_Mesoplodon_europaeus_HiC.softMask
#SBATCH --output=rCode_trim_fasta_Mesoplodon_europaeus_HiC.softMask.log
#SBATCH --nodelist=node15
#SBATCH --mem=90GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript rCode_trim_fasta_Mesoplodon_europaeus_HiC.softMask.R
