#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_GCF_005190385.1_NGI_Narwhal_1_genomic
#SBATCH --output=rCode_trim_fasta_GCF_005190385.1_NGI_Narwhal_1_genomic.log
#SBATCH --nodelist=node13
#SBATCH --mem=90GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R 
Rscript rCode_trim_fasta_GCF_005190385.1_NGI_Narwhal_1_genomic.R
