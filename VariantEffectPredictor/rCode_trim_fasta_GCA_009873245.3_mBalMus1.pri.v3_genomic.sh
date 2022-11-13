#!/bin/bash
#SBATCH --job-name=rCode_trim_fasta_GCA_009873245.3_mBalMus1.pri.v3_genomic
#SBATCH --output=rCode_trim_fasta_GCA_009873245.3_mBalMus1.pri.v3_genomic.log
#SBATCH --nodelist=node11
#SBATCH --mem=90GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R 
Rscript rCode_trim_fasta_GCA_009873245.3_mBalMus1.pri.v3_genomic.R
