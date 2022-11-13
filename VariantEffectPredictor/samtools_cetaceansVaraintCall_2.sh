#!/bin/sh
#SBATCH --job-name=samtools_cetaceansVariantCall_2
#SBATCH --output=samtools_cetaceansVariantCall_2.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load bio/bcftools
bcftools call -m -O v cetaceans.raw.bcf -o cetaceans.var.vcf

