#!/bin/sh
#SBATCH --job-name=bwa_Bison_UMD1.0_HiC_pruned
#SBATCH --output=bwa_Bison_UMD1.0_HiC_pruned.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa Bison_UMD1.0_HiC_pruned.fastq > Bison_UMD1.0_HiC_pruned.sam
