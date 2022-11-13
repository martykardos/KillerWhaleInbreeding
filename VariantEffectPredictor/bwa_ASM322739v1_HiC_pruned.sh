#!/bin/sh
#SBATCH --job-name=bwa_ASM322739v1_HiC_pruned
#SBATCH --output=bwa_ASM322739v1_HiC_pruned.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa ASM322739v1_HiC_pruned.fastq > ASM322739v1_HiC_pruned.sam
