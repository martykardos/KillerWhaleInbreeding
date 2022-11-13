#!/bin/sh
#SBATCH --job-name=bwa_Eubalaena_glacialis_HiC_pruned_reads1
#SBATCH --output=bwa_Eubalaena_glacialis_HiC_pruned_reads1.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa Eubalaena_glacialis_HiC_pruned_reads1.fastq > Eubalaena_glacialis_HiC_pruned_reads1.sam
