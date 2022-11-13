#!/bin/sh
#SBATCH --job-name=bwa_Ovir.te_1.0_HiC_pruned_reads1
#SBATCH --output=bwa_Ovir.te_1.0_HiC_pruned_reads1.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa Ovir.te_1.0_HiC_pruned_reads1.fastq > Ovir.te_1.0_HiC_pruned_reads1.sam
