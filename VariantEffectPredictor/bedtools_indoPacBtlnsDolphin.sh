#!/bin/bash
#SBATCH --job-name=bedtools_indoPacBtlnsDolphin
#SBATCH --output=bedtools_indoPacBtlnsDolphin.log
#SBATCH --mem=20GB
#SBATCH -t 20:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load bio/bedtools
bedtools getfasta -fi ASM322739v1_HiC_pruned.fasta -bed ASM322739v1_HiC_pruned_read.bed > ASM322739v1_HiC_pruned_reads1

