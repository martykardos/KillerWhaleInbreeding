#!/bin/sh
#SBATCH --job-name=bwa_GCF_005190385.1_NGI_Narwhal_1_genomic_pruned
#SBATCH --output=bwa_GCF_005190385.1_NGI_Narwhal_1_genomic_pruned.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1.fastq > GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1.sam
