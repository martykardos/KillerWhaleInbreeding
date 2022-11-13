#!/bin/sh
#SBATCH --job-name=bwa_GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1
#SBATCH --output=bwa_GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load aligners/bwa
bwa mem -B 3 kw_ref.fa GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1.fastq > GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1.sam
