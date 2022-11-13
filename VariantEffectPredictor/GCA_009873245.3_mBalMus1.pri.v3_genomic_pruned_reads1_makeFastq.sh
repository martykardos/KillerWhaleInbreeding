#!/bin/bash
#SBATCH --job-name=GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1_makeFastq
#SBATCH --output=GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1_makeFastq.R
