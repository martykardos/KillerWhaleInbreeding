#!/bin/bash
#SBATCH --job-name=GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1_makeFastq
#SBATCH --output=GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1_makeFastq.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1_makeFastq.R
