#!/bin/sh
#SBATCH --job-name=samtools_cetaceansVariantCall
#SBATCH --output=samtools_cetaceansVariantCall.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load bio/samtools
samtools mpileup -g --positions kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1.hetDepthFiltered.maxMiss0.25_positions_noHeader -f kw_ref.fa ASM367639v1_HiC_pruned_sorted.bam Eschrichtius_robustus_HiC.softMask_pruned_reads1_sorted.bam GCF_005190385.1_NGI_Narwhal_1_genomic_pruned_reads1_sorted.bam Ovir.te_1.0_HiC_pruned_reads1_sorted.bam ASM654740v1_HiC_pruned_sorted.bam Eubalaena_glacialis_HiC_pruned_reads1_sorted.bam GCF_008692025.1_mPhoSin1.pri_genomic_pruned_reads1_sorted.bam Bison_UMD1.0_HiC_pruned_sorted.bam GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1_sorted.bam  Mesoplodon_europaeus_HiC.softMask_pruned_reads1_sorted.bam ASM322739v1_HiC_pruned_sorted.bam  > cetaceans.raw.bcf

