dir=$(pwd)

#########
plink=XX/software/variation/plink_1.9/plink
snpDir=XX/14.151-killer-whales/00.data
#geno maf
$plink --noweb --file $snpDir/kw_151.snp --geno 0.05 --maf 0.005 --make-bed --out kw_151.snp_QC
#LD
$plink --bfile kw_151.snp_QC --indep-pairwise 50 5 0.2 --out pruned
$plink --bfile kw_151.snp_QC --extract pruned.prune.in --make-bed --out kw_151.snp_QC_ld
#NJtree
$plink --bfile kw_151.snp_QC_ld --distance 1-ibs --out kw_151.snp_QC_ld_nj --memory 4000 --threads 4
