dir=$(pwd)

#########
plink=XX/software/variation/plink_1.9/plink
admixture=XX/software/tools/admixture/admixture_linux-1.3.0/admixture
snpDir=XX/14.151-killer-whales/00.data
#geno maf
$plink --noweb --file $snpDir/kw_151.snp --geno 0.05 --maf 0.005 --make-bed --out kw_151.snp_QC
#LD
$plink --bfile kw_151.snp_QC --indep-pairwise 50 5 0.2 --out pruned
$plink --bfile kw_151.snp_QC --extract pruned.prune.in --make-bed --out kw_151.snp_QC_ld
#admixture
for k in 1 2 3 4 5 6 7 8 9 10;
	do echo "$admixture --cv -j20 -s time  -B5 $dir/kw_151.snp_QC_ld.bed $k | tee log${k}.out" > admix.$k.sh
done

