#main steps for using sentieon to call SNP and Indel by population using joint calling method

step1. from fastq to gvcf
sentieon_quickstart_FQtoGVCF.sh
[eg] sh sentieon_quickstart_FQtoGVCF.sh /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/02.snp/01.fqtoGVCF//offshore1 offshore1 /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/00.reads//offshore1/offshore1_0_1_clean.fq.gz /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/00.reads//offshore1/offshore1_0_2_clean.fq.gz offshore1

step2. from gvcf to vcf [optional: you can use by chromosome]
sentieon_quickstart_GVCFtoVCF.sh
[eg] sh sentieon_quickstart_GVCFtoVCF.sh /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/02.snp/02.joint_calling/kw_151 /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/02.snp/02.joint_calling/kw_g.vcf.lst /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/02.snp/02.joint_calling/chr9

step3. [merge raw vcf, sort] split SNP and Indel and get high-quality SNP and Indel
VCF.merge.sort.split.HQ.sh
[eg] sh VCF.merge.sort.split.HQ.sh -g /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/01.genome/kw_ref.fa -v chr_vcf.lst -p kw_151 -t 10

step5. finally filter, do SNP and Indel, seperatly
VCF.final_filter.sh
[eg] sh VCF.final_filter.sh -g /dellfsqd2/P18Z19700N0073/zhangyaolei/project/killer_whale/01.genome/kw_ref.fa -v ../03.merge/kw_151.snp_indel.raw.vcf -p kw_151 -c SNP

