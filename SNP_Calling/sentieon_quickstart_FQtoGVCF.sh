#!/bin/sh

###usage##
#  sh shell outdir SampleName FQ1 FQ2 ReadGroup
###

# ******************************************
# 0. User settings 
# ******************************************

# Number of threads. You can get the number of cpu cores by running nproc  
nt=48

# Sample name and read group name. 
# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names 
# for each lane, when you need to combine results from different lanes from later on.
# in fact, the name of vairants are not meaningful enougth, but i don't mean to change again. by Shuai
name="$1"   # outdir 

# reads information
sample="$2" # sample name. if no sepcial reason, 'group' and 'sample' should be same. By Shuai
group="$5"  # sample ID. By Shuai

# Input pair-ended Illumina fastq files 
fastq_1="$3"
fastq_2="$4"


if [ $name ] || [ $sample ] || [ $fastq_1 ] || [ $fastq_2 ] || [ $group ]
then
	echo "start sentieon ..."
else
	echo "sh sentieon_quickstart_FQtoGVCF.sh [outdir] [sample name] [fq1] [fq2] [reads group id]"
	exit
fi

if [ -d $name ]
then
   echo "$name exists."
else
    mkdir -p $name
fi

# defult 
#   Sequencing platform.
pl="BGISEQ-500"
#   Fasta reference file
fasta=./ref_genome.fa
#   SNP known sites
#dbsnp=/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference/dbsnp/dbsnp-147.hg19.vcf.gz
#dbsnp=/zfsqd1/ST_OCEAN/USRS/shaolibin/database/gnomad.genomes.r2.1.hg19.vcf.gz
#   Log file 
logfile="$name/$sample.quick_start_$(date +"%Y%m%d%H%M").log"
bgzip='=./software/variation/htslib-1.4.1/bin/bgzip'
tabix='./software/variation/htslib-1.4.1/bin/tabix'

# there is no need to modify the rest of the script
# ******************************************
# 0. input locations
# ******************************************

# Test if location of the Sentieon installation directory is set
#pl="ILLUMINA" #platform

#SOURCE="$0"
set -x
exec 2>$logfile



# ******************************************
# 1. Mapping reads with BWA-MEM. Output coordinate-sorted BAM file.
# ******************************************
start_time_mod=$(date +%s)
sentieon=/share/app/st_bigdata/Sentieon/Sentieon
$sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$pl" -t $nt $fasta $fastq_1 $fastq_2 -o $name/$sample.sorted.bam
end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module BWA Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

# ******************************************
# 2. Metrics
#
# Produces Picard standard metrics of the input sequence and alignment. 
# ******************************************
bamlist="-i $name/$sample.sorted.bam"

start_time_mod=$(date +%s)
$sentieon driver -r $fasta -t $nt $bamlist --algo MeanQualityByCycle $name/$sample.mq_metrics.txt --algo QualDistribution $name/$sample.qd_metrics.txt --algo GCBias --summary $name/$sample.gc_summary.txt $name/$sample.gc_metrics.txt --algo AlignmentStat $name/$sample.aln_metrics.txt --algo InsertSizeMetricAlgo $name/$sample.is_metrics.txt 
$sentieon plot metrics -o $name/$sample.metrics_report.pdf gc=$name/$sample.gc_metrics.txt qd=$name/$sample.qd_metrics.txt mq=$name/$sample.mq_metrics.txt isize=$name/$sample.is_metrics.txt 
end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module metrics Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

# ******************************************
# 3. Remove Duplicate Reads
#
# This step will mark and remove duplicates.
# If only marking without removal is desired, remove -rmdup option in the second command.
# ******************************************
start_time_mod=$(date +%s)

$sentieon driver  -t $nt $bamlist --algo LocusCollector --fun score_info $name/$sample.score.txt
$sentieon driver  -t $nt $bamlist --algo Dedup --rmdup --score_info $name/$sample.score.txt --metrics $name/$sample.dedup_metrics.txt $name/$sample.deduped.bam 

# rm sorted bam. By Shuai
rm $name/$sample.sorted.bam $name/$sample.sorted.bam.bai

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module dedup Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

# ******************************************
# 5. Indel realigner
#
# This step produces InDel-realigned BAM file.
# ******************************************
start_time_mod=$(date +%s)

#$sentieon driver -r $fasta -t $nt -i $name/$sample.deduped.bam --algo Realigner -k $dbsnp $name/$sample.realigned.bam
$sentieon driver -r $fasta -t $nt -i $name/$sample.deduped.bam --algo Realigner $name/$sample.realigned.bam

# rm dedup bam. by Shuai
rm $name/$sample.deduped.bam $name/$sample.deduped.bam.bai

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module realign Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

# ******************************************
# 6. BQSR. Abort to do BQSR for non-model speceis By Shuai
#
# NOTES: 
# 1. -k options to known SNP sites are optional, but recommended. 
# 2. Applying the ReadWriter algo is optional, the next step in HC will apply calibration table on the fly.
# ******************************************
start_time_mod=$(date +%s)
#$sentieon driver -r $fasta  -t $nt -i $name/$sample.realigned.bam --algo QualCal $name/$sample.recal_data.table
#$sentieon driver -r $fasta  -t $nt -i $name/$sample.realigned.bam -q $name/$sample.recal_data.table --algo QualCal -k $dbsnp $name/$sample.recal_data.table.post --algo ReadWriter $name/$sample.recaled.bam
#$sentieon driver -t $nt --algo QualCal --plot --before $name/$sample.recal_data.table --after $name/$sample.recal_data.table.post $name/$sample.recal_data.csv
#$sentieon plot bqsr -o $name/$sample.bqsrreport.pdf $name/$sample.recal_data.csv

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module BQSR Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile


# ******************************************
# 7b. HC Variant caller
#
# Notes:
# 1. HC variant caller can take either Non-BQSR-calibrated BAM files and the calibration table as shown in the example below, 
#    or calibrated BAM directly. 
# 2. In the former case, where HC variant caller will apply calibration table on the fly (needs both smoke_realigned.bam and smoke_recal_data.table input)
# 3. In the latter case, use only the calibrated BAM, in this case, smoke_recaled.bam. 
# 4. DO NOT INCLUDE THE CALIBRATION TABLE TOGETHER WITH THE CALIRATED BAM INPUT, otherwise the calibration table will be applied twice.
# ******************************************
#start_time_mod=$(date +%s)
#$sentieon driver -r $fasta  -t $nt -i $name/$sample.realigned.bam -q $name/$sample.recal_data.table --algo Haplotyper -d $dbsnp $name/$sample.output-hc.vcf
#end_time_mod=$(date +%s)
#if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
#if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
#echo "Module HC Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

# ******************************************
# 7c. HC Variant caller generating GVCF
# ******************************************
start_time_mod=$(date +%s)
$sentieon driver -r $fasta  -t $nt -i $name/$sample.realigned.bam --algo Haplotyper --emit_mode gvcf $name/$sample.hc.g.vcf
#$sentieon driver -r $fasta  -t $nt -i $name/$sample.realigned.bam -q $name/$sample.recal_data.table --algo Haplotyper -d $dbsnp --emit_mode gvcf $name/$sample.output-hc.g.vcf

# bgzip. By Shuai
#$sentieon util vcfconvert $name/$sample.hc.g.vcf $name/$sample.hc.g.vcf.gz
$bgzip -@ $nt $name/$sample.hc.g.vcf
$tabix -p vcf $name/$sample.hc.g.vcf.gz

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module HC generating GVCF Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

