#!/bin/sh

###usage##
#  sh shell out_prefix GVCFList [Chr list] 
###

# ******************************************
# 0. User settings 
# ******************************************

prefix="$1"    # out vcf file
in_gvcf="$2"    # 
chr_list="$3"   #

if [ $prefix ] || [ $in_gvcf ]
then
	echo "start sentieon ..."
else
	echo "sh sentieon_quickstart_GVCFtoVCF.sh <ouput prefix> <input gvcf list> [chr list file]"
	exit
fi

if [ ! -f $in_gvcf ]; then
    echo "GVCF not exist"
    exit 1
fi

if [ -f $chr_list ]; then
    chr_list=`readlink -f $chr_list`
    interval="--interval $chr_list"
fi

# deal with path
opt_dirname=`cd \`dirname $prefix\`;pwd`
opt_basname=`basename $prefix`
prefix="$opt_dirname/$opt_basname"
chr_prefix=${chr_list##*/}
chr_prefix=${chr_prefix%.*}

in_gvcf=`readlink -f $in_gvcf`

# user defult 
#   Number of threads. You can get the number of cpu cores by running nproc  
nt=48
#   Fasta reference file
fasta=./ref_genome.fa

# default set
sentieon=/share/app/st_bigdata/Sentieon/Sentieon
bgzip=./software/variation/htslib-1.4.1/bin/bgzip
tabix=./software/variation/htslib-1.4.1/bin/tabix
logfile="$prefix.$chr_prefix.GVCFtoVCF_$(date +"%Y%m%d%H%M").log"

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
# 1. Joint calling
#
# The GVCFtyper algorithm performs the joint variant calling of multiple samples
# Each single sample must be previously processed using the Haplotyper algorithm with the option --emit_mode gvcf
# ******************************************
start_time_mod=$(date +%s)

# 1.1 GVCFtyper
gvcf_argument=""

while read -r line; do
    gvcf_argument="$gvcf_argument -v $line"
done < $in_gvcf

$sentieon driver -r $fasta $interval --algo GVCFtyper $gvcf_argument ${prefix}.${chr_prefix}.raw.vcf

# 1.2 bgzip
$bgzip -@ $nt ${prefix}.$chr_prefix.raw.vcf
$tabix -p vcf ${prefix}.$chr_prefix.raw.vcf.gz

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module GVCFtyper Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

