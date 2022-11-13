#!/bin/bash
#=============================================================================
# Filename: VCF.merge.sort.split.HQ.sh
#=============================================================================
Bin=`cd \`dirname $0\`;pwd`

# help info
help(){
    cat <<HELP
    VCF.merge.sort.split.HQ.sh -- merge VCF from different chromosome, sort by reference, sperate SNP and Indel, fetch HQ SNP and Indel

    [usage]
     VCF.merge.sort.split.HQ.sh -g|--genome <genome.fa> -v|--vcf <in.vcf.list>
                
    [option]
     -p|--prefix   [str]   out prefix [./out]
     -t|--threads  [int]   threads to use [1]
    -h|--help             help info
                        
    [example]

    [note]
HELP
}

# pass parameter
[ -z $1 ] && { help; exit 1;}
while [ -n $1 ];do
    case $1 in
        -h|--help)      help; exit 1;;
        -g|--genome)    opt_genome=$2;  shift 2;;
        -v|--vcf)       opt_vcf=$2;     shift 2;;
        -p|--prefix)    opt_prefix=$2;  shift 2;;
        -t|--threads)   opt_threads=$2; shift 2;;
        -*)             echo "*** no option $1, -h for help"; exit 1;;
        *)              break;; # not to deal
    esac
done

# parameter check and initialize
[ -z "$opt_prefix" ] && opt_prefix='./out'
[ -z "$opt_threads" ] && opt_threads=1

# deal with path
opt_genome=`readlink -f $opt_genome`
opt_dict=${opt_genome%.*}
opt_dict="${opt_dict}.dict"

opt_dirname=`cd \`dirname $opt_prefix\`;pwd`
opt_basname=`basename $opt_prefix`
opt_prefix="$opt_dirname/$opt_basname"

opt_vcf=`readlink -f $opt_vcf`

# default
java8=./java1.8/bin/java
picard=/dellfsqd2/ST_OCEAN/USER/sunshuai/software/variation/picard/picard.jar
gatk3=./software/variation/GATK/GenomeAnalysisTK-3.8-1.jar
bgzip=./software/variation/htslib-1.4.1/bin/bgzip
tabix=./software/variation/htslib-1.4.1/bin/tabix

mem=10  # memory 10G

set -x
# step 1. merge vcf
# firstly, check vcf header same or not

error_info=`perl $Bin/VCF_header.compare.pl $opt_vcf | grep ERROR`

if [ -n "$error_info" ]; then
    echo "#ERROR: the header section of vcf files are not same"
    exit
fi


# secondly, start to merge vcf
## output header section
for i in `cat $opt_vcf`; do 
    suffix=${i##*.}
    if [ $suffix = "gz" ]; then
        header_num=`gzip -dc $i | grep -nv '^#' | head -n 1 | cut -d ':' -f 1`
        header_num=`expr $header_num - 1`
        gzip -dc $i | head -n $header_num > ${opt_prefix}.snp_indel.tmp.vcf
    else
        header_num=`grep -nv '^#' $i | head -n 1 | cut -d ':' -f 1`
        header_num=`expr $header_num - 1`
        head -n $header_num $i > ${opt_prefix}.snp_indel.tmp.vcf
    fi
done
## output body section
for i in `cat $opt_vcf`; do 
    suffix=${i##*.}
    if [ $suffix = "gz" ]; then
        header_num=`gzip -dc $i | grep -nv '^#' | head -n 1 | cut -d ':' -f 1`
        gzip -dc $i | tail -n +$header_num >> ${opt_prefix}.snp_indel.tmp.vcf
    else
        header_num=`grep -nv '^#' $i | head -n 1 | cut -d ':' -f 1`
        tail -n +$header_num $i >> ${opt_prefix}.snp_indel.tmp.vcf
    fi
done

# step 2 sort vcf using picard
$java8 -jar $picard SortVcf I=${opt_prefix}.snp_indel.tmp.vcf SEQUENCE_DICTIONARY=$opt_dict O=${opt_prefix}.snp_indel.raw.vcf

# step 3 split SNP and Indel
mkdir $opt_dirname/JavaTmpDir
$java8 -Xmx${mem}G -XX:+UseSerialGC -Djava.io.tmpdir=$opt_dirname/JavaTmpDir -jar $gatk3 -T SelectVariants -R $opt_genome -V ${opt_prefix}.snp_indel.raw.vcf -selectType SNP -o ${opt_prefix}.snp.raw.vcf &

$java8 -Xmx${mem}G -XX:+UseSerialGC -Djava.io.tmpdir=$opt_dirname/JavaTmpDir -jar $gatk3 -T SelectVariants -R $opt_genome -V ${opt_prefix}.snp_indel.raw.vcf -selectType INDEL -o ${opt_prefix}.indel.raw.vcf &
wait

# step 4 get HQ SNP
perl $Bin/get_HQ_site.pl ${opt_prefix}.snp.raw.vcf ${opt_prefix}.snp.HQ.vcf &
perl $Bin/get_HQ_site.pl ${opt_prefix}.indel.raw.vcf ${opt_prefix}.indel.HQ.vcf &
wait
$bgizp -@ $opt_threads ${opt_prefix}.snp.HQ.vcf
$bgizp -@ $opt_threads ${opt_prefix}.indel.HQ.vcf
$tabix -p vcf ${opt_prefix}.snp.HQ.vcf.gz
$tabix -p vcf ${opt_prefix}.indel.HQ.vcf.gz
