#!/bin/bash
#=============================================================================

# Filename: VCF.final_filter.sh
#=============================================================================
Bin=`cd \`dirname $0\`;pwd`

# help info
help(){
    cat <<HELP
    VCF.final_filter.sh -- filter VCF using GATK3.8

    [usage]
     VCF.final_filter.sh -g|--genome <genome.fa> -v|--vcf <in.vcf>
                
    [option]
     -p|--prefix   [str]   out prefix [./out]
     -c|--class    [str]   SNP/INDEL [SNP]
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
        -c|--class)     opt_class=$2;   shift 2;;
        -p|--prefix)    opt_prefix=$2;  shift 2;;
        -t|--threads)   opt_threads=$2; shift 2;;
        -*)             echo "*** no option $1, -h for help"; exit 1;;
        *)              break;; # not to deal
    esac
done

# parameter check and initialize
[ -z "$opt_prefix" ] && opt_prefix='./out'
[ -z "$opt_threads" ] && opt_threads=1
[ -z "$opt_class" ] && opt_class='SNP'

# deal with file path
opt_genome=`readlink -f $opt_genome`
opt_vcf=`readlink -f $opt_vcf`
opt_dirname=`cd \`dirname $opt_prefix\`;pwd`
opt_basname=`basename $opt_prefix`
opt_prefix="$opt_dirname/$opt_basname"

# default
java8=./java1.8/bin/java
gatk3=./software/variation/GATK/GenomeAnalysisTK-3.8-1.jar
bgzip=./software/variation/htslib-1.4.1/bin/bgzip
tabix=./software/variation/htslib-1.4.1/bin/tabix

mem=10  # memory 10G

set -x
# step 1. merge vcf
mkdir -p $opt_dirname/JavaTmpDir

# for SNP filter
if [ $opt_class = "SNP" ]; then
    $java8 -Xmx${mem}G -XX:+UseSerialGC -Djava.io.tmpdir=$opt_dirname/JavaTmpDir -jar $gatk3 -T VariantFiltration -R $opt_genome --variant $opt_vcf -o ${opt_prefix}.snp.filtered.vcf --filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || ReadPosRankSum < -8.0 || MQRankSum < -12.5 || SOR >3.0" --filterName LowQualFilter
    awk '$1~/^#/ || $7=="PASS"' ${opt_prefix}.snp.filtered.vcf > ${opt_prefix}.snp.final.vcf
    $bgzip -@ $opt_threads -f ${opt_prefix}.snp.final.vcf
    $tabix -p vcf ${opt_prefix}.snp.final.vcf.gz
fi

# for INDEL filter
if [ $opt_class = "INDEL" ]; then
    $java8 -Xmx${mem}G -XX:+UseSerialGC -Djava.io.tmpdir=$opt_dirname/JavaTmpDir -jar $gatk3 -T VariantFiltration -R $opt_genome --variant $opt_vcf -o ${opt_prefix}.indel.filtered.vcf --filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" --filterName LowQualFilter
    awk '$1~/^#/ || $7=="PASS"' ${opt_prefix}.indel.filtered.vcf > ${opt_prefix}.indel.final.vcf
    $bgzip -@ $opt_threads -f ${opt_prefix}.indel.final.vcf
    $tabix -p vcf ${opt_prefix}.indel.final.vcf.gz
fi
