unset R_LIBS PYTHONPATH LD_LIBRARY_PATH PERL5LIB
export PATH=XX/software/conda/bin/:$PATH
export PATH=/share/app/bowtie2-2.2.5/:$PATH
export PATH=/share/app/R-3.4.1/bin/:$PATH
export PATH=/share/biosoft/pipeline/DNA/HiC_Assembly/01.HiCPro_v2/bin:$PATH
export PATH=/share/biosoft/pipeline/DNA/HiC_Assembly/bin:$PATH
export PYTHONPATH=xx/software/conda/lib/python3.6/site-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/share/app/gcc-5.2.0/lib64
export PERL5LIB=/usr/lib64/perl5
python /share/biosoft/pipeline/DNA/HiC_Assembly/bin/digest_genome_v1.py xx/ref_genome.fa xx/ref_genome.fa.mboi.bed mboi
/share/biosoft/pipeline/DNA/HiC_Assembly/bin/fa_stat.pl xx/ref_genome.fa > xx/ref_genome.fa.len
bowtie2-build xx/ref_genome.fa xx/ref_genome.fa
/share/biosoft/pipeline/DNA/HiC_Assembly/bin/HiC-Pro_v3.2 -c xx/hicpro.config -i xx -o xx