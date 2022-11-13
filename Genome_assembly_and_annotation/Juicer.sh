unset R_LIBS PYTHONPATH LD_LIBRARY_PATH PERL5LIB
export PATH=/share/app/bwa-0.7.12/:$PATH
export PATH=./conda/bin/:$PATH
export PERL5LIB=/usr/lib64/perl5
ln -s XX/polish.genome/ref_genome.fa ./genome/ref_genome
ln -s XX/read.1.fastq.gz ./fastq/
ln -s XX/read.2.fastq.gz ./fastq/
bwa index genome/xquas098_genome
python XX/generate_site_positions.py MboI genome/ref_genome genome/ref_genome
XXfastaDeal.pl --attribute ID:len genome/ref_genome > genome/ref_genome.size.txt
mkdir tmp
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=XX/06.Juicer/tmp -XX:ParallelGCThreads=12"
export TMPDIR=XX06.Juicer/tmp
XX/juicer_v3.sh -z genome/ref_genome -y genome/ref_genome_MboI.txt -p genome/ref_genome.size.txt -s MboI -d XX/06.Juicer -D /share/biosoft/pipeline/DNA/HiC_Assembly/04.Juicer -T 12 -S early -m 30G
