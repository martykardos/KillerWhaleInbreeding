unset R_LIBS PYTHONPATH LD_LIBRARY_PATH PERL5LIB
export PATH=/share/app/lastz-1.03.73/:$PATH
export PATH=/zfsqd1/app/parallel-20150322/bin:$PATH
export PATH=/share/biosoft/pipeline/DNA/HiC_Assembly/bin:$PATH
export PERL5LIB=/usr/lib64/perl5
export PYTHONPATH=XX/software/conda/lib/python3.6/site-packages/
ln -s XX/06.Juicer/aligned/merged_nodups.txt .
mkdir tmp
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=XX/07.3ddna/tmp -XX:ParallelGCThreads=10"
export TMPDIR=XX/07.3ddna/tmp
XX/3ddna_v1.sh -m haploid -s 0 -c 22 -j 10 -S init- XX/ref_genome.fa merged_nodups.txt 200G
