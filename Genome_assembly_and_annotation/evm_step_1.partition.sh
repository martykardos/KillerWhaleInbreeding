mkdir 01.Partitioning
cd 01.Partitioning

echo "xxn/03.software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome xx/ref_genome.fa --gene_predictions xx/denovo.augustus.gff --protein_alignments xx/homolog.gff3  --segmentSize 1000000 --overlapSize 200000 --partition_listing partitions_list.out; echo Finished" > partition.sh

qsub -cwd -l vf=1g,num_proc=1 -binding linear:1 -P P18Z19700N0073 -q st.q partition.sh
cd ../
