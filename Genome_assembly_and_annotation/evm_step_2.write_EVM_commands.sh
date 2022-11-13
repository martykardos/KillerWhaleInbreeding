mkdir 02.CommandSet
cd 02.CommandSet
/dellfsqd2/ST_OCEAN/USER/songyue/bin/03.software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --partitions ../01.Partitioning/partitions_list.out --genome xx/ref_genome.fa --gene_predictions xx/denovo.augustus.gff --protein_alignments xx/homolog.gff3  --weights xx/weights.txt --output_file_name evm.out > CommandSet.sh
/dellfsqd2/ST_OCEAN/USER/yangxianwei/bin/SGE/qsub-sge.pl  --queue st.q -P P18Z19700N0073 --lines 1 --convert no --maxjob 120 --resource vf=1g,num_proc=1 CommandSet.sh &

cd ../
