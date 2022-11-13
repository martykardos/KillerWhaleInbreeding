mkdir 04.Convert
cd 04.Convert

/dellfsqd2/ST_OCEAN/USER/songyue/bin/03.software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions ../01.Partitioning/partitions_list.out --output_file_name evm.out --genome xx/ref_genome.fa' > convert_EVM_outputs.sh

qsub -cwd -l vf=1g,num_proc=1 -binding linear:1 -P P18Z19700N0073 -q st.q convert_EVM_outputs.sh

cd ../

