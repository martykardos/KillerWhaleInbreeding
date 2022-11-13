mkdir 03.recombine_EVM_partial
cd 03.recombine_EVM_partial
echo "/dellfsqd2/ST_OCEAN/USER/songyue/bin/03.software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions ../01.Partitioning/partitions_list.out --output_file_name evm.out " > recombine_EVM_partial_output.sh
qsub -cwd -l vf=1g,num_proc=1 -binding linear:1 -P P18Z19700N0073 -q st.q recombine_EVM_partial_output.sh
cd ../


