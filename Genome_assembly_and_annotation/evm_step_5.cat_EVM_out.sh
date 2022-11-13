mkdir 05.cat_evm

find ./01.Partitioning -regex ".*evm.out.gff3" -exec cat {} \; > ./05.cat_evm/evm.all.gff3
perl xx/gff.simple.pl -gff ./05.cat_evm/evm.all.gff3 > ./05.cat_evm/evm.all.gff
p





