#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os,sys
import gzip
def find_digest(seq,digest_seq):
	i = 0
	l_loc = []
	seq_len = len(seq)
	seq = seq.upper()
	while True:
		try:
			index = seq.find(digest_seq)
		except:
			print("The genome file is abnormal,please check!")
		if len(l_loc) < 1:
			i += index
		else:
			i += index + len(digest_seq)
		if index == -1:
			break
		else:
			l_loc.append(str(i))
			seq = seq[index + len(digest_seq):]
	if l_loc == [] or l_loc[0] != "0":
		l_loc.insert(0,"0")
	if i != seq_len:
		l_loc.append(str(seq_len))
	return l_loc

def print_out(chr_id,lis_loc,out):
	for i in range(1,len(lis_loc)):
		out.write(chr_id + "\t" + lis_loc[i-1] + "\t" + lis_loc[i] + "\t" + "HIC_" + chr_id + "_" + str(i) + "\t0\t+\n")

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python digest_genome_v1.py input output digest[mboi,dpnii,hindiii]")
		sys.exit()	
	input_f = sys.argv[1]
	output_f = sys.argv[2]
	digest_name = sys.argv[3].lower()
	d_digest = {"mboi":"GATC","dpnii":"GATC","hindiii":"AGCTT"}
	try:
		digest_seq = d_digest[digest_name]
	except:
		print("Digest maybe unsupported")
		sys.exit()
	out = open(output_f,"w")
	chr_id = None
	
	if input_f.endswith(".gz"):
		try:
			f = gzip.open(input_f)
		except:
			print("Fail to open",input_f)
			sys.exit()
	else:
		try:
			f = open(input_f)
		except:
			print("Fail to open",input_f)
			sys.exit()

	for line_raw in f:
		line = line_raw.strip()
		if line.startswith(">"):
			if chr_id:
				l_return = find_digest(seq,digest_seq)
				print_out(chr_id,l_return,out)	
			chr_id = line.split()[0][1:]
			seq = ""
		else:
			seq += line
	l_return = find_digest(seq,digest_seq)
	print_out(chr_id,l_return,out)
	f.close()
	out.close()

