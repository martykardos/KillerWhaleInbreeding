#!/usr/bin/python

import sys,re,os

###  读取样本id
ID_file = open(r"kw_148.snp_QC_ld_nj.mdist.id","r")
arr = []
for i in ID_file:
    arr.append(i.strip().split()[0])
head = map(str,range(1,len(arr)+1))
ID_file.close()

### 读取样本间IBS距离
dist_file =open(r"kw_148.snp_QC_ld_nj.mdist","r").readlines()
print ("#mega\n!Title: fasta file;\n!Format DataType=Distance DataFormat=LowerLeft NTaxa=%d;\n" % len(arr))
for i,j in enumerate(arr):
    print ('['+str(i+1)+']'+' #'+j)
#print

print ('[    '+'    '.join(head)+' ]')
print ('[1]')
for l,m in enumerate(dist_file):
    tmp = m.strip().split()
    #tmp[-1] = ''
    print ('['+str(l+2)+']    '+'    '.join(tmp))

#dist_file.close()
