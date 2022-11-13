#!/bin/bash
/share/app/bwa-0.7.12/bwa index -a bwtsw genome.fa

/share/app/bwa-0.7.12/bwa mem -t 60 genome.fa wgs_1.clean.fq.gz wgs_2.clean.fq.gz | /share/app/samtools-1.9/bin/samtools view -@ 60 -Sb - | /share/app/samtools-1.9/bin/samtools sort -@ 60 - -o align.bam
/share/app/samtools-1.9/bin/samtools index -@ 60 align.bam
