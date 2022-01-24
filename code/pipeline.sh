#!/bin/bash

R2_suffix=_2.fastq.gz
assemble_suffix=.assembled.fastq
fastq_prefix=../fastq/
output_prefix=../result/
outout_suffix=.tsv
for f in ../fastq/*_1.fastq.gz
do
	name=`basename $f "_1.fastq.gz"`
	r="$fastq_prefix$name$R2_suffix"
	pear -f $f -r $r -o $output_prefix$name -j 20
	python ./fastq2count.py $output_prefix$name$assemble_suffix ../ref/WT_seq.fasta $output_prefix$name$outout_suffix
done