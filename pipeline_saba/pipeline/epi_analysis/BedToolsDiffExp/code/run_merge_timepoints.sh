#!/bin/sh

mkdir -p $4
cat $1 >> $4/$3.narrowPeak
cat $2 >> $4/$3.narrowPeak
sort -k1,1 -k2,2n $4/$3.narrowPeak -o $4/$3.narrowPeak
bedtools merge -i $4/$3.narrowPeak -c 4,5,6,7,8,9,10 -o distinct,collapse,distinct,collapse,collapse,collapse,collapse | sort -k1,1 -k2,2n > $4/$3.merged.narrowPeak
