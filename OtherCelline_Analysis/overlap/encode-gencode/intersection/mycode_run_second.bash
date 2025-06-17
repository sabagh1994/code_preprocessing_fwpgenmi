#!/bin/bash


# you have to be in intersection folder already

module load glibc
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sabag2/anaconda3/envs/bed_env/lib/ ./bedtools closest -a tfs_intergenic_protein_coding.bed10.sorted -b ../gencode.v27lift37.protein_coding.gene.bed4.sorted -d -k 50 > coords/tfs_intergenic_protein_coding_closest_gene_coords.tsv
module unload glibc



