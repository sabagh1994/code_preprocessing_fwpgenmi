#!/bin/bash

# addr= /shared-mounts/sinhas/Saba/software/bedtools2 => where ./bedtools 2.28 could be found :)
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sabag2/anaconda3/envs/bed_env/lib/ ./bedtools closest

thr="200Kb"


# Intersect coords files
../../code/intersect.pl ../coords/tfs_intergenic_protein_coding_bed_coords.tsv > ./tfs_intergenic_protein_coding.bed10
../../code/intersect.pl ../coords/tfs_intron_protein_coding_bed_coords.tsv > ./tfs_intron_protein_coding.bed10

# Sort files
for i in *.bed10; do 
	sort -k1,1 -k2,2n $i -o $i.sorted
done

# Get closest gene for integenic
# you have to be in intersection folder already

module load glibc
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sabag2/anaconda3/envs/bed_env/lib/ ./bedtools closest -a tfs_intergenic_protein_coding.bed10.sorted -b ../gencode.v27lift37.protein_coding.gene.bed4.sorted -d -k 20181 > coords/tfs_intergenic_protein_coding_closest_gene_coords.tsv
module unload glibc

# Get closest gene within 1MB
awk -F $'\t' 'BEGIN{OFS=FS} {if ($14 < 200000) print $1,$2,$3,$4"|"$13,$5,$6,$7,$8,$9}' coords/tfs_intergenic_protein_coding_closest_gene_coords.tsv >> tfs_intergenic_protein_coding_with_closest_gene_$thr.bed10
sort -k1,1 -k2,2n tfs_intergenic_protein_coding_with_closest_gene_$thr.bed10 -o tfs_intergenic_protein_coding_with_closest_gene_$thr.bed10.sorted

# Concatenate intronic and intergenic bed10 files into final file
cat tfs_intergenic_protein_coding_with_closest_gene_$thr.bed10.sorted > ../tfs_protein_coding_$thr.bed10
cat tfs_intron_protein_coding.bed10.sorted >> ../tfs_protein_coding_$thr.bed10
sort -k1,1 -k2,2n ../tfs_protein_coding_$thr.bed10 -o ../tfs_protein_coding_$thr.bed10.sorted

