#!/bin/bash


ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline/gencode/gencode.v27lift37.protein_coding.gene.bed4.sorted .
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline/gencode/gencode.v27lift37.protein_coding.intron.bed4.sorted .
ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline/gencode/gencode.v27lift37.protein_coding.intergenic.bed4.sorted .

# Make dir coords (which contain overlap coords for both files) and intersection (which is just the bed file)
mkdir coords intersection
tfs_dir="/shared-mounts/sinhas-storage1/mayo/offer_project/OtherCelline_Analysis/ENCODE_OtherOrgan_BLOOD/"
intersectBed -a $tfs_dir/tfs.bed10.sorted -b gencode.v27lift37.protein_coding.intron.bed4.sorted -wa -wb > coords/tfs_intron_protein_coding_bed_coords.tsv
intersectBed -a $tfs_dir/tfs.bed10.sorted -b gencode.v27lift37.protein_coding.intergenic.bed4.sorted -wa -wb > coords/tfs_intergenic_protein_coding_bed_coords.tsv

# Intersect coords files
code/intersect.pl coords/tfs_intergenic_protein_coding_bed_coords.tsv > intersection/tfs_intergenic_protein_coding.bed10
code/intersect.pl coords/tfs_intron_protein_coding_bed_coords.tsv > intersection/tfs_intron_protein_coding.bed10

# Sort files
cd intersection
for i in *.bed10; do
        sort -k1,1 -k2,2n $i -o $i.sorted
done
