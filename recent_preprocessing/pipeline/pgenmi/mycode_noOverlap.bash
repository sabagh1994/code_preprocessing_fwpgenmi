#!/bin/bash
set -e

run1=0
if [ $run1 -eq 1 ]
then
	ln -s /shared-mounts/sinhas-storage1/mayo/offer_project/pipeline/gencode/protein_coding_genes.list .
fi

#tag="200Kb"
tag=$1

mkdir -p binary_two_feature_noOverlap/$tag
for i in `cat marks.list`; do mkdir -p binary_two_feature_noOverlap/$tag/$i; done


basedir='/shared-mounts/sinhas-storage1/mayo/offer_project/recent_preprocessing'
evidaddr="$basedir/pipeline/epi_analysis/output_hg19_noOverlap/intersections/$tag/by_mark_$tag/binary/"
datadir="$basedir/updated_data/"

cd $evidaddr

for i in *; do for j in $i/*; do $basedir/pipeline/pgenmi/make_input.pl $basedir/pipeline/pgenmi/protein_coding_genes.list $j $datadir/DESeq2_results/main_comparison_p0vsp6/p0vsp6_DESeq_processed > $basedir/pipeline/pgenmi/binary_two_feature_noOverlap/$tag/$j; done done 
