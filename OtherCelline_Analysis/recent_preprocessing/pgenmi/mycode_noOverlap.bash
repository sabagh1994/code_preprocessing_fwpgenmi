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
basedir_othercelline='/shared-mounts/sinhas-storage1/mayo/offer_project/OtherCelline_Analysis/recent_preprocessing/'
evidaddr="$basedir_othercelline/epi_analysis/output_hg19_noOverlap/intersections/$tag/by_mark_$tag/binary/"
datadir="$basedir/updated_data/"

cd $evidaddr

for i in *; do for j in $i/*; do $basedir_othercelline/pgenmi/make_input.pl $basedir_othercelline/pgenmi/protein_coding_genes.list $j $datadir/DESeq2_results/main_comparison_p0vsp6/p0vsp6_DESeq_processed > $basedir_othercelline/pgenmi/binary_two_feature_noOverlap/$tag/$j; done done 
