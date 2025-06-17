#!/bin/sh


set -e

run1=0;

if [ $run1 -eq 1 ]
then
	addr='/shared-mounts/sinhas-storage1/mayo/offer_project/pipeline_saba/pipeline/overlap/encode_gencode/TFBS-only-analysis-saba/no-nearest-gene/'
	ln -s $addr/tfs_protein_coding_10Kb.bed10 ./
	ln -s $addr/tfs_protein_coding_50Kb.bed10 ./
	ln -s $addr/tfs_protein_coding_200Kb.bed10 ./
	ln -s $addr/tfs_protein_coding_1Mb.bed10 ./
fi

tag=$1
basedir='/shared-mounts/sinhas-storage1/mayo/offer_project/recent_preprocessing/pipeline/ACC_analysis'
mkdir -p intersections/$tag/
mkdir -p intersections/$tag/tfSeparation

# there is no need to sort the data by chromosome and the start position but if you are 
# dealing with large files you should sort them and use -sorted option according to 
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

bedtools intersect -a tfs_protein_coding_$tag.bed10 -b diffAcc/ACC.diff.peak.sorted -wa -wb > intersections/$tag/acc_tfbs_overlap


# filter by (TF, TARGET, DIFF PEAK SIGNAL)
awk -F $'\t' '{x=gsub(/\|/,"\t",$4); print $4"\t"$13}' intersections/$tag/acc_tfbs_overlap > intersections/$tag/acc_tfbs_overlap.separated

# get maximum signal
cd intersections/$tag/
mkdir -p max_by_tf/continuous max_by_tf/binary
$basedir/get_max_by_tf.pl acc_tfbs_overlap.separated > ./max_by_tf/continuous/acc_tfbs_overlap
cd max_by_tf/continuous
#convert to pgenmi bit vector with two fields up/down
for i in *; do awk '{if ($3 < 0) { print $1"\t"$2"\t0\t1"} else { print $1"\t"$2"\t1\t0"} }' $i >> ../binary/$i; done

#you can separate by tfs as well but there is no need and if you want to use the following code make the 
# appropriate changes though
cd $basedir/intersections/$tag/max_by_tf/binary
for i in *; do awk -v odir=$basedir/intersections/$tag/tfSeparation/ '{print $2"\t"$3"\t"$4 > odir"/"$1}' $i; done 
