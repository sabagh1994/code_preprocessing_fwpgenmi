#!/bin/sh

run1=1;
run2=1;
run3=1;
run4=1;

codedir="/home/casey/offer/analysis/BedToolsDiffExp/code";
basedir="/home/casey/offer/analysis/BedToolsDiffExp/K27me3";
r1dir="/home/casey/offer/data/encode_narrowPeak/K27me3/Rep1";
r2dir="/home/casey/offer/data/encode_narrowPeak/K27me3/Rep2";
p0r1="$r1dir/SW480-A-K27me3-1.FCC7C5GACXX_L8_R1_IAACCAG.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p2r1="$r1dir/P2-B-K27me3-2.FCC7C5GACXX_L8_R1_ITGGTGA.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p4r1="$r1dir/P4-C-K27me3-3.FCC7C5GACXX_L8_R1_IAGTGAG.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p6r1="$r1dir/P6-D-K27me3-4.FCC7C5GACXX_L8_R1_IGCACTA.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p0r2="$r2dir/2SW480-A-K27me3-1.FCC7WFRACXX_L5_R1_IAACCAG.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p2r2="$r2dir/2P2-B-K27me3-2.FCC7WFRACXX_L5_R1_ITGGTGA.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p4r2="$r2dir/2P4-C-K27me3-3.FCC7WFRACXX_L5_R1_IAGTGAG.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";
p6r2="$r2dir/2P6-D-K27me3-4.FCC7WFRACXX_L5_R1_IGCACTA.PE.U12.dedup.s1-W200-G600-FDR1E-2.narrowPeak";

## Step 1: Merge Replicates: Produce x/x_intersections.narrowPeak files
if [ $run1 -eq 1 ]
then
	$codedir/run_replicate_intersection.sh $p0r1 $p0r2 p0 $basedir/p0
	$codedir/run_replicate_intersection.sh $p2r1 $p2r2 p2 $basedir/p2
	$codedir/run_replicate_intersection.sh $p4r1 $p4r2 p4 $basedir/p4
	$codedir/run_replicate_intersection.sh $p6r1 $p6r2 p6 $basedir/p6
fi

## Step 2: Merge 0-2 and 4-6: Produce y/y.narrowPeak and y/y.merged.narrowPeak files
p0="$basedir/p0/p0_intersections.narrowPeak";
p2="$basedir/p2/p2_intersections.narrowPeak";
p4="$basedir/p4/p4_intersections.narrowPeak";
p6="$basedir/p6/p6_intersections.narrowPeak";
if [ $run2 -eq 1 ]
then
	$codedir/run_merge_timepoints.sh $p0 $p2 m02 $basedir/m02
	$codedir/run_merge_timepoints.sh $p4 $p6 m46 $basedir/m46
fi

## Step 3: Aggregate Merged Statistics: m02.final.narrowPeak and m46.final.narrowPeak
m02_merged="$basedir/m02/m02.merged.narrowPeak";
m46_merged="$basedir/m46/m46.merged.narrowPeak";
m02_final="$basedir/m02/m02.final.narrowPeak";
m46_final="$basedir/m46/m46.final.narrowPeak";
if [ $run3 -eq 1 ]
then
	$codedir/aggregate_merged.pl $m02_merged > $m02_final;
	$codedir/aggregate_merged.pl $m46_merged > $m46_final;
fi

## Step 4: Differential Expression
if [ $run4 -eq 1 ]
then
	$codedir/run_de.sh $m02_final $m46_final de_m02_m46 $basedir/final;
fi
