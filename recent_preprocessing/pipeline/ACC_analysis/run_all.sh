#!/bin/sh

rm -rf intersections
./tfbs_acc_intersection.sh 200Kb
./tfbs_acc_intersection.sh 10Kb
./tfbs_acc_intersection.sh 50Kb
./tfbs_acc_intersection.sh 1Mb
