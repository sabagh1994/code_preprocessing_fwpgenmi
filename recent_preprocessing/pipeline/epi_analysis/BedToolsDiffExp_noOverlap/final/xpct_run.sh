#!/bin/sh


mkdir -p narrowPeak/10pct
./xpct.sh narrowPeak/all/K4me1.narrowPeak 0.1 narrowPeak/10pct/K4me1.narrowPeak
./xpct.sh narrowPeak/all/K4me3.narrowPeak 0.1 narrowPeak/10pct/K4me3.narrowPeak
./xpct.sh narrowPeak/all/K27ac.narrowPeak 0.1 narrowPeak/10pct/K27ac.narrowPeak
./xpct.sh narrowPeak/all/K27me3.narrowPeak 0.1 narrowPeak/10pct/K27me3.narrowPeak
