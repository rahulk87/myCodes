#!/bin/bash

function run_macs2 {
	~/miniconda2/bin/macs2 callpeak -t ../../../bam/$1.bam -c ../../../bam/$2.bam --nomodel -f BAM -g hs -n $1 -B -q 0.1
}
export -f run_macs2

run_macs2 PC-si PC-TS15 &
run_macs2 PX2-si PX2-TS16 &
run_macs2 CHP18-si CHP18-TS18 &
run_macs2 CHP30-si CHP30-TS19 &
