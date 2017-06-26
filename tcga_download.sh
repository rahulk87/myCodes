#!/bin/bash
##PBS -q jrf.q
#PBS -N tcga_download
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
# -j oe
#PBS -l mem=32G
##PBS -V -m n
##PBS -t 0-127%10

path='/cbio/ski/reis-filho/home/limr/share/data/gct/TCGA/test'

/cbio/ski/reis-filho/home/gularter/.local/bin/gdc-client download -m $path/gdc_manifest.2017-06-25T17-21-36.934493.tsv -t /cbio/ski/reis-filho/home/gularter/.gdc-token/gdc-jrf-token.2017-06-08T10_25_18-04_00.txt
