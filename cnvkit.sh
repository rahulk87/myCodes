#!/bin/bash

~/miniconda3/bin/python cnvkit.py batch ~/share/data/myoep/impact468/bam/AM1T.bam --normal ~/share/data/myoep/impact468/bam/AM1N.bam --targets my_run/impact468.bed --annotate my_run/refFlat_mod.txt --fasta /home/gularter/genomes/homo_sapiens/Ensembl/Grch37.p13/Sequence/WholeGenomeFasta/genome.fa --output-reference my_run/AM1.cnn --output-dir my_run/results --diagram --scatter

~/miniconda3/bin/python cnvkit.py scatter my_run/results/AM5TC.cnr -s my_run/results/AM5TC.cns -o my_run/results/AM5TC-scatter.pdf

~/miniconda3/bin/python cnvkit.py heatmap my_run/results/*.cns -o my_run/results/combined_heatmap.pdf
