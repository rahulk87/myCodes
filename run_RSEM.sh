#! /bin/bash
# set echo
##$ -q jrf.q
#$ -w e -N RSEM
#$ -l h_vmem=16G
#$ -pe smp 5 -j yes
#$ -V -cwd
##$ -t 1-12

#/home/kumarr/software/RSEM/software/bin/rsem-prepare-reference \
#--gtf ~/genomes/Homo_sapiens.GRCh38.91.gtf \
#--bowtie ~/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa ref/human_ensembl

~/software/RSEM/software/bin/rsem-calculate-expression -p 8 --paired-end --strandedness reverse fastqs/GCT10_R1.fastq.gz fastqs/GCT10_R2.fastq.gz ref/human_ensembl GCT10
