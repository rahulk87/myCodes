#!/bin/bash
#BSUB -u rkumar\@icr.ac.uk
#BSUB -J ATF2
#BSUB -e /scratch/DBC/GENFUNC/rkumar/extra/dream/tf/pred_on_dnase/dnase_motif_genic_pred/ATF2/tf.err
#BSUB -o /scratch/DBC/GENFUNC/rkumar/extra/dream/tf/pred_on_dnase/dnase_motif_genic_pred/ATF2/tf.out
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -P DBCDOBZAK
#BSUB -W 96:00
#BSUB -q normal

#gunzip *.gz

input1="list1"
while IFS=$'\t' read -r -a input1
do
tf1="${input1[0]}"
l_cell="${input1[1]}"
t_cell="${input1[2]}"


input="list"
while IFS=$'\t' read -r -a input
do
tf="${input[0]}"
train="${input[1]}"
pred="${input[2]}"
pred_l="${input[3]}"
te_cell="${input[4]}"
pred_te="${input[5]}"


## preparation of training cell line data
#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/motif/train/train.$tf.$train\_motif.bed ../../processed/dnase/train/Tr.$tf.$train\.dnase |paste - ../../processed/motif/train/train.$tf.$train\_motif.bed |cut -f1-7,11,12 >$tf.$train\.dnase.motif

#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/genic/train_genic.bed $tf.$train\.dnase.motif |paste - ../../processed/genic/train_genic.bed |cut -f1-9,13-15 >$tf.$train\.dnase.motif.genic

##../../../../software/bedtools2/bin/bedtools intersect  -a Tr.ATF2.GM12878.dnase.bed -b train.ATF2.GM12878_motif.bed -wao -f 1 >kd

##../../../../software/bedtools2/bin/bedtools intersect  -a kd -b ../../processed/genic/train_genic.bed -wao -f 1 >kd1

#cat header $tf.$train\.dnase.motif.genic >tr.$tf.$train\.dnase.motif.genic

#grep -v $'^\t' tr.$tf.$train\.dnase.motif.genic >trr.$tf.$train\.dnase.motif.genic

#rm $tf.$train\.dnase.motif $tf.$train\.dnase.motif.genic tr.$tf.$train\.dnase.motif.genic


## preparation of leaderboard cell line data
#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/motif/leader/ladder.$tf.$train\_motif.bed ../../processed/dnase/leader/L.$pred\.dnase |paste - ../../processed/motif/leader/ladder.$tf.$train\_motif.bed |cut -f1-7,11,12 >$tf.$pred\.dnase.motif

#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/genic/ladder_genic.bed $tf.$pred\.dnase.motif |paste - ../../processed/genic/ladder_genic.bed |cut -f1-9,13-15 >$tf.$pred\.dnase.motif.genic

#cat header $tf.$pred\.dnase.motif.genic >pred.$tf.$pred.$train\.dnase.motif.genic

#rm $tf.$pred\.dnase.motif $tf.$pred\.dnase.motif.genic


## preparation of test cell line data
#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/motif/test/test.$tf.$train\_motif.bed ../../processed/dnase/test/Te.$t_cell\.dnase |paste - ../../processed/motif/test/test.$tf.$train\_motif.bed |cut -f1-7,11,12 >$tf.$t_cell\.dnase.motif

#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' ../../processed/genic/test_genic.bed $tf.$t_cell\.dnase.motif |paste - ../../processed/genic/test_genic.bed |cut -f1-9,13-15 >$tf.$t_cell\.dnase.motif.genic

#cat header $tf.$t_cell\.dnase.motif.genic >test.$tf.$t_cell.$train\.dnase.motif.genic

#rm $tf.$t_cell\.dnase.motif $tf.$t_cell\.dnase.motif.genic

/home/breakthr/rkumar/R/bin/Rscript glm.R trr.$tf.$train\.dnase.motif.genic pred.$tf.$pred.$train\.dnase.motif.genic $pred_l test.$tf.$t_cell.$train\.dnase.motif.genic $pred_te

done < "$input"

## preparing leader cell line to submit
paste pred_* >prediction_l

cat prediction_l |awk '{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}' >max_pred

cut -f1,2,3 pred.$tf.$pred.$train\.dnase.motif.genic |grep -v 'start' >tmp

paste tmp max_pred >max_pred1

#awk -F'\t' 'NR==FNR{c[$1$2$3]++;next};c[$1$2$3] > 0' max_pred1 ../../../../annotations/ladder_regions.blacklistfiltered.bed |paste - max_pred1 >tmp3

../../../software/bedtools2/bin/bedtools intersect  -a ../../../annotations/ladder_regions.blacklistfiltered.bed -b max_pred1 -wao -f 1 >tmp3

cat tmp3 |awk 'OFS="\t" {print$1,$2,$3,1-$7}' >tmp4

mv tmp4 L.$tf1\.$l_cell\.tab

gzip L.$tf1\.$l_cell\.tab

rm max_pred* tmp tmp3


## preparing test cell line to submit
paste test_* >prediction_te

cat prediction_te |awk '{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i;print m}' >max_pred

cut -f1,2,3 test.$tf.$t_cell.$train\.dnase.motif.genic |grep -v 'start' >tmp

paste tmp max_pred >max_pred1

../../../software/bedtools2/bin/bedtools intersect  -a ../../../annotations/test_regions.blacklistfiltered.bed -b max_pred1 -wao -f 1 >tmp3

cat tmp3 |awk 'OFS="\t" {print$1,$2,$3,1-$7}' >tmp4

mv tmp4 F.$tf1\.$t_cell\.tab

gzip F.$tf1\.$t_cell\.tab

rm max_pred* tmp tmp3

done < "$input1" 
