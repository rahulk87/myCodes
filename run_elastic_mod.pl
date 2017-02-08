#### Run elastic net to select gene features and performance ####
#!/usr/bin/perl
for($i=1;$i<=100;$i++)
{
print "Running for $i....\t";
open(FD,">>perf_$i");
`head -58 input.txt >train`;
`sed '59,69!d' input.txt >test`;
`cut -f3- train |grep -v '[A-Z]' >x`;
`cut -f2 train |grep -v '[A-Z]' >y`;
`cut -f3- test |grep -v '[A-Z]' >x1`;
`cut -f2 test |grep -v '[A-Z]' >y1`;
`cut -f3- input.txt |grep -v '[A-Z]' >x2`;
`cut -f2 input.txt |grep -v '[A-Z]' >y2`;
`head -1 train |perl -pi -e 's/\t/\n/g' |grep -v 'Cell_line' |grep -v 'AUC' >genes`;
##################
	open(FA,">>glmnet.R");
	print FA "library(glmnet)\n";
	print FA "x<-as.matrix(read.csv(\"x\",sep=\"\\t\",header=FALSE))\n";
	print FA "y<-as.matrix(read.csv(\"y\",sep=\"\\t\",header=FALSE))\n";
	print FA "cvfit = cv.glmnet(x, y,nfolds=10)\n";
	print FA "lmin = cvfit\$lambda.min\n"; #### lambda minimum
	print FA "write.table(lmin,file=\"lmin\")\n";
	close FA;
	`/scratch/breakthr/rkumar/R/bin/R CMD BATCH glmnet.R`;
	$lmin=`sed '2!d' lmin |perl -pi -e 's/"1" //g'`;chomp $lmin;
	for($a=0;$a<=1;$a+=0.1)
        {
	print "$a\t";
	open(FE,">>glmnet1.R");
	print FE "x<-as.matrix(read.csv(\"x\",sep=\"\\t\",header=FALSE))\n";
        print FE "y<-as.matrix(read.csv(\"y\",sep=\"\\t\",header=FALSE))\n";
	print FE "x1<-as.matrix(read.csv(\"x1\",sep=\"\\t\",header=FALSE))\n";
	print FE "library(glmnet)\n";
	print FE "fit = glmnet(x, y,lambda=$lmin,alpha=$a)\n";
	print FE "a<-predict(fit, newx = x1, s = \"$lmin\")\n";
	print FE "b<-as.matrix(a)\n";
	print FE "write.table(b,file=\"c\")\n";
        close FE;
        `/scratch/breakthr/rkumar/R/bin/R CMD BATCH glmnet1.R`;
        ##############
	open(FB,">>cor.R");
        print FB "library(pspearman)\n";
        print FB "library(MASS)\n";
        print FB "library(hydroGOF)\n";
        print FB "a=read.table(\"act\")\n";
        print FB "c=read.table(\"pred\")\n";
        print FB "aa=as.matrix(a)\n";
        print FB "cc=as.matrix(c)\n";
        print FB "x=cor.test(aa,cc)\n";
        print FB "write.matrix(x, file=\"cor\")\n";
        print FB "rmse=rmse(aa,cc)\n";
        print FB "write.table(rmse, file=\"rmse\")\n";
        close FB;
	`cp y1 act`;
        `cut -f2 -d ' ' c |grep -v '"' >pred`;
	`/scratch/breakthr/rkumar/R/bin/R CMD BATCH cor.R`;
        $r=`sed 4!d cor`;chomp $r;
        $rmse=`sed '2!d' rmse |perl -pi -e 's/"V1" //g'`;chomp $rmse;
	print FD "$a\t$r\t$rmse\t$lmin\n";
	`rm c rmse cor.R* glmnet1.R* cor act pred`;
	}
$alpha=`sort -nk3 perf_$i |head -1 |cut -f3`;chomp $alpha;
open(FF,">>glmnet2.R");
print FF "library(glmnet)\n";
print FF "x<-as.matrix(read.csv(\"x2\",sep=\"\\t\",header=FALSE))\n";
print FF "y<-as.matrix(read.csv(\"y2\",sep=\"\\t\",header=FALSE))\n";
print FF "fit = glmnet(x, y,lambda=$lmin,alpha=$alpha)\n";
print FF "d<-coef(fit,s=\"$lmin\")\n";
print FF "e<-as.matrix(d)\n\n";
print FF "\n";
print FF "write.table(e,file=\"f\")\n";
close FF;
`/scratch/breakthr/rkumar/R/bin/R CMD BATCH glmnet2.R`;
###### gene selection #####
        `perl -pi -e 's/ 0\$/ #/g' f`;
        `grep -v '#' f |cut -d ' ' -f1 |perl -pi -e 's/"//g' |perl -pi -e 's/V//g' |grep -v -w '^1' |grep -v 'Intercept' >g`;
        `grep -v '#' f |cut -d ' ' -f1,2 |perl -pi -e 's/"//g' |perl -pi -e 's/V//g' |grep -v -w '^1' |grep -v 'Intercept' |perl -pi -e 's/ /\t/g'>g`;
        open(FG,">>genes_elastic_$i");
        open(IN,"g");
        while($line=<IN>)
        {
        chomp $line;
        @sp=split('\t',$line);
        $gene=`sed '$sp[0]!d' genes`;chomp $gene;
        print FG "$gene\t$sp[1]\n";
        }
        close FG;
######################
`rm x y x1 y1 glmnet2.R* glmnet.R* f g train test x2 y2 genes lmin`;
print "\n";
}
