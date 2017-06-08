#!/usr/bin/perl
open(IN, "list");
	while($line=<IN>)
	{
	chomp $line;print "Running for $line....\n";
	@sp=split("_", $line);
	`grep -v '#' vcf/$line |grep -v 'REJECT' |cut -f1,2,4,5 >~/Downloads/tmp.txt`;
	`perl -pi -e 's/^/1\tchr/g' ~/Downloads/tmp.txt`;
	
	open(FA,">>decons.R");
	print FA "library(deconstructSigs)\n";
	print FA "mat <- read.csv(\"~/Downloads/tmp.txt\", header=TRUE, sep=\"\\t\")\n";
	print FA "colnames(mat) <- c(\"Sample\", \"chr\", \"pos\", \"ref\", \"alt\")\n";
	print FA "sigs.input <- mut.to.sigs.input(mut.ref = mat, sample.id = \"Sample\", chr = \"chr\",pos = \"pos\", ref = \"ref\",alt = \"alt\")\n";
	print FA "sample_1 <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = 1,contexts.needed = TRUE,tri.counts.method = 'default')\n";
	print FA "pdf(\"$sp[0].pdf\")\n";
	print FA "plotSignatures(sample_1, sub = '$sp[0]')\n";
	print FA "makePie(sample_1, sub = '$sp[0]')\n";
	print FA "dev.off()\n";
	close FA;
	`R CMD BATCH decons.R`;
	`rm decons.R* ~/Downloads/tmp.txt`;
	}
