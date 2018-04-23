library(stringr)

list <- list.files("./", pattern="TCGA")

i=1
id <- substr(list[i], 1, 12)
tumor <- paste(id, "-01A.gatk_snps.vcf", sep="")
normal <- paste(id, "-10A.gatk_snps.vcf", sep="")
vcf_t <- read.table(file=tumor, sep="\t", header=F)
vcf_n <- read.table(file=normal, sep="\t", header=F)
colnames(vcf_n) <- c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Normal")
colnames(vcf_t) <- c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Tumor")
merge <- merge(vcf_t, vcf_n, by=c("CHROM", "POS", "REF", "ALT"))

## Tumor calculations
merge$new <- str_split_fixed(merge$Tumor, ":", 4)
merge$t_cov <- as.numeric(merge$new[,3])
merge$x <- str_split_fixed(merge$new[,2], ",", 2)
merge$t_ref <- as.numeric(merge$x[,1])
merge$t_alt <- as.numeric(merge$x[,2])
merge$taf <- as.numeric(merge$t_alt/merge$t_cov)

## Tumor calculations
merge$new <- str_split_fixed(merge$Normal, ":", 4)
merge$n_cov <- as.numeric(merge$new[,3])
merge$y <- str_split_fixed(merge$new[,2], ",", 2)
merge$n_ref <- as.numeric(merge$y[,1])
merge$n_alt <- as.numeric(merge$y[,2])
merge$naf <- as.numeric(merge$n_alt/merge$n_cov)

names <- c("CHROM", "POS", "REF", "ALT", "ID.x", "Tumor", "Normal", "t_ref", "t_alt", "taf", "n_cov", "n_ref", "n_alt", "naf")
merge <- merge[ , which(names(merge) %in% names)]
merge$Cancer <- "LUSC"
write.table(merge[1:100,], file="x.txt", sep="\t", quote=F, row.names=F)
