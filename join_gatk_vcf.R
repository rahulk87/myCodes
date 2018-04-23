library(stringr)

list <- list.files("../vcf_ann/", pattern="snps")
list <- unique(substr(list, 1, 12))

for(i in 1:length(list)){

id <- list[i]
samples <- list.files("../vcf_ann/", pattern=id)
t_code <- substr(samples[1], 13,16)
n_code <- substr(samples[3], 13,16)

#print(i)
print(id)
#print(n_code)
#print(t_code)

tumor <- paste(id, t_code, ".gatk_snps.vcf", sep="")
normal <- paste(id, n_code, ".gatk_snps.vcf", sep="")
vcf_t <- read.table(file=paste("../vcf_ann/", tumor, sep=""), sep="\t", header=F)
vcf_n <- read.table(file=paste("../vcf_ann/", normal, sep=""), sep="\t", header=F)
colnames(vcf_n) <- c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Normal")
colnames(vcf_t) <- c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Tumor")
merge <- merge(vcf_t, vcf_n, by=c("CHROM", "POS", "REF", "ALT"))

## Tumor calculations
merge$new <- str_split_fixed(merge$Tumor, ":", 4)
merge$tcov <- as.numeric(merge$new[,3])
merge$x <- str_split_fixed(merge$new[,2], ",", 2)
merge$trc <- as.numeric(merge$x[,1])
merge$tac <- as.numeric(merge$x[,2])
merge$taf <- as.numeric(merge$tac/merge$tcov)

## Tumor calculations
merge$new <- str_split_fixed(merge$Normal, ":", 4)
merge$ncov <- as.numeric(merge$new[,3])
merge$y <- str_split_fixed(merge$new[,2], ",", 2)
merge$nrc <- as.numeric(merge$y[,1])
merge$nac <- as.numeric(merge$y[,2])
merge$naf <- as.numeric(merge$nac/merge$ncov)

## other required columns
merge$y <- str_split_fixed(merge$INFO.y, '\\|', 31)
merge$Variant_Classification <- merge$y[,2]
merge$Effect <- merge$y[,3]
merge$Hugo_Symbol <- merge$y[,4]
merge$cDNA <- merge$y[,10]
merge$AAChange <- merge$y[,11]

names <- c("CHROM", "POS", "REF", "ALT", "ID.x", "Tumor", "Normal", "tcov", "trc", "tac", "taf", "ncov", "nrc", "nac", "naf", "Variant_Classification", "Effect", "Hugo_Symbol", "cDNA", "AAChange")
merge <- merge[ , which(names(merge) %in% names)]
merge$Cancer <- "LUSC"
write.table(merge, file=paste(id, ".vcf", sep=""), sep="\t", quote=F, row.names=F)
}
