pdf("barplot.pdf", width=6, height=4)
mat <- read.table("irs4_rna_exp.txt", sep="\t", header=TRUE)
my_vector=mat$exonRPKM
names(my_vector)=mat$Sample
barplot(my_vector, ylab="RPKM", las=2, col="black", ylim=c(0,300))
dev.off()
