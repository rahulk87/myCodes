mat<-read.table("ct83_expression.txt", sep="\t", header=TRUE)

pdf("CT83_expression_rpkm.pdf")
boxplot(mat$exonRPKM, col="grey", ylab="RPKM", main="CT83 expression in TCGA (n=107)")
stripchart(mat$exonRPKM, add=TRUE, method="jitter", vertical = TRUE, pch=19, col="black")
dev.off()

pdf("CT83_expression_raw.pdf")
boxplot(mat$countsByGene, col="green", ylab="Raw counts", main="CT83 expression in TCGA (n=107)")
stripchart(mat$exonRPKM, add=TRUE, method="jitter", vertical = TRUE, pch=19, col="darkgreen")
dev.off()
