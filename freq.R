mat <- read.csv("combined_maf.txt", sep="\t", header=TRUE)

genes<-read.table("../../../impact410_genes.txt", header=F)

sample_count <- length(as.vector(unique(mat$Tumor_Sample_Barcode)))

data <- list()

for(gene in 1:length(genes$V1)){
l <-  length(unique(subset(mat, mat$Hugo_Symbol==as.vector(genes$V1[gene]))$Tumor_Sample_Barcode))
freq <- l/sample_count
p <-paste(as.vector(genes$V1[gene]), l, freq, sep="\t")
data  <-  rbind(p, data)
}

write.table(data, file="freq.txt", quote=FALSE, row.names=FALSE)
