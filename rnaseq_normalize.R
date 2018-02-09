library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(RColorBrewer)

colors <- brewer.pal(3, "Set2")

mat1 <- read.table("/Users/kumarr/work/sclerosing_stromal_7008/rnaseq/geneCounts.txt", sep="\t", header=TRUE)
rownames(mat1) <- paste(rownames(mat1),mat1$gene, sep="_")
mat1 <- mat1[,2:7]
filter <- apply(mat1, 1, function(x) length(x[x>5])>=2)
filtered <- mat1[filter,]
filtered <- data.matrix(filtered)
print(dim(filtered))

## Plot unnormalized data
pdf("Unnormalized_plots.pdf")
plotRLE(filtered, outline=FALSE, ylim=c(-4, 4), col=colors, ylab = "Relative Log Expression")
plotPCA(filtered, col=colors, cex=1.2)
dev.off()

## Between lane normalization
filtered_norm <- betweenLaneNormalization(filtered, which="upper")
pdf("Normalized_plots.pdf")
plotRLE(filtered_norm, outline=FALSE, ylim=c(-4, 4), col=colors, ylab = "Relative Log Expression")
plotPCA(filtered_norm, col=colors, cex=1.2)
dev.off()

## RUVg normalization using HK genes
## HK gene sets
spikes <- filtered[1:100,] ## change here for using different set of control genes
spikes <- rownames(spikes)
##########################

filtered_ruvg <- RUVg(filtered_norm, spikes, k=1) # giving output of "Between lane normalization" here.
pdf("RUVg_normalized_plots.pdf")
plotRLE(filtered_ruvg$normalizedCounts, outline=FALSE, ylim=c(-4, 4), col=colors, ylab = "Relative Log Expression")
plotPCA(filtered_ruvg$normalizedCounts, col=colors, cex=1.2)
dev.off()

## Write normalize counts
print(dim(filtered_ruvg$normalizedCounts))
filtered_ruvg <- cbind(rownames(filtered), filtered_ruvg$normalizedCounts)
write.table(filtered_ruvg, file="RUVg_normalized_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)
