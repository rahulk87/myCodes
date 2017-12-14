## combat
library("sva")

mat <- read.table("../merge_file_rpkm_filtered.txt", header=TRUE, sep="\t")

mat1 <- mat[,2:187]

batch <- as.numeric(c(rep("1", 180), c(rep("2",6))))

combat_data = ComBat(dat=mat1, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

## heatmap
library(ComplexHeatmap)
library(circlize)
library(MASS)
library(pvclust)

combat_data.pv <- pvclust(combat_data, nboot=100, method.hclust="ward.D2",method.dist="euclidean")

df <- data.frame("Tumor Type" = c(rep("FA", 25), rep("FTC", 30), rep("FV", 48), rep("PTC", 43), rep("PTC_PAX8", 1), rep("PTC", 33), rep("THT", 6)))

colnames(df) <- "Tumor type"

ha = HeatmapAnnotation(df = df, col = list("Tumor type" = c("FA" = "red", "FTC" =  "blue", "FV" = "green", "PTC" = "darkgreen", "PTC_PAX8" = "darkmagenta", "THT" = "orange")), show_annotation_name = TRUE)

pdf("thyroid_samples_clutering.pdf", height=6, width=20)
Heatmap(combat_data, name = "RPKM", cluster_rows = FALSE, top_annotation = ha, show_row_names = FALSE, show_column_names = FALSE, column_dend_height = unit(5, "cm"), col=colorRamp2(c(200, 150, 100, 50, 0), c("brown4", "brown3", "brown2", "brown1", "white")), cluster_columns = combat_data.pv$hclust)
dev.off()

## pvclust plotting
pdf("pvclust0.5.pdf", height=8, width=20)
plot(combat_data.pv, cex=0.5)
pvrect(mat1.pv, alpha=.95)
dev.off()

pdf("pvclust0.8.pdf", height=8, width=20)
plot(combat_data.pv, cex=0.8)
pvrect(mat1.pv, alpha=.95)
dev.off()

pdf("pvclust1.pdf", height=8, width=20)
plot(combat_data.pv, cex=1)
pvrect(mat1.pv, alpha=.95)
dev.off()

## PCA
library(FactoMineR)
combat_data <- t(combat_data)
result <- PCA(combat_data, graph=FALSE)
pdf("facto_plot.pdf")
plot(result)
dev.off()
