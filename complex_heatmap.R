library(ComplexHeatmap)
library(circlize)

pheno <- read.table("pheno.txt", header = TRUE, sep = "\t")

colnames(pheno) <- c("ID", "Tumor Stage", "Histologic Grade", "Histologic Type", "Associated Borderline Tumor")

## Defining color annotations ##
#ha <- HeatmapAnnotation(df = pheno[,2:5], col = list(stage = c("IA" =  "red", "IC" = "blue", "IIB" = "yellow", "IIC" = "green", "IIIA" = "burlywood4", "IIIB" = "orange", "IIIC" = "pink"), grade = c("3" = "cadetblue4", "2" = "cadetblue3", "1" = "cadetblue1"), type = c("Intestinal" = "darkgoldenrod3", "Simple" = "cyan4", "Gastric" = "coral4", "Unknown" = "grey"), border = c("Yes" = "darkolivegreen", "No" = "darkmagenta", "Unknown" = "grey")), na_col = "grey", gap = unit(2, "mm"))

ha <- HeatmapAnnotation(df = pheno[,2:5], col = list("Tumor Stage" = c("IA" =  "red", "IC" = "blue", "IIB" = "yellow", "IIC" = "green", "IIIA" = "burlywood4", "IIIB" = "orange", "IIIC" = "pink"), "Histologic Grade" = c("3" = "cadetblue4", "2" = "cadetblue3", "1" = "cadetblue1"), "Histologic Type" = c("Intestinal" = "darkgoldenrod3", "Simple" = "cyan4", "Gastric" = "coral4", "Unknown" = "grey"), "Associated Borderline Tumor" = c("Yes" = "darkolivegreen", "No" = "darkmagenta", "Unknown" = "grey")), na_col = "grey", gap = unit(2, "mm"), show_annotation_name = TRUE)

## Creating dummy matrix of equal number of columns
mat = matrix(rnorm(120, 2), 8, 15)
rownames(mat) = paste0("R", 1:8)
colnames(mat) = paste0("C", 1:15)
####################################

## Creating null matrix###################
zero_row_mat = matrix(nrow = 0, ncol = 15)
##########################################

pdf("pheno.pdf", height = 8, width = 10)
#png("pheno.png", width = 4000, height = 4000, res = 400)
ht = Heatmap(mat, top_annotation = ha, cluster_columns = FALSE, show_heatmap_legend = FALSE)
draw(ht)
dev.off()


pdf("pheno_only_pheno.pdf", height = 8, width = 10)
#png("pheno.png", width = 4000, height = 4000, res = 400)
ht1 = Heatmap(zero_row_mat, top_annotation = ha, cluster_columns = FALSE, show_heatmap_legend = FALSE)
draw(ht1, show_annotation_legend = FALSE)
dev.off()
