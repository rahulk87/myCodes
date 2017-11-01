library("limma")
options(scipen=999)

## Read files ##

run6 <- read.table("../without_removing_batch/files/run6_log2.txt", sep="\t", header=TRUE, check.names=TRUE)
run7 <- read.table("../without_removing_batch/files/run7_log2.txt", sep="\t", header=TRUE, check.names=TRUE)
run8 <- read.table("../without_removing_batch/files/run8_log2.txt", sep="\t", header=TRUE, check.names=TRUE)
run9 <- read.table("../without_removing_batch/files/run9_log2.txt", sep="\t", header=TRUE, check.names=TRUE)

## remove batch effect ##
run_all <- cbind(run6[,c(4:15)], run7[,c(4:15)], run8[,c(4:9)], run9[,c(4:7)])
batch <- c("A","A", "A","A", "A","A", "A","A", "A","A", "A","A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "D", "D", "D", "D")
run_all_batch <- as.data.frame(removeBatchEffect(run_all, batch))
pdf("batch_correction.pdf")
par(mfrow=c(1,2))
boxplot(as.data.frame(run_all),main="Before")
boxplot(as.data.frame(run_all_batch),main="After")
dev.off()

## Group formation ##

groupA <- cbind(run6[,c(1:3)], run_all_batch$X11F1024, run_all_batch$X15F0915.A, run_all_batch$X10H2046.A, run_all_batch$X10H2046.B, run_all_batch$X09H1645)
groupB <- cbind(run6[,c(1:3)], run_all_batch$X10H2429, run_all_batch$X10H1078, run_all_batch$X12F0155B, run_all_batch$X11H0356)
groupC <- cbind(run6[,c(1:3)], run_all_batch$X12H1089, run_all_batch$X15F1048.A, run_all_batch$X15F1048.B, run_all_batch$X15F1036.A, run_all_batch$X15F1036.B, run_all_batch$X15F1051, run_all_batch$X10H0875, run_all_batch$X09H2044.1T1, run_all_batch$X09H2044.1FT.A, run_all_batch$X09H2044.1FT.B, run_all_batch$X15F1015)
groupD <- cbind(run6[,c(1:3)], run_all_batch$X11H2315, run_all_batch$X11H0264.1E1, run_all_batch$X11H0264.TVA)
groupE <- cbind(run6[,c(1:3)], run_all_batch$X10H1384, run_all_batch$X12H3134.2, run_all_batch$X12F0402, run_all_batch$X12F0403, run_all_batch$X12H1860, run_all_batch$X10H091, run_all_batch$X12H2273, run_all_batch$X12H2237, run_all_batch$X10H0362, run_all_batch$X12F0271, run_all_batch$X11H2784_27.30)
groupAB <- cbind(run6[,c(1:3)], run_all_batch$X11F1024, run_all_batch$X15F0915.A, run_all_batch$X10H2046.A, run_all_batch$X10H2046.B, run_all_batch$X09H1645, run_all_batch$X10H2429, run_all_batch$X10H1078, run_all_batch$X12F0155B, run_all_batch$X11H0356)

## Calculate group median##
# GroupA
grpA <- groupA[,4:8]
rownames(grpA) <- groupA$class
grpA_median <- as.data.frame(apply(grpA, 1, median))
colnames(grpA_median) <- "grpA_median"
#test <- cbind(grpA, grpA_median)
#write.table(test, "test.txt", sep="\t", quote=FALSE)
## Median plotting
grpA_median$genes <- rownames(grpA_median)
grpA_median <- grpA_median[!grepl("POS", grpA_median$gene),]
grpA_median <- grpA_median[!grepl("NEG", grpA_median$gene),]
x <- as.data.frame(grpA_median[order(grpA_median$grpA_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpA_median
names(my_vector)=y$genes
pdf(paste("grpA_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

# Group B
grpB <- groupB[,4:7]
rownames(grpB) <- groupB$class
grpB_median <- as.data.frame(apply(grpB, 1, median))
colnames(grpB_median) <- "grpB_median"
## Median plotting
grpB_median$genes <- rownames(grpB_median)
grpB_median <- grpB_median[!grepl("POS", grpB_median$gene),]
grpB_median <- grpB_median[!grepl("NEG", grpB_median$gene),]
x <- as.data.frame(grpB_median[order(grpB_median$grpB_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpB_median
names(my_vector)=y$genes
pdf(paste("grpB_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

# Group C
grpC <- groupC[,4:14]
rownames(grpC) <- groupC$class
grpC_median <- as.data.frame(apply(grpC, 1, median))
colnames(grpC_median) <- "grpC_median"
## Median plotting
grpC_median$genes <- rownames(grpC_median)
grpC_median <- grpC_median[!grepl("POS", grpC_median$gene),]
grpC_median <- grpC_median[!grepl("NEG", grpC_median$gene),]
x <- as.data.frame(grpC_median[order(grpC_median$grpC_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpC_median
names(my_vector)=y$genes
pdf(paste("grpC_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

# Group D
grpD <- groupD[,4:6]
rownames(grpD) <- groupD$class
grpD_median <- as.data.frame(apply(grpD, 1, median))
colnames(grpD_median) <- "grpD_median"
## Median plotting
grpD_median$genes <- rownames(grpD_median)
grpD_median <- grpD_median[!grepl("POS", grpD_median$gene),]
grpD_median <- grpD_median[!grepl("NEG", grpD_median$gene),]
x <- as.data.frame(grpD_median[order(grpD_median$grpD_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpD_median
names(my_vector)=y$genes
pdf(paste("grpD_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

# Group E
grpE <- groupE[,4:14]
rownames(grpE) <- groupE$class
grpE_median <- as.data.frame(apply(grpE, 1, median))
colnames(grpE_median) <- "grpE_median"
## Median plotting
grpE_median$genes <- rownames(grpE_median)
grpE_median <- grpE_median[!grepl("POS", grpE_median$gene),]
grpE_median <- grpE_median[!grepl("NEG", grpE_median$gene),]
x <- as.data.frame(grpE_median[order(grpE_median$grpE_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpE_median
names(my_vector)=y$genes
pdf(paste("grpE_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

# Group AB
grpAB <- groupAB[,4:12]
rownames(grpAB) <- groupAB$class
grpAB_median <- as.data.frame(apply(grpAB, 1, median))
colnames(grpAB_median) <- "grpAB_median"
## Median plotting
grpAB_median$genes <- rownames(grpAB_median)
grpAB_median <- grpAB_median[!grepl("POS", grpAB_median$gene),]
grpAB_median <- grpAB_median[!grepl("NEG", grpAB_median$gene),]
x <- as.data.frame(grpAB_median[order(grpAB_median$grpAB_median, decreasing=TRUE),])
y <- x[1:80,]
my_vector=y$grpAB_median
names(my_vector)=y$genes
pdf(paste("grpAB_median.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="log2 expression", las=2, col="black", cex.names=0.6, main="")
dev.off()

## limma ##
## for group A
## A_vs_C
this_analysis <- "A_vs_C"
R <- grpA_median$grpA_median
G <- grpC_median$grpC_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpA_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpA_median$grpA_median, grpC_median$grpC_median))
colnames(result) <- c("Gene", "Normalized_score", "Group A", "Group C")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## A_vs_D
this_analysis <- "A_vs_D"
R <- grpA_median$grpA_median
G <- grpD_median$grpD_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpA_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpA_median$grpA_median, grpD_median$grpD_median))
colnames(result) <- c("Gene", "Normalized_score", "Group A", "Group D")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## A_vs_E
this_analysis <- "A_vs_E"
R <- grpA_median$grpA_median
G <- grpE_median$grpE_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpA_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpA_median$grpA_median, grpE_median$grpE_median))
colnames(result) <- c("Gene", "Normalized_score", "Group A", "Group E")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## for group B
## B_vs_C
this_analysis <- "B_vs_C"
R <- grpB_median$grpB_median
G <- grpC_median$grpC_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpB_median$grpB_median, grpC_median$grpC_median))
colnames(result) <- c("Gene", "Normalized_score", "Group B", "Group C")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## B_vs_D
this_analysis <- "B_vs_D"
R <- grpB_median$grpB_median
G <- grpD_median$grpD_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpB_median$grpB_median, grpD_median$grpD_median))
colnames(result) <- c("Gene", "Normalized_score", "Group B", "Group D")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## B_vs_E
this_analysis <- "B_vs_E"
R <- grpB_median$grpB_median
G <- grpE_median$grpE_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpB_median$grpB_median, grpE_median$grpE_median))
colnames(result) <- c("Gene", "Normalized_score", "Group B", "Group E")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## A_vs_B
this_analysis <- "A_vs_B"
R <- grpA_median$grpA_median
G <- grpB_median$grpB_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpA_median$grpA_median, grpB_median$grpB_median))
colnames(result) <- c("Gene", "Normalized_score", "Group A", "Group B")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## AB_vs_C
this_analysis <- "AB_vs_C"
R <- grpAB_median$grpAB_median
G <- grpC_median$grpC_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpAB_median$grpAB_median, grpC_median$grpC_median))
colnames(result) <- c("Gene", "Normalized_score", "Group AB", "Group C")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## AB_vs_D
this_analysis <- "AB_vs_D"
R <- grpAB_median$grpAB_median
G <- grpD_median$grpD_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpAB_median$grpAB_median, grpD_median$grpD_median))
colnames(result) <- c("Gene", "Normalized_score", "Group AB", "Group D")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()

## AB_vs_E
this_analysis <- "AB_vs_E"
R <- grpAB_median$grpAB_median
G <- grpE_median$grpE_median
RG <- list(R, G)
names(RG) <- c("R", "G")
MA <- normalizeWithinArrays(RG, method="robustspline")
genes <- rownames(grpB_median)
result <- cbind(genes, MA$M)
result[is.na(result)] <- 0 ## replace all NAs to 0
result <- as.data.frame(cbind(result, grpAB_median$grpAB_median, grpE_median$grpE_median))
colnames(result) <- c("Gene", "Normalized_score", "Group AB", "Group E")
write.table(result, paste(this_analysis, ".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

## plotting
result <- result[!grepl("POS", result$Gene),]
result <- result[!grepl("NEG", result$Gene),]
result<-result[order(result$Normalized_score, decreasing=TRUE),]
a <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] < -1)
a <- a[order(a$Normalized_score, decreasing=FALSE),]
b <- subset(result, as.numeric(levels(result$Normalized_score))[result$Normalized_score] > 1)
b <- b[order(b$Normalized_score, decreasing=TRUE),]
c <- rbind(b, a)
my_vector=c$Normalized_score
my_vector=as.numeric(levels(my_vector))[my_vector]
names(my_vector)=c$Gene
pdf(paste(this_analysis, ".bar.pdf", sep=""), height=3, width=12)
barplot(my_vector, ylab="Normalized_score", las=2, col="black", cex.names=0.6, main=this_analysis, ylim=c(-6,6))
dev.off()
