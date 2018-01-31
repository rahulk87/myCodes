## All biallelic cases
library(data.table)
library(ggplot2)

load("combined_snps_exac_LOH.RData")
combined_snps_exac_LOH$aa_gene <- paste(combined_snps_exac_LOH$Entrez_Gene_Id, combined_snps_exac_LOH$AAChange, sep="_")

master_score <- fread("../../Master_Score.txt")

snp14<-read.table("snp_14.txt", sep="\t", header=FALSE)
rahul=fread("../../more_signatures.txt")
    rahul$Id <- substr(rahul$Id, 1, 12)
    rahul <- rahul[!duplicated(rahul$Id),]
    setkey(master_score, 'Id')
    setkey(rahul, 'Id')
    overlapping_Ids <- master_score$Id[master_score$Id %in% rahul$Id ]

    replace_columns  <- c("dominant", "number.of.mutations")
    for (nn in seq(1,30)) {
        strname=paste("Signature.", nn, sep="")
        replace_columns <- c(replace_columns, strname)
    }
    for (nn in replace_columns) {
        master_score[overlapping_Ids][[nn]] <- rahul[overlapping_Ids][[nn]]
    }

combined_snps_biallelic <- merge(combined_snps_exac_LOH, master_score, by.x="Sample", by.y="Id", all.x=TRUE)

pdf("all_biallelic.pdf", height=8, width=8)
par(mfrow=c(2,2))
for(i in 1:length(snp14$V1)){
  x <- subset(combined_snps_biallelic, combined_snps_biallelic$aa_gene==as.vector(snp14$V1[i]))
  x$Signature.3[is.na(x$Signature.3)] <- 0
  with(x, plot(Signature.3, LST, pch=20, main=paste(as.vector(snp14$V1[i]), ", n=", nrow(x), sep="")))
  with(subset(x, is.na(x$dominant)==TRUE), points(Signature.3, LST, pch=20, col="orange"))
  }
dev.off()

## GGplot
pdf("all_biallelic_ggplot.pdf")
for(i in 1:length(snp14$V1)){
  x <- subset(combined_snps_biallelic, combined_snps_biallelic$aa_gene==as.vector(snp14$V1[i]))
 x$Signature.3[is.na(x$Signature.3)] <- 0
 x$Signatures <- ifelse(is.na(x$dominant)==TRUE, "Absent", "Present")
 plot <- ggplot(x, aes(Signature.3, LST, size=Taf, colour=Signatures)) + geom_point() + scale_colour_manual(values = c("Blue", "Red")) + labs(title = paste(as.vector(snp14$V1[i]), ", n=", nrow(x), sep="")) + labs(size = "Tumor VAF") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))
 print(plot)
#do.call("grid.arrange", plot)
}
dev.off()


## Either high LST or dominant sig 3
load("../../second_phase/or/final_OR_after_exclusion.RData")

finalOR$aa_gene <- paste(finalOR$Entrez_Gene_Id, finalOR$AAChange, sep="_")

pdf("hig_LST_or_dom_sig3.pdf", height=8, width=8)
par(mfrow=c(2,2))
for(i in 1:length(snp14$V1)){
  x <- subset(finalOR, finalOR$aa_gene==as.vector(snp14$V1[i]))
  x$Signature.3[is.na(x$Signature.3)] <- 0
  with(x, plot(Signature.3, LST, pch=20, main=paste(as.vector(snp14$V1[i]), ", n=", nrow(x), sep="")))
  with(subset(x, is.na(x$dominant)==TRUE), points(Signature.3, LST, pch=20, col="orange"))
  }
dev.off()

## GGplot
pdf("high_LST_or_dom_sig3_ggplot.pdf")
for(i in 1:length(snp14$V1)){
x <- subset(finalOR, finalOR$aa_gene==as.vector(snp14$V1[i]))
x$Signature.3[is.na(x$Signature.3)] <- 0
x$Signatures <- ifelse(is.na(x$dominant)==TRUE, "Absent", "Present")
plot <- ggplot(x, aes(Signature.3, LST, size=Taf, colour=Signatures)) + geom_point() + scale_colour_manual(values = c("Blue", "Red")) + labs(title = paste(as.vector(snp14$V1[i]), ", n=", nrow(x), sep="")) + labs(size = "Tumor VAF") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))
print(plot)
#do.call("grid.arrange", plot)
}
dev.off()
