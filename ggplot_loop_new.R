rm(list=ls())
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

####
load("final_OR_after_exclusion.RData")
#write.table(finalOR, file="finalOR.txt", sep="\t", row.names=FALSE)
finalOR$aa_gene <- with(finalOR, paste(Hugo_Symbol, AAChange, sep="_"))
x <- as.data.frame(table(finalOR$aa_gene))
x <- x[order(x$Freq, decreasing=TRUE),]
snps <- subset(x, x$Freq > 1)$Var1
#snps <- unique(finalOR$aa_gene)
####

load("../combined_snps_exac_LOH.RData")
combined_snps_exac_LOH$aa_gene <- with(combined_snps_exac_LOH, paste(Hugo_Symbol, AAChange, sep="_"))

master_score <- fread("../../../../Master_Score.txt")

rahul=fread("../../../../more_signatures.txt")
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

combined_path_monoallelic <- merge(combined_snps_exac_LOH, master_score, by.x="Sample", by.y="Id", all.x=TRUE)

combined_path_monoallelic$X <- ifelse(combined_path_monoallelic$LOH_sample_counts > 0, "X", 0) ## for plotting purpose
write.table(combined_path_monoallelic, file="combined_path_monoallelic.txt", sep="\t", row.names=FALSE)

## All samples
grid.arrange(rectGrob(), rectGrob())
pl <- lapply(1:length(snps), function(i) {
  x <- subset(combined_path_monoallelic, combined_path_monoallelic$aa_gene==snps[i])
  x$Signature.3[is.na(x$Signature.3)] <- 0
  x$Signatures <- ifelse(is.na(x$dominant)==TRUE, "Absent", "Present")
  cols <- c("Present" = "Blue", "Absent" = "Red")
  return(ggplot(x, aes(Signature.3, LST, size=Taf)) +  geom_point(aes(colour = Signatures), alpha=1/1.5) + scale_colour_manual(values = cols) + labs(title = paste(snps[i], ", n=", nrow(x), sep="")) + labs(size = "Tumor VAF") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"))+ theme(plot.title = element_text(size=8)) + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) + guides(size = guide_legend(order = 1)))
    }
  )

ml <- marrangeGrob(pl, nrow=1, ncol=2)
ggsave("monoallelic_path_ggplot.pdf", ml)

## Either high LST or dominant sig 3
finalOR$X <- ifelse(finalOR$biallelic_N > 0, "X", 0) ## for plotting purpose
write.table(finalOR, file="finalOR.txt", sep="\t", row.names=FALSE)

grid.arrange(rectGrob(), rectGrob())
pl <- lapply(1:length(snps), function(i) {
  x <- subset(finalOR, finalOR$aa_gene==snps[i])
  x$Signature.3[is.na(x$Signature.3)] <- 0
  x$Signatures <- ifelse(is.na(x$dominant)==TRUE, "Absent", "Present")
  cols <- c("Present" = "Blue", "Absent" = "Red")
  return(ggplot(x, aes(Signature.3, LST, size=Taf)) +  geom_point(aes(colour = Signatures), alpha=1/1.5) + scale_colour_manual(values = cols) + labs(title = paste(snps[i], ", n=", nrow(x), sep="")) + labs(size = "Tumor VAF") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + theme(plot.title = element_text(size=8)) + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) + guides(size = guide_legend(order = 1)))
    }
  )

ml <- marrangeGrob(pl, nrow=1, ncol=2)
ggsave("high_LST_or_dom_sig3_monoallelic_path_ggplot.pdf", ml)
