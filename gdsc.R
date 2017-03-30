ci_data <- read.csv("/Users/kumarr/work/common_files/common_things_icr/gdsc/drug_profiles/CI1040.txt",sep="\t",header=FALSE)

ci_data_sub <- subset(ci_data,ci_data$V2 != '')

col <- paste(ci_data_sub$V2, ci_data_sub$V3, sep = "_")

new_data <- cbind(ci_data_sub,col)

rownames(new_data) <- new_data$col

colnames(new_data) <- c("drug", "cell_line", "tissue", "cell_id", "drug_id", "AA", "logIC50", "AUC", "AB", "cell_tissue")

subset(new_data,new_data$cell_line=="REH")$logIC50

cells <- c("REH","SNU81", "OC314", "M14", "FARAGE", "PANC0403")
mut_ic <- new_data[grepl(paste(cells, collapse="|"), new_data$cell_line), ]$logIC50
wt_ic <- new_data[grepl(paste(cells, collapse="|"), new_data$cell_line), ]$logIC50
p <- t.test(wt_ic,mut_ic)$p.value


#cells <- "cells.txt"
#conn <- file(cells,open="r")
#linn <-readLines(conn)
#sink("mut_ic.txt")
#for (i in 1:length(linn)){
#	cell=linn[i]
	#print (cell)
	#print (subset(new_data,new_data$cell_line==cell)$logIC50)
#	ic = subset(new_data,new_data$cell_line==cell)$logIC50
#	cat(ic)
#	cat("\n")
#	}
#sink()
