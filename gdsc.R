###
# drugs.txt
# 17AAG.txt
# 5FLUOROURACIL.txt
###

drugs <- "drugs.txt"
conn <- file(drugs,open="r")
linn <-readLines(conn)
sink("drug_pValue.txt")
pdf("box.pdf")
	for (i in 1:length(linn)){
	drug = linn[i]
	ci_data <- read.csv(paste("/Users/kumarr/work/common_files/common_things_icr/gdsc/drug_profiles/", drug,sep=""),sep="\t",header=FALSE)
	ci_data_sub <- subset(ci_data,ci_data$V2 != '')
	col <- paste(ci_data_sub$V2, ci_data_sub$V3, sep = "_")
	new_data <- cbind(ci_data_sub,col)
	rownames(new_data) <- new_data$col
	colnames(new_data) <- c("drug", "cell_line", "tissue", "cell_id", "drug_id", "AA", "logIC50", "AUC", "AB", "cell_tissue")
	subset(new_data,new_data$cell_line=="REH")$logIC50
	cells <- c("REH","SNU81", "OC314", "M14", "FARAGE", "PANC0403")
	mut_ic <- new_data[grepl(paste(cells, collapse="|"), new_data$cell_line), ]$logIC50
	wt_ic <- new_data[!grepl(paste(cells, collapse="|"), new_data$cell_line), ]$logIC50
	l <- length(mut_ic)
	if(l > 3){
		cat (drug)
		cat ("\t")
		p <- t.test(wt_ic,mut_ic)$p.value
		p <- round(p, digits = 4)
		cat (p)
		cat ("\n")
		n <- max(length(mut_ic), length(wt_ic))
		length(mut_ic) <- n 
		length(wt_ic) <- n
		profile<-as.data.frame(cbind(wt_ic,mut_ic))
		drug <- lapply(drug, function(x) {sub(".txt", "", x)})
		#pdf(paste(drug,".pdf",sep="")) ## For individual plots
		boxplot(profile,main=paste(drug, "P:", p, sep=" "),names=c("WT","Mutant"))
		stripchart(profile,vertical=TRUE,method = "jitter", pch = 19,add = TRUE,col=c("blue","red"))
		#dev.off()
	}
}
sink()
dev.off()
