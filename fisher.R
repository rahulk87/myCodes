impact <- read.table("../impact/freqImpact_samples_impact341.txt", sep = "\t", header= FALSE, stringsAsFactors = FALSE)
tcga <- read.table("../tcga_mucinous_colon_fireBrowse/freq341.txt", sep = "\t", header= FALSE, stringsAsFactors = FALSE)
impact1 <- cbind(impact, (15 - impact[,2]))
tcga1 <- cbind(tcga, (32 - tcga[,2]))
colnames(impact1) <- c("gene", "count", "freq", "diff")
colnames(tcga1) <- c("gene", "count", "freq", "diff")

sink("output.txt")

for (i in 1:length(tcga1$gene)){
	gene <- tcga1[i, "gene"]
	count_impact <- impact1[i, "count"]
	freq_impact <- impact1[i, "freq"]
	count_tcga <- tcga1[i, "count"]
        freq_tcga <- tcga1[i, "freq"]
	mat <- matrix(c(impact1[i, "count"], impact1[i, "diff"], tcga1[i, "count"], tcga1[i, "diff"]), nrow = 2)
	p <- fisher.test(mat)$p.value
	cat (gene)
	cat ("\t")
	cat (count_impact)
	cat ("\t")
	cat (freq_impact)
	cat ("\t")
	cat (count_tcga)
	cat ("\t")
	cat (freq_tcga)
	cat ("\t")
	cat (p)
	cat ("\n")
	#print(p)
	}
sink()

out <- read.table("output.txt", sep = "\t", header=FALSE)
colnames(out) <- c("gene", "count_impact", "freq_impact", "count_tcga", "freq_tcga", "pvalue")
out$p.adjust <- p.adjust(out$pvalue, method = "BH")
write.table(out, file = "result_fisher.txt", sep = "\t")
system("rm output.txt")
