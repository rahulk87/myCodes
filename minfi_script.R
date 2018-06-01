library('minfi')
library('IlluminaHumanMethylationEPICmanifest')
library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
library('IlluminaHumanMethylation450kmanifest')
library('IlluminaHumanMethylation450kanno.ilmn12.hg19')

base = "../rawdata/"
targets = read.metharray.sheet(base)
RGSet = read.metharray.exp(base=base, targets=targets)
RGSet@annotation = c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b2.hg19")
gSet850 = preprocessFunnorm(RGSet)
snps = getSnpInfo(gSet850)
gSet850 = addSnpInfo(gSet850)
gSet850 = dropLociWithSnps(gSet850, snps=c("SBE","CpG"), maf=0)
beta850 = getBeta(gSet850)
phenoData = pData(gSet850)
colnames(beta850) = as.character(phenoData[,2])
annotation = getAnnotation(gSet850)
index = annotation[,"chr"]=="chrX" | annotation[,"chr"]=="chrY"
beta850 = beta850[!index,,drop=FALSE]
beta850_ann <- merge(beta850, annotation, by="row.names")
write.table(beta850_ann, file="beta850_ann.txt", sep="\t", quote=F, row.names=F)
