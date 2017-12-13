# adopted from Raymond's code
suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(Rsamtools));
suppressPackageStartupMessages(library(GenomicAlignments));
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(org.Hs.eg.db))

count <- read.table("geneCounts.txt", sep="\t", header=TRUE)

getExprs <- function( features, featureCounts, feature = 'gene' ){
        numBases <- sum(width(features))
        numKBases <- numBases / 1000
        millionsMapped <- sum(featureCounts) / 10^6
        rpm <- featureCounts / millionsMapped
        rpkm <- rpm / numKBases
        return( list( 'rpm' = rpm, 'rpkm' = rpkm) )
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txByGene <- transcriptsBy( txdb, 'gene' )
exonsByGene <- exonsBy(txdb, 'gene')
intronsByTx <- intronsByTranscript(txdb, use.names = TRUE )
introns <- unlist(intronsByTx)
introns <- introns[ !duplicated(introns) ]
newSeqNames <- sub('chr', '', seqlevels(txByGene))
names(newSeqNames) <- seqlevels(txByGene)
txByGene <- renameSeqlevels( txByGene, newSeqNames )

newSeqNames <- sub('chr', '', seqlevels(exonsByGene))
names(newSeqNames) <- seqlevels(exonsByGene)
exonsByGene <- renameSeqlevels( exonsByGene, newSeqNames )

newSeqNames <- sub('chr', '', seqlevels(introns))
names(newSeqNames) <- seqlevels(introns)
introns <- renameSeqlevels( introns, newSeqNames )

txnames <- names(introns)
map <- select(txdb, keys=txnames, keytype='TXNAME', columns='GENEID')
idx <- map$GENEID[!is.na(map$GENEID)]
intronsByGene <- split(introns[!is.na(map$GENEID)], idx)
names(intronsByGene) <- unique(idx)
genes <- c( names(txByGene), names(exonsByGene), names(intronsByGene) )
genes <- unique( genes )

exonExprsVals <- getExprs( exonsByGene, count[,2:4]) ## here provide number of columns with count
geneSymbols <- sapply(mget(genes, org.Hs.egSYMBOL, ifnotfound = NA), function (x) x[1])
x <- cbind(count$gene, exonExprsVals$rpkm)
write.table(x, file="RPKM.txt", sep="\t", row.names=FALSE, quote=FALSE)
