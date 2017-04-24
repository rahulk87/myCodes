library(DEXSeq)

## Make GFF file
#python ~/software/R-3.2.0/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
#/cbio/ski/reis-filho/home/kumarr/work/genomes/Homo_sapiens.GRCh37.75.gtf \
#$path/Homo_sapiens.GRCh37.75.gff

## Exon count
#python ~/software/R-3.2.0/lib64/R/library/DEXSeq/python_scripts/dexseq_count.py \
#-f bam \
#-p yes \
#$path/Homo_sapiens.GRCh37.75.gff \
#/cbio/ski/reis-filho/home/kumarr/work/tht/THT01T/THT01T.accepted_hits.rgfixed.reorder.bam \
#$path/THT01T_dex_my.txt

inDir = ("/Users/kumarr/work/Trabecular_adenoma_030317/Trabecular_Adenomas/dexseq/diff_analysis_run/")

countFiles = list.files(inDir, pattern="THT*", full.names=TRUE)

flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)

sampleTable = data.frame(row.names = c( "untreated1", "treated2", "untreated2", "treated3", "treated4", "treated5"), condition = c("control", "fuse", "control", "fuse", "fuse", "fuse"), libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))

dxd = DEXSeqDataSetFromHTSeq(
countFiles,
sampleData=sampleTable,
design= ~ sample + exon + condition:exon,
flattenedfile=flattenedFile )

genesForSubset = read.table(
file.path(inDir, "geneIDsinsubset.txt"),
stringsAsFactors=FALSE)[[1]]

## Make subset
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

## Accessing data
#colData(dxd)
#head( counts(dxd), 5 )
#head( featureCounts(dxd), 5 )
#head( rowRanges(dxd), 3 )
#sampleAnnotation( dxd )

## Normalisation (for different sequencng depth)
dxd = estimateSizeFactors( dxd )

## Estimate dispersion
dxd = estimateDispersions( dxd )

plotDispEsts( dxd )

## Differential analysis
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

## Save results
dxr1 = DEXSeqResults( dxd )

#mcols(dxr1)$description
table ( dxr1$padj < 0.1 )
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )

## plotting
plotMA( dxr1, cex=0.8 )

pdf ("dexse_glis3.pdf", width=15, height=5)
plotDEXSeq( dxr1, "ENSG00000107249", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, expression=FALSE, norCounts=TRUE)
dev.off()

pdf ("dexseq_pax8.pdf", width=15, height=5)
plotDEXSeq( dxr1, "ENSG00000125618", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, expression=FALSE, norCounts=TRUE)
dev.off()
