## This will take the results file that was created using the combineMAF code and create 2 mutation files
setwd('./')
results      = read.csv("mucinous_ovarian_cancer_results_for_fig.csv", header=T, stringsAsFactors = F, sep="\t")

# If you are doing a sample set where there are multiple components you can use this to set the overall samples to get individual plots
# Otherwise this is just set to make all samples from one patient to make one plot.
results$Patient = "A" #substr(results$NORMAL_SAMPLE, 1, nchar(results$NORMAL_SAMPLE)-1)
patients  = unique(results$Patient)

samples = unique(results$TUMOR_SAMPLE)
#### If not all of the samples have mutations then this won't add those samples into the figure.
# Alternative method of getting sample list: 
# files = list.files(pattern = "\\cncf.txt$")
# samples = c(as.data.frame(strsplit(files, "_"))[1,])
# samples = as.vector(unname(samples))
# samples = as.character(unlist(samples))
b=1
#for (b in 1:length(patients)){
  muts = results[results$Patient == patients[b],]
  muts = muts[!is.na(muts$Cancer_Cell_Fraction),]
plot_file     = paste('Mutation_Heatmap', patients[b], ".pdf", sep="")
plot_file_CCF = paste('CCF_Heatmap', patients[b], ".pdf", sep="")
########################################################
muts$CCF_group = "0"
for (i in 1:nrow(muts)){
  if(muts$Cancer_Cell_Fraction[i] <= 1 & muts$Cancer_Cell_Fraction[i] > 0.8){muts$CCF_group[i] = 5}
  if(muts$Cancer_Cell_Fraction[i] <= 0.8 & muts$Cancer_Cell_Fraction[i] > 0.6){muts$CCF_group[i] = 4}
  if(muts$Cancer_Cell_Fraction[i] <= 0.6 & muts$Cancer_Cell_Fraction[i] > 0.4){muts$CCF_group[i] = 3}
  if(muts$Cancer_Cell_Fraction[i] <= 0.4 & muts$Cancer_Cell_Fraction[i] > 0.2){muts$CCF_group[i] = 2}
  if(muts$Cancer_Cell_Fraction[i] <= 0.2 & muts$Cancer_Cell_Fraction[i] > 0.05){muts$CCF_group[i] = 1}
}


muts                 = muts[order(muts$SYMBOL, decreasing=F),]
muts                 = muts[!is.na(muts$Cancer_Cell_Fraction),]
muts                 = muts[order(muts$CCF_group, decreasing=T),]
muts                 = muts[order(muts$TUMOR_SAMPLE, decreasing=F),]

###### Do any necessary filtering of mutation dataframe ########
#muts = subset(muts, subset = muts$patient == "RG6T") #| muts$TUMOR_SAMPLE == "AM2" | muts$TUMOR_SAMPLE == "AM3" | muts$TUMOR_SAMPLE == "AM4" | muts$TUMOR_SAMPLE == "AM5 - Primary_AME" | muts$TUMOR_SAMPLE == "AM6" | muts$TUMOR_SAMPLE == "AM7" | muts$TUMOR_SAMPLE == "AM32AME1")# & muts$TUMOR_SAMPLE != "AM8-Ips-Breast-Rec-Carc" & muts$TUMOR_SAMPLE != "AM5 - Axillary_Lymph_Node_Metastasis" & muts$TUMOR_SAMPLE != "AM46T2" & muts$TUMOR_SAMPLE != "AM5 - Primary_Carcarcinoma" & muts$TUMOR_SAMPLE != "AM8-LNE" & muts$TUMOR_SAMPLE != "AM8-LNM")
muts = subset(muts, subset = muts$Variant_Classification != "synonymous_variant")# & muts$ANN....EFFECT != "splice_region_variant&synonymous_variant" & muts$ANN....EFFECT != "downstream_ANN....GENE_variant" 
#           & muts$ANN....EFFECT != "intron_variant" & muts$ANN....EFFECT != "frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant"
#& muts$ANN....EFFECT != "non_coding_exon_variant|synonymous_variant" & muts$ANN....EFFECT != "SYNONYMOUS_CODING" & muts$ANN....EFFECT != "splice_acceptor_variant&splice_region_variant&intron_variant" & muts$ANN....EFFECT != "Silent")
#muts = subset(muts, subset = muts$ANN....EFFECT != "downstream_gene_variant|synonymous_variant" & muts$ANN....EFFECT != "upstream_gene_variant|5_prime_UTR_variant")
#muts = subset(muts, subset = muts$lawrence=="TRUE" | muts$kandoth =="TRUE" | muts$cancer_gene_census == "TRUE")

#########
muts$CGC = ''
for (i in 1:nrow(muts)) {
  if (muts$cancer_gene_census[i] == "TRUE" || muts$kandoth[i] == "TRUE" || muts$lawrence[i] == "TRUE") { muts$CGC[i] = '*'
  }
}
muts$newGene <- paste(muts$CGC, muts$SYMBOL, sep="")
#########

sample_names   =  as.list(sort(unique(muts$TUMOR_SAMPLE))) 
if(length(sample_names) == 1){ sample_names = c(sample_names, "Spacer")}
mutation_genes = unique(muts$newGene)
if(length(mutation_genes) == 1){ mutation_genes = c(mutation_genes, "Spacer")}
rownames(muts) = 1:nrow(muts)

###########################################################################
#impact410_list = read.csv('/Users/selenicp/Documents/Code/IMPACT410_genes.csv', header = T, stringsAsFactors = F)
# impact341 = read.csv('/Users/burkek/Documents/Code/IMPACT341_genes.csv', header =T, stringsAsFactors = F)
#impact410_list = unique(c(impact410_list$Approved.Symbol, impact410_list$HGNC.symbol))
#mutation_genes = unique(c(intersect(mutation_genes,impact410_list)))
#muts = muts[muts$ANN....GENE %in% mutation_genes,]
###########################################################################

return_gene_order   = T
sort_samples        = T
sort_genes          = T
show_sample_names   = T
TCGA                = F
remove_genes_with_no_mutation = F
width               = NULL 
height              = NULL
include_percentages = T
sample_name_col     = "TUMOR_SAMPLE"

#### Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples

mutation_heatmap <- matrix(9, nrow=sum(unlist(lapply(sample_names, length))), ncol=sum(unlist(lapply(mutation_genes, length))))
rownames(mutation_heatmap) <- unlist(sample_names)
colnames(mutation_heatmap) <- paste(mutation_genes)

#### Make sure the sample and mutations are both in the list of gene mutations and gene samples
if (!TCGA) { smallmaf <- muts[which(muts$newGene %in% unlist(mutation_genes) & muts$TUMOR_SAMPLE %in% unlist(sample_names)),]
}else { 
  muts$id <- unlist(lapply(muts$TUMOR_SAMPLE, function(x){substr(x, 1, 12)}))
  print(head(muts$id))
  print(head(unlist(sample_names)))
  smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),] 
}

### Define categories for different Effects
cat1 <- c('missense_variant_hotspot', 'stop_gained_hotspot','missense_variant&splice_region_variant_hotspot')
cat2 <- c("STOP_GAINED", "Nonsense_Mutation", "stop_gained&splice_region_variant", "stop_gained", "stop_gained&inframe_insertion", "stop_gained&disruptive_inframe_insertion", "stop_gained&disruptive_inframe_insertion&splice_region_variant", "stop_gained&inframe_insertion|stop_gained&inframe_insertion", "stop_gained&inframe_insertion&splice_region_variant", "stop_gained|stop_gained"  )
cat3 <- c("frameshift_indel","FRAME_SHIFT", "FRAME_SHIFT", "frameshift_variant&stop_lost", "frameshift_variant&stop_lost|frameshift_variant&stop_lost", "frameshift_variant&start_lost","frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant", "Frame_Shift_Del", "Frame_Shift_Ins", "frameshift_variant", "frameshift_variant|frameshift_variant",  "frameshift_variant&stop_gained", "frameshift_variant&splice_region_variant", "frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&stop_gained&splice_region_variant", "frameshift_variant&stop_gained|frameshift_variant&stop_gained","frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant" )
cat4 <- c("missense_snv_recurrent","missense_snv","NON_SYNONYMOUS_CODING", "Missense_Mutation", "missense_variant", "missense_variant&splice_region_variant", "missense_variant|missense_variant", "missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant", "missense_variant|missense_variant|missense_variant", "missense_variant&splice_region_variant|missense_variant&splice_region_variant","missense_variant|protein_protein_contact|protein_protein_contact|protein_protein_contact|protein_protein_contact|protein_protein_contact|protein_protein_contact|protein_protein_contact")
cat5 <- c("inframe_indel_recurrent","inframe_indel","CODON_CHANGE_PLUS_CODON_DELETION", "CODON_DELETION", "CODON_INSERTION", "In_Frame_Ins", "In_Frame_Del", "inframe_deletion|inframe_deletion", "disruptive_inframe_deletion", "disruptive_inframe_insertion", "inframe_deletion", "inframe_insertion", "disruptive_inframe_deletion&splice_region_variant", "inframe_deletion&splice_region_variant", "inframe_insertion&splice_region_variant", "disruptive_inframe_insertion|disruptive_inframe_insertion", "disruptive_inframe_insertion&splice_region_variant", "inframe_insertion|inframe_insertion", "disruptive_inframe_deletion|disruptive_inframe_deletion" )
cat6 <- c("splice_donor_variant&intron_variant|splice_donor_variant&intron_variant","splice_acceptor_variant&intron_variant|splice_acceptor_variant&intron_variant","splice_acceptor_variant&splice_donor_variant&intron_variant","splice_mut","SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_REGION", "Splice_Site", "splice_donor_variant&intron_variant", "splice_acceptor_variant&intron_variant", "splicing", "splice_donor_variant&splice_region_variant&intron_variant", "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "splice_donor_variant&inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "splice_acceptor_variant&inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "Splice_Region", "splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|start_lost&inframe_deletion&splice_region_variant", "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant","splice_acceptor_variant&splice_region_variant&intron_variant&non_coding_exon_variant", "splice_acceptor_variant&splice_region_variant&intron_variant")
cat7 <- c("STOP_LOST", "START_LOST", "START_GAINED", "UTR_5_PRIME", "start_lost", "stop_lost", "start_lost&inframe_deletion", "stop_lost&splice_region_variant", "stop_lost&disruptive_inframe_insertion")
cat8 <- c("upstream_gene_variant","downstream_gene_variant|synonymous_variant", "upstream_gene_variant|5_prime_UTR_variant","upstream_gene_variant|synonymous_variant")
cat9 <- c("nonsense_snv","synonymous_variant","Silent", "splice_region_variant&synonymous_variant", "downstream_gene_variant", "intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant",
         "non_coding_exon_variant|synonymous_variant", "SYNONYMOUS_CODING", "synonymous_variant|synonymous_variant", "splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant", "intergenic_region", "intron_variant","intron_variant|downstream_gene_variant","intron_variant|intron_variant","intergenic_region|downstream_gene_variant","intron_variant|upstream_gene_variant")


#### For each row read the Effect and create the type based on which category it fits in
#### If there is an error because the mutation type is unknow just add it to the correct category above
smallmaf$mut_type = '0'
smallmaf$types = '0'
for (i in 1:nrow(smallmaf)) {
  if(!TCGA) { type = smallmaf$Variant_Classification[i] } else { type = smallmaf$Variant_Classification[i] }
  if (smallmaf$HOTSPOT[i] == "TRUE" || type %in% cat1) { type = 1; smallmaf$mut_type[i] = 'A'; smallmaf$types[i] = 1
  } else if (type %in% cat2) { type = 2; smallmaf$mut_type[i] = 'B'; smallmaf$types[i] = 2
  } else if (type %in% cat3) { type = 3; smallmaf$mut_type[i] = 'C'; smallmaf$types[i] = 3
  } else if (type %in% cat4) { type = 4; smallmaf$mut_type[i] = 'D'; smallmaf$types[i] = 4
  } else if (type %in% cat5) { type = 5; smallmaf$mut_type[i] = 'E'; smallmaf$types[i] = 5
  } else if (type %in% cat6) { type = 6; smallmaf$mut_type[i] = 'F'; smallmaf$types[i] = 6
  } else if (type %in% cat7) { type = 7; smallmaf$mut_type[i] = 'G'; smallmaf$types[i] = 7
  } else if (type %in% cat8) { type = 8; smallmaf$mut_type[i] = 'H'; smallmaf$types[i] = 8
  } else if (type %in% cat9) { type = 9; smallmaf$mut_type[i] = 'I'; smallmaf$types[i] = 9
  } else {print(paste(i,type,sep="_"))
    stop("Mutation type not found")}
    print(paste(i,type,sep="_"))
  
  if (!TCGA) { 
    if (mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap)==smallmaf$newGene[i])] > type) {
      mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap)==smallmaf$newGene[i])] <- type
    }else { mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$TUMOR_SAMPLE[i]), which(colnames(mutation_heatmap)==smallmaf$Hugo[i])] <- type }
  }
}

###############################
# Order CCF-Group
CCFgroup_heatmap <-matrix(0, nrow=sum(unlist(lapply(sample_names, length))), ncol=sum(unlist(lapply(mutation_genes, length))))
rownames(CCFgroup_heatmap) <- unlist(sample_names)
colnames(CCFgroup_heatmap) <- mutation_genes


for (i in 1:nrow(smallmaf)) {
  type = smallmaf$CCF_group[i] 
  if (CCFgroup_heatmap[which(rownames(CCFgroup_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(CCFgroup_heatmap)==smallmaf$newGene[i])] < type) {
    CCFgroup_heatmap[which(rownames(CCFgroup_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(CCFgroup_heatmap)==smallmaf$newGene[i])] <- type}
}
i=1
for (i in 1 :nrow(mutation_heatmap)){
  mutation_heatmap = mutation_heatmap[, order(CCFgroup_heatmap[nrow(CCFgroup_heatmap)-(i-1),], decreasing = T) ]
  CCFgroup_heatmap = CCFgroup_heatmap[, order(CCFgroup_heatmap[nrow(CCFgroup_heatmap)-(i-1),], decreasing = T) ]
}

#######################################


if (sort_samples) {
  mutation_heatmap <- do.call(rbind, lapply(sample_names, function(x) { m <- mutation_heatmap[which(rownames(mutation_heatmap) %in% x),, drop=F]; m[do.call(order, transform(m)),]}))
}
rownames(mutation_heatmap) <- unlist(sample_names)

#### Sort genes by the number of mutations that do not equal the blank category
if(sort_genes) {
  print("Sorting genes")
  oo <- unlist(apply(mutation_heatmap,2, function(x){length(which(x!=9))}))
  print(oo)
  oo <- oo[which(!duplicated(names(oo)))]
  
  oo <- names(sort(oo, decr=T))
  mutation_heatmap <- mutation_heatmap[,match(unlist(oo), colnames(mutation_heatmap))]	
  
}


if(remove_genes_with_no_mutation) { mutation_heatmap <- mutation_heatmap[,which(unlist(apply(mutation_heatmap,2, function(x){length(which(x!=10))}))!=0)] }
#order_heatmap = ifelse(mutation_heatmap==9, 2, 1)
i=1
#for (i in 1 :ncol(mutation_heatmap)){
#  mutation_heatmap = mutation_heatmap[ order(order_heatmap[,ncol(mutation_heatmap)-(i-1)]) , ]
#  order_heatmap = order_heatmap[ order(order_heatmap[,ncol(mutation_heatmap)-(i-1)]) , ]
#}

## 
mutation_heatmap = t(mutation_heatmap)

#### Choose color palette for table - set to colkey
library(RColorBrewer)
colkey <- cbind(1:9, c(brewer.pal(8, "Set1"), "white"))
colkey[1,2]="#CD1719" 
colkey[2,2]="#984EA3"
colkey[3,2]="#377EB8"
colkey[4,2]="#4DAF4A"
colkey[5,2]="#FF7F00"
colkey[6,2]="#FFFF33"
colkey[7,2]="#A65628"
colkey[8,2]="#808080"
colkey[9,2]="gray85"


width=height=NULL
if (is.null(width)) { width = 1+(length(unlist(mutation_genes))/1.5) }
if (width<4){width=4}
if (is.null(height)) { height = 1+(length(unlist(sample_names)))/2 }
if (height<4){height=4}


mutation_heatmap <- t(mutation_heatmap)


#### Create empty pdf
height = height+1
pdf(plot_file, width=height, height=width)
if (show_sample_names) { top=8 } else {top = 2}
if (include_percentages) { right= 4 } else { right=1 }



#### Plot figure
par(oma=c(2,8,1,1), mar=c(2,5,top,right))
image(mutation_heatmap, xaxt='n', yaxt='n', col=colkey[,2], zlim=c(1,9), xlab="", ylab="")
axis(2, at=seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels=colnames(mutation_heatmap), las=2, tick=F, cex.axis=2, font = 3, family = 'sans')
axis(3, at=seq(0, 1, 1/(nrow(mutation_heatmap)-1)), labels=rownames(mutation_heatmap), las=2, tick=F, cex.axis=2, font = 1, family = 'sans')
abline(v=(0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)+(1/(2*(nrow(mutation_heatmap)-1)))), col = "white")
abline(h=(0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)+(1/(2*(ncol(mutation_heatmap)-1)))), col = "white")
dev.off()
mutation_heatmap <- t(mutation_heatmap)

##########################
## Mutation type ggplot ##
##########################
library(ggplot2)

smallmaf$CGC = ''
for (i in 1:nrow(smallmaf)) {
  if (smallmaf$cancer_gene_census[i] == "TRUE" || smallmaf$kandoth[i] == "TRUE" || smallmaf$lawrence[i] == "TRUE") { smallmaf$CGC[i] = '*'
  }
}
smallmaf$newGene <- paste(smallmaf$CGC, smallmaf$SYMBOL, sep="")

df11 <-as.data.frame(cbind(smallmaf$TUMOR_SAMPLE, smallmaf$newGene, smallmaf$Cancer_Cell_Fraction, smallmaf$Clonal_Status, smallmaf$loh, smallmaf$mut_type, smallmaf$types))

colnames(df11) <- c("samples", "genes", "ccf", "clone", "loh", "mut_type", "types")

df11 <- df11[order(df11$samples, df11$types, decreasing=TRUE),]

df11 <- as.data.frame(lapply(df11, function(x) {sub("DL-moc-0", "MOC", x)})) ## Modify samples names (if required)
df11 <- as.data.frame(lapply(df11, function(x) {sub("DL-moc-00", "MOC", x)})) ## Modify samples names (if required)

df11$genes <- factor(df11$genes, levels = dput(rev(rownames(mutation_heatmap))))

cols <- c('A'="#CD1719",'B'="#984EA3",'C'="#377EB8",'D'="#4DAF4A",'E'="#FF7F00",'F'="#FFFF33",'G'="#A65628", 'H' = "#808080", 'I' = "gray85")

comb<-expand.grid(samples = unique(df11$samples), genes = unique(df11$genes), ccf=".", clone = ".", loh=".", mut_type="I", types = 9)

df33 <- as.data.frame(rbind(comb, df11))

df22<-subset(df11,clone=="Clonal")

pdf("Mut_type_ggplot.pdf", width=12, height=89)

#p <- ggplot(df33, aes(x=samples,y=genes,fill=mut_type,group=loh)) + geom_tile(width=0.8,height=0.8) + geom_tile(data=df22,aes(samples, genes),size=1,fill=NA,width=0.8,height=0.8,color="gold3") + scale_fill_manual(values = cols, name="CCF",guide = FALSE); p + geom_segment( aes(x=xmin,xend=xmax,y=ymin,yend=ymax), subset(ggplot_build(p)$data[[1]],group==2), inherit.aes=F,color="white",size=0.5) + theme(axis.text.x=element_text(size=20,  family = "sans", angle = 90, hjust=0, vjust=0)) + theme(axis.text.y=element_text(size=20, face="italic", family = "sans")) + scale_x_discrete(name="",position="top",drop=FALSE) + scale_y_discrete(name="",drop=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank()) + theme(legend.position="none") + theme(axis.ticks = element_blank()) + theme(axis.text=element_text(colour = "black"))

p <- ggplot(df33, aes(x=samples,y=genes,fill=mut_type,group=loh)) + geom_tile(width=0.8,height=0.8) + scale_fill_manual(values = cols, name="CCF",guide = FALSE); p + geom_segment( aes(x=xmin,xend=xmax,y=ymin,yend=ymax), subset(ggplot_build(p)$data[[1]],group==2), inherit.aes=F,color="white",size=1) + theme(axis.text.x=element_text(size=20,  family = "sans", angle = 90, hjust=0, vjust=0)) + theme(axis.text.y=element_text(size=20, face="italic", family = "sans")) + scale_x_discrete(name="",position="top",drop=FALSE) + scale_y_discrete(name="",drop=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank()) + theme(legend.position="none") + theme(axis.ticks = element_blank()) + theme(axis.text=element_text(colour = "black"))

dev.off()
########################################
######## More than two ####
genes <- c("KEAP1", "TGFBR1", "TET1", "TERT", "TBX3", "SPOP", "SMO", "SF3B1", "PTPRT", "PTEN", "PRKAR1A", "PIK3C2G", "PDPK1", "PAK1", "NFE2L2", "MYOD1", "MRE11A", "MITF", "MEN1", "MED12", "MDM2", "LATS1", "KDM5A", "JAK2", "IL7R", "IDH2", "HGF", "GRIN2A", "GATA1", "FLT1", "FGFR1", "EZH2", "ETV1", "ESR1", "ERBB4", "DCUN1D1", "CTCF", "CRLF2", "CBL", "BTK", "BRIP1", "BMPR1A", "BLM", "BCL2", "YES1", "XPO1", "TP63", "TOP1", "TMPRSS2", "SUFU", "STAG2", "SOX2", "ROS1", "RET", "RASA1", "PTPRS", "PBRM1", "NTRK3", "NTRK2", "NSD1", "NOTCH3", "MYCN", "MET", "MAP3K13", "KIT", "KDM5C", "JAK1", "IGF1", "FAT1", "EPHB1", "EPHA5", "CARD11", "BCOR", "ATRX", "ASXL2", "ARID2", "APC", "AMER1", "KDM6A", "NTRK1", "PMS1", "HNF1A", "SHQ1", "PIK3R2", "EED", "MSH6", "U2AF1", "HRAS", "HIST1H1C", "CIC", "PHOX2B", "CCND1", "CDKN2A", "BAP1", "NKX2-1", "RB1", "KMT2A", "PTPRD", "PDGFRA", "E2F3", "PTCH1", "ALOX12B")

mutation_heatmap_red <- mutation_heatmap[!grepl(paste(genes, collapse="|"), rownames(mutation_heatmap)), ]
 
df11_red <- df11[!grepl(paste(genes, collapse="|"), df11$genes), ]

df11_red$genes <- factor(df11_red$genes, levels = dput(rev(rownames(mutation_heatmap_red))))

cols <- c('A'="#CD1719",'B'="#984EA3",'C'="#377EB8",'D'="#4DAF4A",'E'="#FF7F00",'F'="#FFFF33",'G'="#A65628", 'H' = "#808080", 'I' = "gray85")

comb<-expand.grid(samples = unique(df11_red$samples), genes = unique(df11_red$genes), ccf=".", clone = ".", loh=".", mut_type="I", types = 9)

df33_red <- as.data.frame(rbind(comb, df11_red))

df22_red <-subset(df11_red,clone=="Clonal")

pdf("Mut_type_ggplot_2.pdf", width=4, height=8)

p <- ggplot(df33_red, aes(x=samples,y=genes,fill=mut_type,group=loh)) + geom_tile(width=0.8,height=0.8) + scale_fill_manual(values = cols, name="CCF",guide = FALSE); p + geom_segment( aes(x=xmin,xend=xmax,y=ymin,yend=ymax), subset(ggplot_build(p)$data[[1]],group==2), inherit.aes=F,color="white",size=0.5) + theme(axis.text.x=element_text(size=6,  family = "sans", angle = 90, hjust=0, vjust=0)) + theme(axis.text.y=element_text(size=6, face="italic", family = "sans")) + scale_x_discrete(name="",position="top",drop=FALSE) + scale_y_discrete(name="",drop=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank()) + theme(legend.position="none") + theme(axis.ticks = element_blank()) + theme(axis.text=element_text(colour = "black"))

dev.off()
###########################


gene_order = rownames(mutation_heatmap)#[length(rownames(mutation_heatmap)):1]
samp_order = colnames(mutation_heatmap)

if(return_gene_order) { list(rev(colnames(mutation_heatmap))) }

return_gene_order   = F
sort_samples        = T
sort_genes          = T
show_sample_names   = T
TCGA                = F
remove_genes_with_no_mutation = F
# width               = NULL 
# height              = NULL
include_percentages = F
sample_name_col     = "TUMOR_SAMPLE"
mutation_genes=unique(muts$newGene)

#### Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples
mutation_heatmap <- matrix(1, nrow=sum(unlist(lapply(samp_order, length))), ncol=sum(unlist(lapply(gene_order, length))))
rownames(mutation_heatmap) <- samp_order
colnames(mutation_heatmap) <- gene_order

#### Make sure the sample and mutations are both in the list of gene mutations and gene samples
#if (!TCGA) { smallmaf <- muts[which(muts$ANN....GENE %in% unlist(mutation_genes) & muts$TUMOR_SAMPLE %in% unlist(sample_names)),]
#}else { 
#  muts$id <- unlist(lapply(muts$TUMOR_SAMPLE, function(x){substr(x, 1, 12)}))
#  print(head(muts$id))
# print(head(unlist(sample_names)))
#  smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),] 
#}

#### For each row read the Effect and create the type based on which category it fits in
for (i in 1:nrow(smallmaf)) {
  if(!TCGA) { type = smallmaf$Cancer_Cell_Fraction[i] } else { type = smallmaf$Variant_Classification[i] }
  if (type == 0) { type = 1
  } else if (type >0    & type<=0.05) { type = 2
  } else if (type >0.05 & type<=0.2)  { type = 3
  } else if (type >0.2  & type<=0.4)  { type = 4
  } else if (type >0.4  & type<=0.6)   { type = 5 
  } else if (type >0.6  & type<=0.8) { type = 6
  } else if (type >0.8  & type<=1) { type = 7
  } else { print("Mutation type not found") }
  print(paste(i,type,sep="_"))
  
  if (mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap)==smallmaf$newGene[i])] < type) {
    mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap)==smallmaf$newGene[i])] <- type}
}

if(remove_genes_with_no_mutation) { mutation_heatmap <- mutation_heatmap[,which(unlist(apply(mutation_heatmap,2, function(x){length(which(x!=10))}))!=0)] }

df=data.frame(matrix(ncol = 0, nrow = nrow(mutation_heatmap)))
for (i in 1:length(gene_order)){
  print(i)
  df = cbind.data.frame(df, mutation_heatmap[ , colnames(mutation_heatmap)== gene_order[i]])
  colnames(df)[i]=gene_order[i]                         
}
rownames(df)=rownames(mutation_heatmap)
newdf=data.frame(matrix(ncol = 0, nrow = ncol(mutation_heatmap)))
for (i in 1:length(samp_order)){
  print(i)
  newdf= rbind.data.frame(newdf, df[rownames(mutation_heatmap) == samp_order[i],])
  rownames(newdf)[i]=samp_order[i]                         
}
colnames(newdf)=colnames(df)
mutation_heatmap <- as.matrix(newdf)


#### Choose color palette for table - set to colkey
library(RColorBrewer)
colkey <- cbind(1:7, c("grey90",  "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5", "#08519C", "#08306B"))

mutation_heatmap = t(mutation_heatmap)


###################

#### Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples
LOH_heatmap <- matrix(0, ncol=sum(unlist(lapply(gene_order, length))), nrow=sum(unlist(lapply(samp_order, length))))
colnames(LOH_heatmap) <- gene_order
rownames(LOH_heatmap) <- samp_order

for (i in 1:nrow(smallmaf)) {
  type = smallmaf$loh[i]
  if (type == "loh") { type = 1
  # } else if (type >0    & type<=0.05) { type = 2
  } else { type = 0 }
  print(paste(i,type,sep="_"))
  
  if (LOH_heatmap[which(rownames(LOH_heatmap)==smallmaf$TUMOR_SAMPLE[i]), which(colnames(LOH_heatmap)==smallmaf$newGene[i],)] < type) {
    LOH_heatmap[which(rownames(LOH_heatmap)==smallmaf$TUMOR_SAMPLE[i]), which(colnames(LOH_heatmap)==smallmaf$newGene[i],)] <- type}
}

##########################################################################################################################
##########################################################################################################################

#### Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples
Clonal_heatmap <- matrix(0, ncol=sum(unlist(lapply(gene_order, length))), nrow=sum(unlist(lapply(samp_order, length))))
colnames(Clonal_heatmap) <- gene_order 
rownames(Clonal_heatmap) <- samp_order

for (i in 1:nrow(smallmaf)) {
  # type = smallmaf$clonality[i]
  # if (type == "clonal") { type = 1
  type = smallmaf$Clonal_Status[i]
  if (type == "Clonal") { type = 1
  # } else if (type >0    & type<=0.05) { type = 2
  } else { type = 0 }
  print(paste(i,type,sep="_"))
  
  if (Clonal_heatmap[which(rownames(Clonal_heatmap)==smallmaf$TUMOR_SAMPLE[i]), which(colnames(Clonal_heatmap)==smallmaf$newGene[i],)] < type) {
    Clonal_heatmap[which(rownames(Clonal_heatmap)==smallmaf$TUMOR_SAMPLE[i]), which(colnames(Clonal_heatmap)==smallmaf$newGene[i],)] <- type}
}

##########################################################################################################################
##########################################################################################################################
#### Create empty pdf
pdf(plot_file_CCF, height = height, width=width)
if (show_sample_names) { top=8 } else {top = 2}
if (include_percentages) { right= 4 } else { right=1 }



#Clonal_heatmap <- Clonal_heatmap[,rowSums(mutation_heatmap != 1) > 1]
#LOH_heatmap <- LOH_heatmap[,rowSums(mutation_heatmap != 1) > 1]
#mutation_heatmap <- mutation_heatmap[rowSums(mutation_heatmap != 1) > 1,]
##########################################################################################################################
##########################################################################################################################
#### Create empty pdf
#width=height=NULL
#if (is.null(width)) { width = 1+nrow(mutation_heatmap)/1.5 }
#if (width<4){width=4}
#if (is.null(height)) { height = 1+(length(unlist(sample_names)))/2 }
#if (height<4){height=4}

#pdf(plot_file_CCF, height = height, width=width)
#if (show_sample_names) { top=8 } else {top = 2}
#if (include_percentages) { right= 4 } else { right=1 }



#### Plot figure
par(oma=c(2,8,1,1), mar=c(2,5,top,right))
#par(mar=c(2,4,2,2))
image(mutation_heatmap, xaxt='n', yaxt='n', col=colkey[,2], zlim=c(1,7), xlab="", ylab="")
#axis(1, at=seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels=colnames(mutation_heatmap), las=2, tick=F)
axis(2, at=seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels=colnames(mutation_heatmap), las=2, tick=F, cex.axis=2, font = 3, family = 'sans')
axis(3, at=seq(0, 1, 1/(nrow(mutation_heatmap)-1)), labels=rownames(mutation_heatmap), las=2, tick=F, cex.axis=1.2, font = 2, family = 'sans')
if(include_percentages) {
  axis(4, at=seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels=unlist(apply(mutation_heatmap, 2, function(x) { paste(round(100*length(which(x!=1))/length(x), 0), "%")})), las=2, tick=F)}
# abline(v=(0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)+(1/(2*(nrow(mutation_heatmap)-1)))), col = "white")
# abline(h=(0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)+(1/(2*(ncol(mutation_heatmap)-1)))), col = "white")

row_pos = 0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)
row_pos = row_pos[!row_pos>1]
col_pos = 0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)
col_pos = col_pos[!col_pos>1]
# LOH_heatmap = LOH_heatmap[,ncol(LOH_heatmap):1]
# LOH_heatmap= t(LOH_heatmap)
LOH_heatmap = LOH_heatmap[nrow(LOH_heatmap):1,]
for (i in 1:nrow(LOH_heatmap)){
  row_plot = row_pos[as.vector(LOH_heatmap[i,]==1)]
  points( row_plot, rep(col_pos[(length(samp_order)+1)-i],length(row_plot)), pch="/", col="white", cex = 3)
}
Clonal_heatmap = Clonal_heatmap[nrow(Clonal_heatmap):1,]
for (i in 1:nrow(Clonal_heatmap)){
  row_plot = row_pos[as.vector(Clonal_heatmap[i,]==1)]
  points( row_plot, rep(col_pos[(length(samp_order)+1)-i],length(row_plot)), pch=16, col="yellow", cex = 2)
}

dev.off()

##########################
## CCF heatmap (ggplot) ##
##########################
library(ggplot2)

#df1<-as.data.frame(cbind(muts$TUMOR_SAMPLE, muts$ANN....GENE, muts$Cancer_Cell_Fraction, muts$Clonal_Status, muts$loh, muts$CCF_group))

df1 <-as.data.frame(cbind(smallmaf$TUMOR_SAMPLE, smallmaf$newGene, smallmaf$Cancer_Cell_Fraction, smallmaf$Clonal_Status, smallmaf$loh, smallmaf$CCF_group, smallmaf$types))

colnames(df1) <- c("samples", "genes", "ccf", "clone", "loh", "ccf_groupA", "types")

df1 <- df1[order(df1$samples, df1$types, decreasing=TRUE),]

df1$ccf_groupB = "0"
for (i in 1:nrow(df1)){
  if(df1$ccf_groupA[i] == 5){df1$ccf_groupB[i]= 'A'}
  if(df1$ccf_groupA[i] == 4){df1$ccf_groupB[i]= 'B'}
  if(df1$ccf_groupA[i] == 3){df1$ccf_groupB[i]= 'C'}
  if(df1$ccf_groupA[i] == 2){df1$ccf_groupB[i]= 'D'}
  if(df1$ccf_groupA[i] == 1){df1$ccf_groupB[i]= 'E'}
}

df1 <- as.data.frame(lapply(df1, function(x) {sub("DL-moc-0", "MOC", x)})) ## Modify samples names (if required)
df1 <- as.data.frame(lapply(df1, function(x) {sub("DL-moc-00", "MOC", x)})) ## Modify samples names (if required)

cols <- c('A'="#08306B",'B'="#08519C",'C'="#2171B5",'D'="#6BAED6",'E'="#9ECAE1",'F'="#C6DBEF",'G'="gray85")

df2<-subset(df1,clone=="Clonal")

df1$genes <- factor(df1$genes, levels = dput(rev(rownames(mutation_heatmap))))

comb<-expand.grid(samples = unique(df1$samples), genes = unique(df1$genes), ccf=".", clone = ".", loh=".", ccf_groupA = ".", types = ".", ccf_groupB="G")

df3 <- as.data.frame(rbind(comb, df1))

height=height+2
print (width)
print (height)
pdf("CCF_Heatmap_ggplot.pdf", width=height, height=width) ## Change size of the plot

#label_color <- c("black","red","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black")

p <- ggplot(df3, aes(x=samples,y=genes,fill=ccf_groupB, group=loh)) + geom_tile(width=0.8,height=0.8) + geom_tile(data=df2,aes(samples, genes),size=1,fill=NA,width=0.8,height=0.8,color="gold3") + scale_fill_manual(values = cols, name="CCF",guide = FALSE); p + geom_segment( aes(x=xmin,xend=xmax,y=ymin,yend=ymax), subset(ggplot_build(p)$data[[1]],group==2), inherit.aes=F,color="white",size=1) + theme(axis.text.x=element_text(size=20,  family = "sans", angle = 90, hjust=0, vjust=0)) + theme(axis.text.y=element_text(size=20, face="italic", family = "sans")) + scale_x_discrete(name="",position="top",drop=FALSE) + scale_y_discrete(name="",drop=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank()) + theme(legend.position="none") + theme(axis.ticks = element_blank()) + theme(axis.text=element_text(colour = "black"))

dev.off()

######## More than two ####
genes <- c("KEAP1", "TGFBR1", "TET1", "TERT", "TBX3", "SPOP", "SMO", "SF3B1", "PTPRT", "PTEN", "PRKAR1A", "PIK3C2G", "PDPK1", "PAK1", "NFE2L2", "MYOD1", "MRE11A", "MITF", "MEN1", "MED12", "MDM2", "LATS1", "KDM5A", "JAK2", "IL7R", "IDH2", "HGF", "GRIN2A", "GATA1", "FLT1", "FGFR1", "EZH2", "ETV1", "ESR1", "ERBB4", "DCUN1D1", "CTCF", "CRLF2", "CBL", "BTK", "BRIP1", "BMPR1A", "BLM", "BCL2", "YES1", "XPO1", "TP63", "TOP1", "TMPRSS2", "SUFU", "STAG2", "SOX2", "ROS1", "RET", "RASA1", "PTPRS", "PBRM1", "NTRK3", "NTRK2", "NSD1", "NOTCH3", "MYCN", "MET", "MAP3K13", "KIT", "KDM5C", "JAK1", "IGF1", "FAT1", "EPHB1", "EPHA5", "CARD11", "BCOR", "ATRX", "ASXL2", "ARID2", "APC", "AMER1", "KDM6A", "NTRK1", "PMS1", "HNF1A", "SHQ1", "PIK3R2", "EED", "MSH6", "U2AF1", "HRAS", "HIST1H1C", "CIC", "PHOX2B", "CCND1","CDKN2A", "BAP1", "NKX2-1", "RB1", "KMT2A", "PTPRD", "PDGFRA", "E2F3", "PTCH1", "ALOX12B")

mutation_heatmap_red <- mutation_heatmap[!grepl(paste(genes, collapse="|"), rownames(mutation_heatmap)), ]

df1_red <- df1[!grepl(paste(genes, collapse="|"), df1$genes), ]

df1_red$genes <- factor(df1_red$genes, levels = dput(rev(rownames(mutation_heatmap_red))))

cols <- c('A'="#08306B",'B'="#08519C",'C'="#2171B5",'D'="#6BAED6",'E'="#9ECAE1",'F'="#C6DBEF",'G'="gray85")

comb<-expand.grid(samples = unique(df1_red$samples), genes = unique(df1_red$genes), ccf=".", clone = ".", loh=".", ccf_groupA = ".", types = ".", ccf_groupB="G")

df3_red <- as.data.frame(rbind(comb, df1_red))

df2_red <-subset(df1_red,clone=="Clonal")

pdf("CCF_Heatmap_ggplot_2.pdf", width=5, height=10)

p <- ggplot(df3_red, aes(x=samples,y=genes,fill=ccf_groupB, group=loh)) + geom_tile(width=0.8,height=0.8) + geom_tile(data=df2_red,aes(samples, genes),size=0.5,fill=NA,width=0.8,height=0.8,color="gold3") + scale_fill_manual(values = cols, name="CCF",guide = FALSE); p + geom_segment( aes(x=xmin,xend=xmax,y=ymin,yend=ymax), subset(ggplot_build(p)$data[[1]],group==2), inherit.aes=F,color="white",size=0.5) + theme(axis.text.x=element_text(size=10,  family = "sans", angle = 90, hjust=0, vjust=0)) + theme(axis.text.y=element_text(size=10, face="italic", family = "sans")) + scale_x_discrete(name="",position="top",drop=FALSE) + scale_y_discrete(name="",drop=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank()) + theme(legend.position="none") + theme(axis.ticks = element_blank()) + theme(axis.text=element_text(colour = "black"))

dev.off()
###########################
