a=read.csv("data.txt",sep="\t",header=TRUE)
breast=subset(a,Tissue=="breast")

fileName <- "tissue"
conn <- file(fileName,open="r")
linn <-readLines(conn)

pdf("box.pdf")
for (i in 1:length(linn)){
	tiss=linn[i]
	print (tiss)
	tissue=subset(a,Tissue==tiss)
	l=length(a)
	write.table(tiss, file = paste(tiss,"csv",sep="."), row.names = FALSE,append = TRUE, col.names = FALSE,sep="\t")
		for (i in 3:l){
		profile=data.frame(cbind(breast[,i],tissue[,i]))
		p=t.test(breast[,i],tissue[,i])$p.value
		write.table(p, file = paste(tiss,"csv",sep="."), row.names = FALSE,append = TRUE, col.names = FALSE,sep="\t")
		boxplot(profile,main=tiss,names=c("Breast",paste(tiss)))
		stripchart(profile,vertical=TRUE,method = "jitter", pch = 19,add = TRUE,col=c("tomato4","seagreen3"))
		#stripchart(act,vertical = TRUE, method = "jitter", pch = 19,add = TRUE,col=c("tomato4","seagreen3"))
		}

	}
	close(conn)
dev.off()	
system("head -1 data.txt |cut -f2- |perl -pi -e 's/\t/\n/g' > header.csv", intern = TRUE, ignore.stderr = TRUE)
system("paste *.csv >final.txt", intern = TRUE, ignore.stderr = TRUE)
system("rm *.csv", intern = TRUE, ignore.stderr = TRUE)
