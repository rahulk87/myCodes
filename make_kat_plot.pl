#!/usr/bin/perl
for($i=1;$i<=7;$i++)
{
open(FA,">>kataegis.R");
print FA "png(\"kat$i.png\",height=1500,width=4000,res=300)\n";
print FA "x=read.csv(\"file$i\",sep=\"\\t\",header=TRUE)\n";
print FA "with(x,plot(Position,log10(Dist),main=\"ICRGFTJN$i\",xlab=\"Genomic Position\",ylab=\"Intermutation Distance (log10bp)\",pch=19,cex=0.4,cex.axis=.8))\n";
print FA "with(subset(x, Type == \"C>A\"), points(Position,log10(Dist), pch=20, col=\"blue\", border=\"blue\"))\n";
print FA "with(subset(x, Type == \"C>G\"), points(Position,log10(Dist), pch=20, col=\"yellow\", border=\"yellow\"))\n";
print FA "with(subset(x, Type == \"C>T\"), points(Position,log10(Dist), pch=20, col=\"red\", border=\"red\"))\n";
print FA "with(subset(x, Type == \"T>A\"), points(Position,log10(Dist), pch=20, col=\"green\", border=\"green\"))\n";
print FA "with(subset(x, Type == \"T>C\"), points(Position,log10(Dist), pch=20, col=\"black\", border=\"black\"))\n";
print FA "with(subset(x, Type == \"T>G\"), points(Position,log10(Dist), pch=20, col=\"purple\", border=\"purple\"))\n";
print FA "with(x,legend(\"topleft\",legend=c(\"C>A\",\"C>G\",\"C>T\",\"T>A\",\"T>C\",\"T>G\"),col=c(\"blue\",\"yellow\",\"red\",\"green\",\"black\",\"purple\"),pch=19,cex=.5))\n";
print FA "dev.off()\n";
close FA;
`R CMD BATCH kataegis.R`;
`rm kataegis.R kataegis.Rout`;
}
