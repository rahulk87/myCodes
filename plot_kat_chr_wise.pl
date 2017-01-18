#!/usr/bin/perl
for($i=1;$i<=7;$i++)
{
open(FA,">>kataegis.R");
print FA "png(\"ICRGFTJN$i.png\",height=1500,width=4000,res=300)\n";
print FA "x=read.csv(\"file$i\",sep=\"\\t\",header=TRUE)\n";
print FA "with(x,plot(Position,log10(Dist),main=\"ICRGFTJN$i\",ylab=\"Intermutation Distance (log10bp)\",xlab=\"Chromosomes\",pch=19,cex=0.4,cex.axis=.8,xaxt=\"n\"))\n";
print FA "with(x,axis(1, at=c(124625310, 370850307, 591461209, 786049562, 972084330, 1148099493, 1313226358, 1465977701, 1609766427, 1748140516, 1883411148, 2017840353, 2142351240, 2253610949, 2358551415, 2454994487, 2540769469, 2620405698, 2689008813, 2750086065, 2805663772, 2855381003, 2958668566,3120000000),labels=c(\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\",\"11\",\"12\",\"13\",\"14\",\"15\",\"16\",\"17\",\"18\",\"19\",\"20\",\"21\",\"22\",\"X\",\"Y\")))\n";
print FA "with(x,abline(v=1,col = 4,lty=2))\n";
print FA "with(x,abline(v=249250621,col = 4,lty=2))\n";
print FA "with(x,abline(v=492449994,col = 4,lty=2))\n";
print FA "with(x,abline(v=690472424,col = 4,lty=2))\n";
print FA "with(x,abline(v=881626700,col = 4,lty=2))\n";
print FA "with(x,abline(v=1062541960,col = 4,lty=2))\n";
print FA "with(x,abline(v=1233657027,col = 4,lty=2))\n";
print FA "with(x,abline(v=1392795690,col = 4,lty=2))\n";
print FA "with(x,abline(v=1539159712,col = 4,lty=2))\n";
print FA "with(x,abline(v=1680373143,col = 4,lty=2))\n";
print FA "with(x,abline(v=1815907890,col = 4,lty=2))\n";
print FA "with(x,abline(v=1950914406,col = 4,lty=2))\n";
print FA "with(x,abline(v=2084766301,col = 4,lty=2))\n";
print FA "with(x,abline(v=2199936179,col = 4,lty=2))\n";
print FA "with(x,abline(v=2307285719,col = 4,lty=2))\n";
print FA "with(x,abline(v=2409817111,col = 4,lty=2))\n";
print FA "with(x,abline(v=2500171864,col = 4,lty=2))\n";
print FA "with(x,abline(v=2581367074,col = 4,lty=2))\n";
print FA "with(x,abline(v=2659444322,col = 4,lty=2))\n";
print FA "with(x,abline(v=2718573305,col = 4,lty=2))\n";
print FA "with(x,abline(v=2781598825,col = 4,lty=2))\n";
print FA "with(x,abline(v=2829728720,col = 4,lty=2))\n";
print FA "with(x,abline(v=2881033286,col = 4,lty=2))\n";
print FA "with(x,abline(v=3036303846,col = 4,lty=2))\n";
print FA "with(subset(x, Type == \"C>A\"), points(Position,log10(Dist), pch=20, col=\"blue\", border=\"blue\"))\n";
print FA "with(subset(x, Type == \"C>G\"), points(Position,log10(Dist), pch=20, col=\"yellow\", border=\"yellow\"))\n";
print FA "with(subset(x, Type == \"T>A\"), points(Position,log10(Dist), pch=20, col=\"green\", border=\"green\"))\n";
print FA "with(subset(x, Type == \"T>C\"), points(Position,log10(Dist), pch=20, col=\"black\", border=\"black\"))\n";
print FA "with(subset(x, Type == \"T>G\"), points(Position,log10(Dist), pch=20, col=\"purple\", border=\"purple\"))\n";
print FA "with(subset(x, Type == \"C>T\"), points(Position,log10(Dist), pch=20, col=\"red\", border=\"red\"))\n";
print FA "with(x,legend(\"topleft\",legend=c(\"C>A\",\"C>G\",\"C>T\",\"T>A\",\"T>C\",\"T>G\"),col=c(\"blue\",\"yellow\",\"red\",\"green\",\"black\",\"purple\"),pch=19,cex=.5))\n";
print FA "dev.off()\n";
close FA;
`R CMD BATCH kataegis.R`;
`rm kataegis.R`;
}
