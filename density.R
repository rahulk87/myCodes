library(ggplot2)

brca1 <- read.table("BRCA1_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
d_sig3 <- density(brca1$Signature.3, na.rm=TRUE)
d_LST <- density(brca1$LST, na.rm=TRUE)
pdf("brca1.pdf")
plot(d_sig3, main="BRCA1 - Signature 3 density", xlim=c(0,1), ylim=c(0,20))
polygon(d_sig3, col="red", border="blue")
plot(d_LST, main="BRCA1 - LST density", xlim=c(0,60))
polygon(d_LST, col="red", border="blue")
abline(v=15, lty=2, col="green")
dev.off()

pdf("brca1_ggplot.pdf")
ggplot(brca1, aes(Signature.3, LST)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  #geom_point() +
  scale_x_continuous(limits = c(0,0.3)) + scale_y_continuous(limits = c(0,30)) +
  labs(color="Density",fill="Density")
dev.off()

brca2 <- read.table("BRCA2_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
d_sig3 <- density(brca2$Signature.3, na.rm=TRUE)
d_LST <- density(brca2$LST, na.rm=TRUE)
pdf("brca2.pdf")
plot(d_sig3, main="BRCA2 - Signature 3 density", xlim=c(0,1), ylim=c(0,20))
polygon(d_sig3, col="red", border="blue")
plot(d_LST, main="BRCA2 - LST density", xlim=c(0,60))
polygon(d_LST, col="red", border="blue")
abline(v=15, lty=2, col="green")
dev.off()

pdf("brca2_ggplot.pdf")
ggplot(brca2, aes(Signature.3, LST)) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') +
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  #geom_point() +
  scale_x_continuous(limits = c(0,0.3)) + scale_y_continuous(limits = c(0,30))
  labs(color="Density",fill="Density")
dev.off()

atm <- read.table("ATM_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
d_sig3 <- density(atm$Signature.3, na.rm=TRUE)
d_LST <- density(atm$LST, na.rm=TRUE)
pdf("atm.pdf")
plot(d_sig3, main="ATM - Signature 3 density", xlim=c(0,1), ylim=c(0,20))
polygon(d_sig3, col="red", border="blue")
plot(d_LST, main="ATM - LST density", xlim=c(0,60))
polygon(d_LST, col="red", border="blue")
abline(v=15, lty=2, col="green")
dev.off()

pdf("atm_ggplot.pdf")
ggplot(atm, aes(Signature.3, LST)) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') +
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") + 
  #geom_point() +
  scale_x_continuous(limits = c(0,0.3)) + scale_y_continuous(limits = c(0,30))
  labs(color="Density",fill="Density")
dev.off()
