rm(list=ls())
library(ggplot2)

load("samples_atm_BRCA12_germline_vus_loh1_noBiallelicHits.RData")

#brca1 <- read.table("BRCA1_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
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
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', size=0.1, bins=50) + 
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  geom_point() +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) +
  labs(color="Density",fill="Density")
dev.off()

#brca2 <- read.table("BRCA2_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
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
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', size=0.1, bins=50) +
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  geom_point() +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) +
  labs(color="Density",fill="Density")
dev.off()

#atm <- read.table("ATM_vus_loh1_noBiallelic_hits.tsv", sep="\t", header=T)
d_sig3 <- density(atm$Signature.3, na.rm=TRUE)
d_LST <- density(atm$LST, na.rm=TRUE)
pdf("atm.pdf")
plot(d_sig3, main="ATM - Signature 3 density", xlim=c(0,1), ylim=c(0,20))
polygon(d_sig3, col="red", border="blue")
plot(d_LST, main="ATM - LST density", xlim=c(0,60))
polygon(d_LST, col="red", border="blue")
abline(v=15, lty=2, col="green")
dev.off()

pdf("atm_ggplot.pdf", )
ggplot(atm, aes(Signature.3, LST)) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', size=0.1, bins=50) +
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") + 
  geom_point() +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) +
  labs(color="Density",fill="Density") + facet_grid(. ~ Hugo_Symbol)
dev.off()

## combined plots
comb <- rbind(brca1, brca2, atm)
pdf("combined_plots.pdf")
ggplot(comb, aes(x = Signature.3)) + geom_density(fill = "gold1", colour = "goldenrod2", alpha=0.6) + facet_grid(. ~ Hugo_Symbol) + scale_x_continuous(limits = c(0,1)) + theme_bw()
ggplot(comb, aes(x = LST)) + geom_density(fill = "gold1", colour = "goldenrod2", alpha=0.6) + facet_grid(. ~ Hugo_Symbol) + geom_vline(xintercept = 15, size = 0.5, colour = "#FF3721", linetype = "dashed") + theme_bw()
dev.off()

pdf("combined_density.pdf", width=15, height=5)
ggplot(comb, aes(Signature.3, LST)) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black', size=0.1, bins=50) +
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  geom_point(size=0.5) +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,50)) +
  labs(color="Density",fill="Density") + facet_grid(. ~ Hugo_Symbol)
dev.off()
