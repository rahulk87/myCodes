#### A<-c("ATOH8","NOX1","SSH1","TMEM219","TMUB2"); ##
library(VennDiagram)
venn.plot <- venn.diagram(x = list(A = A,B = B,C = C),filename = NULL,col = "transparent",fill = c("cornflowerblue", "green", "orange"),alpha = 0.50,cex = 1.5,fontfamily = "serif",fontface = "bold",cat.col = c("cornflowerblue", "green", "orange"),cat.cex = 1.5,cat.pos =  c(-210, -130, -50),cat.dist = c(0.02, 0.02, 0.02),cat.fontfamily = "serif",rotation.degree =200);
jpeg("test.jpg");
grid.draw(venn.plot);
dev.off();
