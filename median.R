grpA <- groupA[,4:8]
rownames(grpA) <- groupA$class
grpA_median <- as.data.frame(apply(grpA, 1, median))
colnames(grpA_median) <- "grpA_median"
