install.packages("ggplot2")
install.packages("gplots")
install.packages("RColorBrewer")

library(ggplot2)
library(gplots)
library(RColorBrewer)


setwd("C:/Users/ayoungk/Desktop")

AY <- read.table("Duodenum.txt", sep="\t", header=T, row.names=1) # Name of the file in between " "
colnames(AY)
pdf("Duodenum_Cholesterol efflux_.pdf", width=10, height=10)
#One of data
data <- as.matrix(AY)
hmcol <- colorRampPalette (c("navyblue", "white", "orange"))
#One of below heatmap.2
heatmap.2(data, col=hmcol, scale="row", density.info="none", trace="none", Rowv = "NA", Colv="NA", dendrogram="none", key=T, cexRow=0.5, cexCol=1)

dev.off()



#Sort specific lane information
data1 <- read.table("Gene.txt", sep="\t", header=T)
results <- read.table("Sham_SBR.txt", sep="\n", header=T)
pos <- match(toupper(results[,1]), toupper(data1[,1]))
data2 <- data1[pos,]
write.table(data2, "Test.txt", sep="\t", row.names=F, col.names=T, quote=F)

install.packages("cluterProfiler")
