install.packages("ggplot2")
install.packages("gplots")
install.packages("RColorBrewer")

library(ggplot2)
library(gplots)
library(RColorBrewer)


setwd("C:/Users/ayoungk/Desktop/Heatmap")

AY <- read.table("File name.txt", sep="\t", header=T, row.names=1) #File name in " "
colnames(AY)
pdf("2.pdf", width=8, height=8)
data <- as.matrix(AY)
hmcol <- colorRampPalette (c("navyblue", "white", "orange"))
heatmap.2(data, col=hmcol, scale="row", density.info="none", trace="none", Colv="NA", dendrogram="row", key=T, cexRow=0.3, cexCol=1)
heatmap.2(data, col=hmcol, scale="row", density.info="none", trace="none", dendrogram="none", key=T, cexRow=0.5, cexCol=1)
dev.off()
