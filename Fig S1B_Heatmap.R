install.packages("ggplot2")
install.packages("gplots")
install.packages("RColorBrewer")

library(ggplot2)
library(gplots)
library(RColorBrewer)


setwd("C:/Users/ayoungk/Desktop/Heatmap_LXR targets")

AY <- read.table("LXR targets.txt", sep="\t", header=T, row.names=1)
colnames(AY)
pdf("2.pdf", width=8, height=8)
data <- as.matrix(AY)
hmcol <- colorRampPalette (c("navyblue", "white", "orange"))
heatmap.2(data, col=hmcol, scale="row", density.info="none", trace="none", Colv="NA", dendrogram="row", key=T, cexRow=0.3, cexCol=1)

