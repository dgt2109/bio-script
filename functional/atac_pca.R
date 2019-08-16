
INPUT <- "C:/Users/David/Desktop/Sequencing/AD002/atac_qc_metrics/atac-counts_macpeak-byrep-input.txt"
PLOT <- "C:/Users/David/Desktop/Sequencing/AD002/atac_qc_metrics/atac-counts_macpeak_PCA.eps"

library(FactoMineR)
library(ggplot2)

title <- "PCA of ATAC replicates"
data <- read.delim(INPUT,header=TRUE)
data2 <- t(data/colSums(data)*1e6)
pca = PCA(data2,graph=FALSE)

pc1 <- pca$ind$coord[,1]
pc2 <- pca$ind$coord[,2]
cl <- c(rep("Veh",4),rep("T0",4),rep("LPS",4),rep("T0+LPS",4))
clf <- factor(cl,levels=c("Veh","T0","LPS","T0+LPS"))
pd <- data.frame(pc1,pc2,cl)
col_pick <- c("steelblue","yellowgreen","violetred1","plum1")

setEPS()
postscript(PLOT, width=6,height=4)
p <- ggplot(data=pd, aes(x=pc1,y=pc2,fill=clf)) + geom_point(shape=21,size=6) +
  theme_bw(base_size=20) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5,size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size_area(guide=FALSE) + scale_fill_manual(values=col_pick) + 
  guides(fill=guide_legend(title=NULL),colour=guide_legend(override.aes=list(size=12))) #+
print(p)
dev.off()