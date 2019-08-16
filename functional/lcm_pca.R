
INPUT <- "C:/Users/David/Desktop/analysis/allrpkms.txt"
PLOT <- "C:/Users/David/Desktop/allrpkms_PCA1.png"

library(FactoMineR)
library(ggplot2)

title <- "PCA of LCM & BMDM"
data2 <- read.delim(INPUT,header=TRUE)
data2 <- t(data2[,-1])
#data2 <- t(data/colSums(data)*1e6)
pca = PCA(data2,graph=FALSE)

pc1 <- pca$ind$coord[,1]
pc2 <- pca$ind$coord[,2]
cl <- c(rep("lcm",6),rep("bmdm1",3),rep("bmdm2",3),rep("lps10",3),rep("t0lps10",3),
        rep("lps100",3),rep("cl.lps100",3),rep("liver",10),rep("adipose",5),rep("spleen",6))
clf <- factor(cl,levels=c("lcm","bmdm1","bmdm2","lps10","t0lps10","lps100","cl.lps100","liver","adipose","spleen"))
pd <- data.frame(pc1,pc2,cl)
col_pick <- c("steelblue","yellowgreen","violetred1","plum1","cyan2","brown4","black","chartreuse","deeppink3","gray")

png(PLOT, width=6,height=4,units="in",res=600)
p <- ggplot(data=pd, aes(x=pc1,y=pc2,fill=clf)) + geom_point(shape=21,size=6) +
  theme_bw(base_size=20) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5,size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size_area(guide=FALSE) + scale_fill_manual(values=col_pick) + 
  guides(fill=guide_legend(title=NULL),colour=guide_legend(override.aes=list(size=12))) #+
print(p)
dev.off()