library(ggplot2)
library(grid)

INPUT <- "C:/Users/David/Desktop/Sequencing/AD005/092818_go/Rfiles/GO_down-vf_nl-input.txt"
OUTPUT <- "C:/Users/David/Desktop/Sequencing/AD005/092818_go/Rfiles/GO_down-vf_nl.svg"

num_cats=16
#br=seq(1,12,1); lim=c(0.5,12); ms=10; yl=c(0,25); title="T0-repressed genes GO"
#br=seq(1,3,1); lim=c(0.5,8); ms=20; yl=c(0,25); title="T0-induced genes GO" 
#br=seq(1,10,1); lim=c(0.5,10); ms=18; yl=c(0,5); title="T0-closed nearest genes GO"
#br=seq(1,10,1); lim=c(0.5,10); ms=10; yl=c(0,10); title="LPS-opened nearest genes GO"
#br=seq(1,10,1); lim=c(0.5,10); ms=13; yl=c(0,10); title="T0-closed-in-LPS nearest genes GO"
#br=seq(1,12,1); lim=c(0.5,10); ms=15; yl=c(0,5); title="Unfilt. T0-closed nearest genes GO"

#br=seq(1,16,2); lim=c(0.5,16); ms=11; yl=c(0,12); title="VF-induced genes in normolipidemia"
br=seq(1,16,2); lim=c(0.5,16); ms=9; yl=c(0,16); title="VF-repressed genes in normolipidemia"
#br=seq(1,10,1); lim=c(0.5,10); ms=12; yl=c(0,20); title="VF-induced genes in hyperlipidemia"
#br=seq(1,12,1); lim=c(0.5,12); ms=12; yl=c(0,10); title="VF-repressed genes in hyperlipidemia"
#br=seq(1,8,1); lim=c(0.5,8); ms=12; yl=c(0,15); title="rHDL-induced genes"
#br=seq(1,8,1); lim=c(0.5,8); ms=12; yl=c(0,18); title="rHDL-repressed genes"
#br=seq(1,8,1); lim=c(0.5,8); ms=12; yl=c(0,15); title="rHDL-induced genes, splenic Cd11b+ cells"
#br=seq(1,10,1); lim=c(0.5,10); ms=12; yl=c(0,30); title="Hyperlipidemia-induced genes"
#br=seq(1,10,1); lim=c(0.5,10); ms=12; yl=c(0,220); title="T0-induced genes, splenic Cd11b+ cells"
#br=seq(1,6,1); lim=c(0.5,10); ms=12; yl=c(0,12); title="T0-repressed genes, splenic Cd11b+ cells"
#br=seq(1,12,1); lim=c(0.5,12); ms=12; yl=c(0,20); title="LIG induced in lesion vs BMDM"
#br=seq(1,12,1); lim=c(0.5,12); ms=12; yl=c(0,50); title="LIG repressed in lesion vs BMDM"
#br=seq(0,180,30); lim=c(0,180); ms=12; yl=c(0,28); title="T39ASO-induced in adipose tissue"
#br=seq(1,9,1); lim=c(0.5,9.5); ms=12; yl=c(0,180); title="T39ASO-induced in adipose tissue"

data <- read.delim(INPUT,header=TRUE)
data <- cbind(as.numeric(as.character(rownames(data))),data)
colnames(data)[1] <- "rank"
#colnames(data)[3] <- "Number\nof Genes"
data <- head(data,n=num_cats)

#setEPS()
#postscript(OUTPUT, width=8,height=5)
svg(OUTPUT, width=8.5,height=5)
p <- ggplot(data, aes(x = rank, y = FE, size = num_genes)) + 
  geom_point(shape = 21,colour="#000999",fill="#000999") + 
  scale_x_continuous(breaks = br, limits=lim) + 
  ylim(yl) + ggtitle(title) + labs(x="Rank",y="Fold Enrichment") + 
  scale_size_area(max_size=ms,name="Number\nof Genes") + 
  theme_bw(base_size=24) + 
  theme(plot.title=element_text(hjust=0.5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
print(p)
dev.off()
