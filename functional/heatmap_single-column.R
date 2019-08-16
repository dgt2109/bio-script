library(gplots)
library(RColorBrewer)
library(DESeq2) 

#INPUT2="C:/Users/David/Desktop/Sequencing/AD005/hld_counts.txt"
#OUTPUT2="C:/Users/David/Desktop/Sequencing/AD005/vf-hld_hld-effect.png"
#INPUT2="C:/Users/David/Desktop/Sequencing/AD001/T0-counts.txt"
#OUTPUT2="C:/Users/David/Desktop/Sequencing/AD001/LPS_column.eps"
INPUT2="C:/Users/David/Desktop/Sequencing/AD003/ifnar_counts.txt"
OUTPUT2="C:/Users/David/Desktop/Sequencing/AD003/ifnar_sc.svg"

#GROUPS=c("wt_wt","wt_ldlr")
#GROUPS=c("Ctrl","LPS","LPSHDL")
#GROUPS=c("Veh","LPS","T0.LPS")
#GROUPS=c("Ctrl","LPS","LPSHDL")
GROUPS=c("WT","LPS","IFNAR_LPS")
COMP2=c(2,3)
REPL=2
FDR=5

data<-read.delim(INPUT2,header=TRUE,stringsAsFactors=TRUE)
rownames(data) = make.names(data$Gene, unique=TRUE)
bp <- data[,6]
counts <- data[,-c(1:6)]
rpm <- sweep(counts,2,colSums(counts),'/')*1e6
rpkm <- sweep(rpm,1,bp,'/')*1e3

r3 <- ((COMP2[1]-1)*REPL+1):(COMP2[1]*REPL)
r4 <- ((COMP2[2]-1)*REPL+1):(COMP2[2]*REPL)
fc2 <- apply(rpkm[,r4],1,mean)/apply(rpkm[,r3],1,mean)

# pass rlist by running pairwise_rna_analysis.R with HEATMAP on
if(exists("rlist") && exists("sig_genes")) {
  fc2_sig <- data.frame(log2(data.matrix(fc2[names(fc2) 
                                             %in% sig_genes])+0.001))
  fc2_sig_rlist <- merge(rlist,fc2_sig,by.x=1,by.y=0)
  colnames(fc2_sig_rlist) <- c("name","h1_order","fc_comp2")
  fc2_sig_ordered <- fc2_sig_rlist[order(as.numeric(as.character(
    fc2_sig_rlist$h1_order))),]
  x <- fc2_sig_ordered$fc_comp2
#  x[265] = min(x[-c(265)])
  x[x == Inf] <- max(x[-c(which(x == Inf))])
  
  # ususally rev(colorRampPalette)
  hmcol <- colorRampPalette(brewer.pal(6,"PRGn"))(100)
  svg(OUTPUT2,width=3,height=8)
  h2 <- heatmap.2(cbind(x,x),Rowv=FALSE,Colv=FALSE,
                  scale="column", trace="none", dendrogram="none",
                  col=hmcol, labCol=FALSE,density.info="none",
                  labRow=FALSE,lhei=c(1,6),lwid=c(1.5,1),
#                  breaks=(0:100-50)/(100/4),
                  key.par=list(cex=0.6))
  dev.off()
}