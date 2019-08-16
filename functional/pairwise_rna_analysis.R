library(gplots)
library(RColorBrewer)
library(DESeq2) 

INPUT="C:/Users/David/Desktop/Sequencing/AK001/ku_counts_refname.txt"
BASEDIR="C:/Users/David/Desktop/Sequencing/AK001/"
OUTPUT="C:/Users/David/Desktop/Sequencing/AK001/HM_FDR10.svg"

#GROUPS=c("Lean_HDL","Lean_HDL_TNF","MHO_HDL","MHO_HDL_TNF",
#         "MUO_HDL","MUO_HDL_TNF","UTC","UTC_TNF")
#GROUPS=c("UTC","Lean_HDL","MHO_HDL","MUO_HDL",
#         "UTC_TNF","Lean_HDL_TNF","MHO_HDL_TNF","MUO_HDL_TNF")
#GROUPS=c("LPS10","LCM")
#GROUPS=c("Veh","LPS","T0.LPS")
#GROUPS=c("LPS","HDL.LPS","CL.LPS","CL.HDL.LPS")
#GROUPS=c("WT_CtrlASO","WT_T39ASO","KO_CtrlASO","KO_T39ASO")
GROUPS=c("Ctrl","MylKO")
COMPS=matrix(c(1,2),nrow=1)
#COMPS=matrix(c(7,7,7,7,8,8,8,8,1,3,5,2,4,6),7,2)
FILES=paste0("Ctrl_MylKO","_FDR10.csv")
#FILES=c("utc_muoHDL_minus1.csv")
HEATMAP=TRUE
TABLE=TRUE
REPL=c(5,5)
FDR=10
THRESHOLD=0

data<-read.delim(INPUT,header=TRUE,stringsAsFactors=TRUE)
rownames(data) = make.names(data$Gene, unique=TRUE)
bp <- data[,6]
counts <- data[,-c(1:6)]
rpm <- sweep(counts,2,colSums(counts),'/')*1e6
rpkm <- sweep(rpm,1,bp,'/')*1e3

# calculate RPKM and FC_RPKM
for (i in 1:nrow(COMPS)) {
  COMP = COMPS[i,]
  print(COMP)
  STATS = paste0(BASEDIR,FILES[i])
  print(STATS)
  
  r1 <- (sum(REPL[1:COMP[1]-1])+1):(sum(REPL[1:COMP[1]]))
  r2 <- (sum(REPL[1:COMP[2]-1])+1):(sum(REPL[1:COMP[2]]))
  fc_rpkm <- apply(rpkm[,r2],1,mean)/apply(rpkm[,r1],1,mean)
  
  # run DESeq2 on the relevant count matrix
  counts_comp <- cbind(counts[,r1],counts[,r2])
  columndata <- data.frame(row.names = colnames(counts_comp),
                            condition=c(rep(GROUPS[COMP[1]],REPL[COMP[1]]),
                                        rep(GROUPS[COMP[2]],REPL[COMP[2]])),
                            libtype=rep("single-end",REPL[COMP[1]]+REPL[COMP[2]]))
  dds <- DESeqDataSetFromMatrix(countData = counts_comp, 
                                colData = columndata, 
                                design =~ condition)
  dds$condition <- factor(dds$condition, levels = c(GROUPS[COMP[1]],
                                                      GROUPS[COMP[2]]))
  dds <- DESeq(dds,fitType="local",minReplicatesForReplace = Inf)
  res <- results(dds)
  
  # process DESeq2 output and append FC RPKM
  fullres <- merge(res,fc_rpkm,by=0)
  fullres <- fullres[!is.na(fullres$padj),]
  rownames(fullres) <- fullres[[1]]
  filterres <- fullres[fullres$padj<(FDR/100),][-1]
  colnames(filterres)[7] <- "fc_rpkm"
  fcfilterres <- filterres[abs(log(filterres$fc_rpkm))>log(THRESHOLD), ]
  sortres <- fcfilterres[order(-fcfilterres$fc_rpkm),]
  if (TABLE) { write.csv(as.data.frame(sortres),file=STATS) }
}

# plot heatmap
sig_genes <- row.names(sortres)
down_genes <- row.names(sortres[sortres$fc_rpkm<1,])
up_genes <- row.names(sortres[sortres$fc_rpkm>1,])

if (HEATMAP) {
  plotdata <- log2(data.matrix(rpkm[row.names(rpkm)
                                    %in% sig_genes, ])+0.001)
# plotdata <- cbind(plotdata[,r1],plotdata[,r2])
  hmcol <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(100))
  
  svg(OUTPUT,width=8,height=8)
  h <- heatmap.2(plotdata,Rowv=TRUE,Colv=FALSE,
                  distfun = function(x) as.dist(1-cor(t(plotdata))),
                  scale="row", trace="none", dendrogram="none",
                  col=hmcol, labCol=FALSE, density.info="none",
                  labRow=FALSE, lhei=c(1,6), lwid=c(1,4),
                  breaks=(0:100-50)/(100/4),
                  key.par=list(cex=0.6))
  dev.off()
  rlist <- data.frame(rownames(plotdata[rev(h$rowInd),h$colInd]))
  rlist <- cbind(rlist,rownames(rlist))
}