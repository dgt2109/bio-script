library(gplots)
library(RColorBrewer)

#genes_LCCA <- c("Icam1","Itgb2","Cx3cr1","Il1b","Ccl5","Vcam1","Olr1")
#genes_CRC <- c("Cxcl3","Dock8","Cxcl1","Ccl3","Ccl5","Zc3h12a","Cxcl2")
#genes_GC <- c("Cxcl3","Itgb2","Cxcl1","Il1b","Ccl3","Anxa1","Ccl5","Cxcl2")

#hdlup_cholsynth <- c("Sc5d","Lss","Cyp51","Fdps","Insig1","Dhcr7","Fdft1","Pmvk","Hsd17b7","Dhcr24")
#hdlup_chemokine <- c("Cxcl3","Ccl17","Ccl7","Cxcl1","Cxcl5","Ccl3","Ccl22","Cxcl2","Ccl2")
#hdlup_nchemotaxis <- c("Cxcl3","Ccl17","Ccl7","Cxcl1","Cxcl5","Ccl3","Ccl22","Cxcl2","Ccl2")

#hdldown_ifnb <- c("Gbp3.1","Tgtp1","Tgtp2","Gm4951","Gm12185","Pyhin1","Gm4841","Ifit3","F830016B08Rik")
#hdldown_regoflymph <- c("Pde5a","Sdc4","Tarm1","Spta1","Fap","Cd86","Ccr2","Cd38","Il2ra","Slfn1","Cd74","Ido1","Aif1","Tlr9","Il6ra","Bst1")

#hdldown_erstress <- c("Wfs1","Creb3l2","Dnajc3","Dnajb9","Sdf2l1")
#lcm_up <- c("Hbegf","Cxcl16","Ch25h","Cx3cl1","Trem1","Cxcl5","Il17ra","Cxcl9","Itga1","Cxcl12.2","Vcam1","Spp1.4")
#lcm_down <- c("Cxcl3","Pf4","Ccrl2","Ccl9","Abcc1","Pip5k1a","Il1rn.1","Cxcl10","Il1b","Ccl3","Ccl5","Ccl4","Sbds","Ednrb.1","Rhog","Ccl22","Cxcl2","Ccl24","Pde4b.4")

#aT39_ifnb <- c("Gbp2","Gm4951","Ifi204","Irgm1","Tgtp2","Igtp","Ifit3","Ifit1","Tgtp1","Irgm2","F830016B08Rik","Mnda")
#aT39_drvirus <- c("Trim30a","Zbp1","Oasl2","Oas2","Ifitm3","Isg15","Oas1g","Mx1","Dhx58","Trim34a","Rsad2","Ifit3","Ddx60","Ifit1","Bst2","Ifit2","Rtp4","Oas3","Stat2","Slfn8.1","Oasl1","Oas1a")
#aT39_rtbact <- c("Gbp2","Gbp10","Hist2h2be","Trim30a","Herc6","Ifi204","Oas2","Irgm1","Isg15","Nlrc5","Dhx58","Gbp5","Usp18","Ifit3","Ifit1","Tgtp1","Irgm2","Cmpk2","Slfn4","Ifi44","Psmb9")

vf_kinetochore <- c("Smc2","Cenpe","Cenph","Cenpn","Ndc80","Cenpa","Nuf2","Cenpw")
vf_g2m <- c("Hspa2","Recql4","Cdk1","Cdc25a","Cdc25c","Ccnb1","Atad5","Dtl","Dyrk3","Fbxo5")
vf_il2 <- c("Cd28","Cd3e","Prkcq","Irf4","Cd1d1","Il1a","Car11")
vf_calcium <- c("Bhlha15","Ccr4","Ppp3cc","Ccr3","Vcam1","Lat","Ank2","Cxcr5","Cxcr3","Ppp1r9a","Tnfsf11","Ccr6.2","Nr5a1","Xcr1","Ksr2","Dmd","Ptgdr2","Ryr2")

INPUT="C:/Users/David/Desktop/Sequencing/AD005/jak2_raw-counts.txt"
OUTPUT="C:/Users/David/Desktop/Sequencing/AD005/il2.svg"
OUTPUT2="C:/Users/David/Desktop/Sequencing/AD005/calcium.svg"
#GROUPS=c("LPS","HDL.LPS","CL.LPS","CL.HDL.LPS")
#GROUPS=c("Veh","LPS","T0L")
#GROUPS=c("BMDM 10ngLPS","Lesional Mac (LCM)")
#GROUPS=c("WT_CtrlASO","WT_T39ASO","KO_CtrlASO","KO_T39ASO")
GROUPS=c("wt_wt","vf_wt","wt_ldlr","vf_ldlr")
COMP=c(1,2)
REPL=c(6,7,6,4)

# calculate RPKM and FC_RPKM
data<-read.delim(INPUT,header=TRUE,stringsAsFactors=TRUE)
rownames(data) = make.names(data$Gene, unique=TRUE)
bp <- data[,6]
counts <- data[,-c(1:6)]
rpm <- sweep(counts,2,colSums(counts),'/')*1e6
rpkm <- sweep(rpm,1,bp,'/')*1e3

r1 <- (sum(REPL[1:COMP[1]-1])+1):(sum(REPL[1:COMP[1]]))
r2 <- (sum(REPL[1:COMP[2]-1])+1):(sum(REPL[1:COMP[2]]))
pd1 <- log2(data.matrix(rpkm[row.names(rpkm) %in% vf_il2, ])+0.001)
#rownames(pd1)[rownames(pd1)=="Gm4951"]="Ifgga2"
# rownames(pd1)[rownames(pd1)=="Gm12185"]="Irgb2b1"
# rownames(pd1)[rownames(pd1)=="Gm4841"]="Ifgga3"
#rownames(pd1)[rownames(pd1)=="F830016B08Rik"]="Ifgga4"
# rownames(pd1)[rownames(pd1)=="Gbp3.1"]="Gbp3"
pd1 <- pd1[,c(r1,r2)]

svg(OUTPUT, width=20,height=4)
hmcol <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(100))
h <- heatmap.2(pd1,Rowv=TRUE,Colv=FALSE,
               distfun = function(x) as.dist(1-cor(t(pd1))),
               scale="row", trace="none", dendrogram="none",
               col=hmcol, labCol=FALSE,density.info="none",
               lhei=c(1,3), lwid=c(1,2.5),cexRow=2.2,
               breaks=(0:100-50)/(100/4),
               key.par=list(cex=0.6),margins=c(2,12))
dev.off()

pd2 <- log2(data.matrix(rpkm[row.names(rpkm) %in% vf_calcium, ])+0.001)
pd2 <- pd2[,c(r1,r2)]

svg(OUTPUT2, width=16,height=6)
hmcol <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(100))
h <- heatmap.2(pd2,Rowv=TRUE,Colv=FALSE,
               distfun = function(x) as.dist(1-cor(t(pd2))),
               scale="row", trace="none", dendrogram="none",
                col=hmcol, labCol=FALSE,density.info="none",
                lhei=c(1,3), lwid=c(1,2.5),cexRow=1.8,
                breaks=(0:100-50)/(100/4),
                key.par=list(cex=0.6),margins=c(2,8))
dev.off()