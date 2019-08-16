
INPUT <- "C:/Users/David/Desktop/Sequencing/AD002/atac_qc_metrics/atac-counts_macpeak-input.txt"
PLOT1 <- "C:/Users/David/Desktop/Sequencing/AD002/atac_qc_metrics/atac-counts_macpeak_Ctrl.tiff"
PLOT2 <- "C:/Users/David/Desktop/Sequencing/AD002/atac_qc_metrics/atac-counts_macpeak_LPS.tiff"

data <- read.delim(INPUT,header=TRUE)
data2<-data[-c(which(data$ATC_D==0),which(data$ATC_DL==0),which(data$X27_K==0),  
               which(data$X27_NT==0),which(data$X4_K==0),which(data$X4_NT==0)),]
#data2<-data[-c(which(data$ATC_D==0),which(data$X27_NT==0)),]

a <- log(data2$ATC_D)
b <- log(data2$X27_NT)
c <- log(data2$ATC_DL)
d <- log(data2$X27_K)
 
tiff(PLOT1, width=4,height=4, units="in", res=600)
plot(a,b,pch=20,cex=0.1,xlab="log ATAC counts/peak Ctrl",ylab="log H3K27ac counts/peak Ctrl",
     main="ATAC vs H3K27ac signal, Ctrl")
dev.off()

setEPS()
tiff(PLOT2, width=4,height=4, units="in", res=600)
plot(c,d,pch=20,cex=0.1,xlab="log ATAC counts/peak LPS",ylab="log H3K27ac counts/peak KLA",
     main="ATAC vs H3K27ac signal, TLR4-stim.")
dev.off()

#res <- lm(y~x+0)
#abline(0,res$coefficients[1],lty=2,col=2)