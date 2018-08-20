library(MendelianRandomization)
library(matrixStats)

# i/o
base_dir="C:/Users/David/Desktop/genetics/data/FDR/"
table_id="table_547"
base_data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)

# set outlier lists from running calc_MRPRESSO.R
outlier_FDR_HDL <- -c(2,3,4,5,8,9,17,18,27,42)
exlist <- outlier_FDR_HDL
PRESSOHDL=TRUE
base_data=base_data[base_data$Common=="TRUE", ]
if(PRESSOHDL) {
  base_data <- base_data[exlist, ]
}

# calculate t-values for each first-stage HDL-C association
T <- cbind(base_data$LDL_beta/base_data$LDL_se,
           base_data$HDL_beta/base_data$HDL_se,
           base_data$TG_beta/base_data$TG_se)
thdl <- base_data$HDL_beta/base_data$HDL_se
tmax <- rowMaxs(T)
data_t <- cbind(base_data,thdl)
dorder <- data_t[order(-abs(thdl)), ]

###################################################################

track <- data.frame(matrix(ncol=8))
colnames(track) <- c("IVW beta","IVW se","IVW P",
                     "Egger beta","Egger se", "Egger P",
                     "Egger intercept", "Egger intercept P")

# iterate over expanding variant sets and calculate MR parameters
for (i in 1:511) {
  data <- head(dorder,i)
  
  W <- cbind(data$LDL_beta,data$TG_beta)
  err <- cbind(data$LDL_se,data$TG_se)
  h <- data$HDL_beta
  herr <- data$HDL_se
  y <- data$CAD_beta
  yerr <- data$CAD_se
  
  MRMVInputObject_y <- mr_mvinput(bx = W,
                                  bxse = err,
                                  by = y,
                                  byse = yerr)
  MRMVInputObject_h <- mr_mvinput(bx = W, 
                                  bxse = err, 
                                  by = h, 
                                  byse = herr)
  
  MRMVObject_y<- mr_mvivw(MRMVInputObject_y)
  MRMVObject_h<- mr_mvivw(MRMVInputObject_h)
  
  adjusted_y <- y - MRMVObject_y$Estimate[1]*W[,1] - MRMVObject_y$Estimate[2]*W[,2]
  adjusted_h <- h - MRMVObject_h$Estimate[1]*W[,1] - MRMVObject_h$Estimate[2]*W[,2]
  
  adjusted_y[adjusted_h < 0] <- -adjusted_y[adjusted_h < 0]
  adjusted_h[adjusted_h < 0] <- -adjusted_h[adjusted_h < 0]
  
  # set up IVW and EGGER objects on lipid covariate adjusted CAD
  MRInputObject <- mr_input(bx = adjusted_h,
                            bxse = herr,
                            by = adjusted_y,
                            byse = yerr,
                            exposure="HDL-C",
                            outcome="CAD")
  
  IVWObject <- mr_ivw(MRInputObject,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple",
                      psi = 0,
                      distribution = "normal",
                      alpha = 0.05)
  
  
  EggerObject <- mr_egger(MRInputObject,
                          robust = FALSE,
                          penalized = FALSE,
                          correl = FALSE,
                          distribution = "normal",
                          alpha = 0.05)
  
  iter_stats <- c(IVWObject$Estimate, IVWObject$StdError, 
                  IVWObject$Pvalue, EggerObject$Estimate,
                  EggerObject$StdError.Est, EggerObject$Pvalue.Est,
                  EggerObject$Intercept, EggerObject$Pvalue.Int)
  track <- rbind(track,iter_stats)
}
track <- track[-1,]

# i/o
se_file=paste(base_dir,table_id,"subsets_se.png",sep="")
p_file=paste(base_dir,table_id,"subsets_p.png",sep="")
beta_file=paste(base_dir,table_id,"subsets_beta.png",sep="")
intercept_file=paste(base_dir,table_id,"subsets_int.png",sep="")
plotheight=3
plotwidth=6
altcol=4

# plot results of iterated regressions
png(se_file,width=plotwidth,height=plotheight,units="in",res=600)
plot(track$`Egger se`,pch=2,col=altcol,xlab="# Variants, ordered by first-stage P-value",
     ylab="HDL-C coefficient SE",main="Standard Error",
     ylim=c(0,0.15),cex=0.5)
points(track$`IVW se`,pch=1,col=1,cex=0.5)
#abline(h=-tail(track$`Egger beta`,1)/2,lty=3)
legend("topright",legend=c("MR-IVW","MR-EGGER"),col=c(1,altcol),pch=c(1,2))
dev.off() 

png(p_file,width=plotwidth,height=plotheight,units="in",res=600)
plot(log(track$`Egger P`),pch=2,col=altcol,xlab="# Variants, ordered by first-stage P-value",
     ylab="log HDL-C coefficient P",main="Log P-value",
     ylim=c(-9,1),cex=0.5)
abline(h=log(0.05),lty=3)
points(log(track$`IVW P`),pch=1,col=1,cex=0.5)
legend("bottomleft",legend=c("MR-IVW","MR-EGGER"),col=c(1,altcol),pch=c(1,2))
dev.off() 

png(beta_file,width=plotwidth,height=plotheight,units="in",res=600)
plot(track$`Egger beta`,pch=2,col=altcol,xlab="# Variants, ordered by first-stage P-value",
     ylab="HDL-C coefficient estimate",main="HDL-C Coefficient",
     ylim=c(-0.15,0.15),cex=0.5)
points(track$`IVW beta`,pch=1,col=1,cex=0.5)
legend("topright",legend=c("MR-IVW","MR-EGGER"),col=c(1,altcol),pch=c(1,2))
dev.off()

png(intercept_file,width=plotwidth,height=plotheight,units="in",res=600)
plot(track$`Egger intercept`,pch=2,col=altcol,xlab="# Variants, ordered by first-stage P-value",
     ylab="MR-EGGER intercept estimate",main="MR-EGGER Intercept",cex=0.5,
     ylim=c(-0.015,0.002))
abline(h=0)
dev.off()