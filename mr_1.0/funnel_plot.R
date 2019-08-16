library(MendelianRandomization)
outlier185 <- c(-2,-3,-4,-7,-15,-33,-34,-41,-54,-65)
outlier563 <- c(-2,-3,-4,-5,-8,-9,-17,-18,-27,-42,-55,-86)

base_dir="C:/Users/David/Desktop/genetics/data/"
table_id="563/table_547"
exlist <- outlier563

data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
PRESSOHDL=TRUE

data=data[data$Common=="TRUE", ]
if(PRESSOHDL) {
  data <- data[exlist, ]
}

###################################################################

W <- cbind(data$LDL_beta,data$TG_beta)
err <- cbind(data$LDL_se,data$TG_se)
h <- data$HDL_beta
herr <- data$HDL_se
y <- data$CAD_beta
yerr <- data$CAD_se

MRMVInputObject2 <- mr_mvinput(bx = W,
                               bxse = err,
                               by = y,
                               byse = yerr,
                               exposure = c("LDL-C", "TG"),
                               outcome = "CAD")

MRMVObject2 <- mr_mvivw(MRMVInputObject2)

##################################################################

adjusted_CAD <- y - MRMVObject2$Estimate[1]*W[,1] - MRMVObject2$Estimate[2]*W[,2]
#adjusted_CAD[h < 0] <- -adjusted_CAD[h < 0]
#h[h < 0] <- -h[h < 0]

vaf <- data$VAF
vaf[vaf<0.5] <- 1 - vaf[vaf<0.5]
vaf_correction <- lm(h ~ vaf)

adjusted_h <- h - vaf_correction$coefficients[1] - vaf_correction$coefficients[2]*vaf
adjusted_CAD[adjusted_h < 0] <- -adjusted_CAD[adjusted_h < 0]
adjusted_h[adjusted_h < 0] <- -adjusted_h[adjusted_h < 0]

plot(adjusted_CAD/h,adjusted_h,xlim=c(-20,20),
     xlab="Indiv. variant HDL-C:CAD assoc.",
     ylab="Indiv. variant precision estimate",
     main="563 variants",cex=0.5)
#abline(v=median(adjusted_CAD/h))