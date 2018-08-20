library(MendelianRandomization)
library(MRPRESSO)

# these variables store results of the analysis below
#outlier_GW_HDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
#outlier_GW_LDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
#outlier_GW_TG <- -c(3,4,7,15,23,33,34,41,54,65)

#outlier_FDR_HDL <- -c(2,3,4,5,8,9,17,18,27,42)
#outlier_FDR_LDL <- -c(2,3,4,5,8,9,17,18,27,42)
#outlier_FDR_TG <- -c(2,3,4,5,8,9,17,18,27,42,55)

datasets=c("GW","FDR")
lipid_list=c("LDL","HDL","TG")

for (i in 1:2) {
  dataset=datasets[i]
  
  for (j in 1:3) {
    lipid=lipid_list[j]
    print(paste(dataset," + ",lipid,sep=""))
    
    # i/o
    if(dataset=="GW") {
      base_dir="C:/Users/David/Desktop/genetics/data/GW/"
      exclusion_file=paste("C:/Users/David/Desktop/genetics/data/GW/",
                           lipid,"_",dataset,"_exlist.txt",sep="")
      table_id="table_179"
      Nsim=10000
    } else if(dataset=="FDR") {
      base_dir="C:/Users/David/Desktop/genetics/data/FDR/"
      exclusion_file=paste("C:/Users/David/Desktop/genetics/data/FDR/",
                           lipid,"_",dataset,"_exlist.txt",sep="")
      table_id="table_547"
      Nsim=100000
    }
    
    data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
    data=data[data$Common=="TRUE", ]
    
    # covariate adjustment
    if(lipid == "LDL") {
      data1=cbind(data$HDL_beta,data$HDL_se)
      data2=cbind(data$TG_beta,data$TG_se)
      data3=cbind(data$LDL_beta,data$LDL_se)
    } else if(lipid == "HDL") {
      data1=cbind(data$LDL_beta,data$LDL_se)
      data2=cbind(data$TG_beta,data$TG_se)
      data3=cbind(data$HDL_beta,data$HDL_se)
    } else if(lipid == "TG") {
      data1=cbind(data$LDL_beta,data$LDL_se)
      data2=cbind(data$HDL_beta,data$HDL_se)
      data3=cbind(data$TG_beta,data$TG_se)
    }
    
    W <- cbind(data1[,1],data2[,1])
    err <- cbind(data1[,2],data2[,2])
    
    x <- data3[,1]
    xerr <- data3[,2]
    
    y <- data$CAD_beta
    yerr <- data$CAD_se
    
    MRMVInputObject_y <- mr_mvinput(bx = W,
                                    bxse = err,
                                    by = y,
                                    byse = yerr)
    MRMVInputObject_x <- mr_mvinput(bx = W, 
                                    bxse = err, 
                                    by = x, 
                                    byse = xerr)
    
    MRMVObject_y<- mr_mvivw(MRMVInputObject_y)
    MRMVObject_x<- mr_mvivw(MRMVInputObject_x)
    
    adjusted_y <- y - MRMVObject_y$Estimate[1]*W[,1] - MRMVObject_y$Estimate[2]*W[,2]
    adjusted_x <- x - MRMVObject_x$Estimate[1]*W[,1] - MRMVObject_x$Estimate[2]*W[,2]
    
    # run MR-PRESSO
    SummaryStats <- data.frame(cbind(adjusted_x,xerr,adjusted_y,yerr))
    colnames(SummaryStats) <- c(lipid,paste(lipid,"_err",sep=""),"CAD","CAD_err")
    res_presso <- mr_presso(BetaOutcome = "CAD",
                            BetaExposure = lipid,
                            SdOutcome = "CAD_err",
                            SdExposure = paste(lipid,"_err",sep=""),
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE,
                            data = SummaryStats,
                            NbDistribution = Nsim, 
                            SignifThreshold = 0.05)
    ot <- res_presso$"Outlier Test"
    ex <- ot[ot$Pvalue < 0.05, ]
    
    write.csv(ex,exclusion_file)
  }
}