library(MendelianRandomization)
library(matrixStats)

# set parameters
dataset="FDR" # set dataset to "GW" = genome-wide sig. or "FDR" = FDR-sig.
TABLES = TRUE

# outlier lists derived from MR-PRESSO runs (too time consuming to run here)
outlier_GW_HDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
outlier_GW_LDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
outlier_GW_TG <- -c(3,4,7,15,23,33,34,41,54,65)

outlier_FDR_HDL <- -c(2,3,4,5,8,9,17,18,27,42)
outlier_FDR_LDL <- -c(2,3,4,5,8,9,17,18,27,42)
outlier_FDR_TG <- -c(2,3,4,5,8,9,17,18,27,42,55)

outlierGW_superset <- unique(c(outlier_GW_HDL,outlier_GW_LDL,outlier_GW_TG))
outlierFDR_superset <- unique(c(outlier_FDR_HDL,outlier_FDR_LDL,outlier_FDR_TG))

# set input and exclusions per parameters
if(dataset == "FDR") {
  base_dir="C:/Users/David/Desktop/genetics/data/FDR/"
  outlist <- outlierFDR_superset
  table_id="table_547"
} else if (dataset == "GW") {
  base_dir="C:/Users/David/Desktop/genetics/data/GW/"
  outlist <- outlierGW_superset
  table_id="table_179"
}

# set output filenames
mvivw_summary <- paste(base_dir,table_id,"_mvivw",".csv",sep="")

regression_params <- data.frame(matrix(ncol=4))
colnames(regression_params) <- c("MVIVW beta", "CILower", "CIUpper", "MVIVW P")
for (j in 1:2) {
  
  # data input
  data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
  data=data[data$Common=="TRUE", ]
  if(j == 2) { data <- data[outlist, ] }
  
  # set up MVIVW object
  X <- cbind(data$LDL_beta,data$HDL_beta,data$TG_beta)
  Xerr <- cbind(data$LDL_se,data$HDL_se,data$TG_se)
  z <- data$CAD_beta
  zerr <- data$CAD_se
  
  MRMVInputObject <- mr_mvinput(bx = X,
                                 bxse = Xerr,
                                 by = z,
                                 byse = zerr,
                                 exposure = c("LDL-C", "HDL-C", "TG"),
                                 outcome = "CAD")
  
  MRMVObject <- mr_mvivw(MRMVInputObject)
  iter_data <- cbind(exp(MRMVObject$Estimate), 
                     exp(MRMVObject$CILower),
                     exp(MRMVObject$CIUpper),
                     MRMVObject$Pvalue)
  
  colnames(iter_data) <- c("MVIVW beta", "CILower", "CIUpper", "MVIVW P")
  rownames(iter_data) <- c("LDL-C","HDL-C","TG")
  regression_params <- rbind(regression_params,iter_data)
}

# regression summary output
if(TABLES) {
  regression_params <- regression_params[-1,]
  rownames(regression_params) <- cbind("LDL-C","HDL-C","TG",
                                       "LDL-C MR-PRESSO","HDL-C MR-PRESSO",
                                       "TG MR-PRESSO")
  write.csv(regression_params,file=mvivw_summary)
}