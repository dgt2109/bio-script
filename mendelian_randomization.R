library(MendelianRandomization)
library(matrixStats)

# set parameters
dataset="FDR" # set dataset to "GW" = genome-wide sig. or "FDR" = FDR-sig.
EXCLUDE_WEAK = TRUE
PLOTTING = FALSE
TABLES = TRUE

# outlier lists derived from MR-PRESSO runs (too time consuming to run here)
outlier_GW_HDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
outlier_GW_LDL <- -c(2,3,4,7,15,23,33,34,41,54,65)
outlier_GW_TG <- -c(3,4,7,15,23,33,34,41,54,65)

outlier_FDR_HDL <- -c(2,3,4,5,8,9,17,18,27,42)
outlier_FDR_LDL <- -c(2,3,4,5,8,9,17,18,27,42)
outlier_FDR_TG <- -c(2,3,4,5,8,9,17,18,27,42,55)

lipid_list=c("LDL-C","HDL-C","TG")
regression_params <- data.frame(matrix(ncol=12))
colnames(regression_params) <- c("IVW beta","IVW CILower","IVW CIUpper","IVW P",
                                 "Egger beta","Egger CILower","Egger CIUpper","Egger P",
                                 "Intercept", "Intercept CILower","Intercept CIUpper",
                                 "Intercept P")

for (i in 1:3) {
  for (j in 1:2) {
  
    lipid=lipid_list[i]
    
    # set input and exclusions per parameters
    if(dataset == "FDR") {
      base_dir="C:/Users/David/Desktop/genetics/data/FDR/"
      if(lipid == "HDL-C") { outlist <- outlier_FDR_HDL
      } else if(lipid == "LDL-C") { outlist <- outlier_FDR_LDL
      } else if(lipid == "TG") { outlist <- outlier_FDR_TG }
      table_id="table_547"
      point_thickness=0.35
    } else if (dataset == "GW") {
      base_dir="C:/Users/David/Desktop/genetics/data/GW/"
      if(lipid == "HDL-C") { outlist <- outlier_GW_HDL
      } else if(lipid == "LDL-C") { outlist <- outlier_GW_LDL
      } else if(lipid == "TG") { outlist <- outlier_GW_TG }
      table_id="table_179"
      point_thickness=0.5
    }
    
    # data input
    input_data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
    input_data=input_data[input_data$Common=="TRUE", ] # exclude rare & low-freq
    
    # set output filenames
    weak_app=""
    if(EXCLUDE_WEAK) { weak_app="_noweak" }
    regression_summary=paste(base_dir,table_id,weak_app,".csv",sep="")
    
    # on second iteration, remove MR-PRESSO outliers from outlist
    if(j == 2) { data <- input_data[outlist, ] } else { data <- input_data}
    
    # if parameter set, exclude variants with first-stage F-statistic < 11
    if(EXCLUDE_WEAK) {
      # calculate F-statistics
      F = cbind(input_data$LDL_beta^2/input_data$LDL_se^2, 
                input_data$HDL_beta^2/input_data$HDL_se^2, 
                input_data$TG_beta^2/input_data$TG_se^2)
      Fmax=rowMaxs(F) 
      if(length(Fmax[Fmax<11]) > 0) {
        input_data <- input_data[-match(Fmax[Fmax<11],Fmax), ]
      }
    }
    
    ###################################################################
    
    # set up intial MV-IVW on lipid covariates for adjustment
    if(lipid == "LDL-C") {
      data1=cbind(data$HDL_beta,data$HDL_se)
      data2=cbind(data$TG_beta,data$TG_se)
      data3=cbind(data$LDL_beta,data$LDL_se)
      adj_cov="HDL-C- and TG-"
    } else if(lipid == "HDL-C") {
      data1=cbind(data$LDL_beta,data$LDL_se)
      data2=cbind(data$TG_beta,data$TG_se)
      data3=cbind(data$HDL_beta,data$HDL_se)
      adj_cov="LDL-C- and TG-"
    } else if(lipid == "TG") {
      data1=cbind(data$LDL_beta,data$LDL_se)
      data2=cbind(data$HDL_beta,data$HDL_se)
      data3=cbind(data$TG_beta,data$TG_se)
      adj_cov="LDL-C- and HDL-C-"
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

    adjusted_y[adjusted_x < 0] <- -adjusted_y[adjusted_x < 0]
    adjusted_x[adjusted_x < 0] <- -adjusted_x[adjusted_x < 0]
    
    # set up IVW and EGGER objects on lipid covariate adjusted CAD
    MRInputObject <- mr_input(bx = adjusted_x,
                              bxse = xerr,
                              by = adjusted_y,
                              byse = yerr,
                              exposure=lipid,
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
    
    # regression plots
    if(PLOTTING && j == 2 && EXCLUDE_WEAK == FALSE) {
      plot_file=paste(base_dir,table_id,lipid,".png",sep="")
      png(plot_file,width=6,height=4,units="in",res=600)
      plot(adjusted_x,adjusted_y,
           xlab=paste("SD change in ",lipid," per allele",sep=""),
           ylab=paste(adj_cov,"adjusted logOR CAD",sep=""),
           main=paste(lipid," vs. lipid-adjusted CAD risk",sep=""),
           pch=20,cex=point_thickness)
      abline(0,IVWObject$Estimate,lty=3)
      abline(h=0)
      dev.off()
    }  
    
    iter_stats <- c(exp(IVWObject$Estimate), exp(IVWObject$CILower), 
                    exp(IVWObject$CIUpper), IVWObject$Pvalue,
                    exp(EggerObject$Estimate), exp(EggerObject$CILower.Est),
                    exp(EggerObject$CIUpper.Est), EggerObject$Pvalue.Est,
                    EggerObject$Intercept, EggerObject$CILower.Int,
                    EggerObject$CIUpper.Int, EggerObject$Pvalue.Int)
    regression_params <- rbind(regression_params,iter_stats)
  }
}

if (TABLES) {
  regression_params <- regression_params[-1,]
  rownames(regression_params) <- cbind("LDL-C","LDL-C MR-PRESSO",
                                       "HDL-C","HDL-C MR-PRESSO",
                                       "TG", "TG MR-PRESSO")
  write.csv(regression_params,file=regression_summary)
}

