library("ggpubr")
sim_output <- "C:/Users/David/Desktop/genetics/revision/plot/sim_"

# sim fit output: a (instrument strength vector), rf_se (instrument error vector)
str_label <- c(rep("Simulated",length(a)), rep("Observed",length(data$HDL_beta)))
instr_str <- c(a,abs(data$HDL_beta))
simfit_str <- data.frame(instr_str,str_label)
svg(paste0(sim_output,"fit_1.svg"),width=4,height=3)
p <- gghistogram(simfit_str,x="instr_str",y="..density..",fill="str_label",bins=15,binwidth=0.01)
p <- ggpar(p, legend.title="",xlab="Instrument Strength",ylab="Frequency")
print(p)
dev.off()

se_label <- c(rep("Simulated",length(rf_se)),rep("Observed",length(data$HDL_se)))
instr_se <- c(rf_se,data$HDL_se)
simfit_se <- data.frame(instr_se,se_label)
svg(paste0(sim_output,"fit_2.svg"),width=4,height=3)
p2 <- gghistogram(simfit_se,x="instr_se",y="..density..",fill="se_label",bins=10,binwidth=0.0015)
p2 <- ggpar(p2, legend.title="",xlab="Instrument Standard Error",ylab="Frequency")
print(p2)
dev.off()

# rptsim output: TPR, FPR, CRCI

for (i in 1:2) { # iterate over variation in oY vs U
  if (i==1) { param_list <- c("0.5","0.6","0.7")
    xaxis_label <- "Measurement Error" }
  if (i==2) { param_list <- c("0.0","0.1","0.2")
    xaxis_label <- "Confounder Strength" }
  param_labels <- rep(param_list,1,each=3)
  
  ri <- 1:3 + (i-1)*3
  rates_list <- list(TPR[ri,], FPR[ri,], CRCI[ri,])
  beta_list <- sim_betas[ri]
  for (j in 1:2) { # iterate over plotting for IVW vs EGGER
    cols <- c(1,2:3+(j-1)*2)
    if (j == 1) { conditions <- c("Obs.","MR-IVW,Bonf.","MR-IVW,FDR")
    regr_palette = c("#1f78b4","#ff7f00","#6a3d9a") }
    if (j == 2) { conditions <- c("Obs.","MR-Egger,Bonf.","MR-Egger,FDR") 
    regr_palette = c("#1f78b4","#33a02c","#e31a1c") }
    
    for (k in 1:3) { # iterate over rate parameter to plot
      if (k == 1) { yaxis_label <- "True Positive Rate" }
      if (k == 2) { yaxis_label <- "False Positive Rate" }
      if (k == 3) { yaxis_label <- "Coverage Rate of CI" }
    
      rate_matrix <- rates_list[[k]][,cols]
      rate_vector <- c(rate_matrix[1,],rate_matrix[2,],rate_matrix[3,])
      regr_type <- rep(conditions,3)
      rate_data <- data.frame(param_labels,regr_type,rate_vector)
      rate_data$regr_type <- factor(rate_data$regr_type,levels=conditions)
      svg(paste0(sim_output,k+(j-1)*3+(i-1)*6,".svg"),width=4,height=3)
      p3 <- ggbarplot(rate_data, x="param_labels", y="rate_vector", fill="regr_type",
                      position=position_dodge(0.7),palette=regr_palette)
      p3 <- ggpar(p3, legend.title="",xlab=xaxis_label,ylab=yaxis_label,ylim=c(0,1))
      print(p3)
      dev.off()
    }
    beta_vector <- c(beta_list[[1]][,cols][,1],beta_list[[1]][,cols][,2],beta_list[[1]][,cols][,3],
                      beta_list[[2]][,cols][,1],beta_list[[2]][,cols][,2],beta_list[[2]][,cols][,3],
                      beta_list[[3]][,cols][,1],beta_list[[3]][,cols][,2],beta_list[[3]][,cols][,3])
    beta_label1 <- c(rep(conditions,3,each=250))
    beta_label2 <- c(rep(param_list,1,each=750))
    simbetas_data <- data.frame(beta_label1,beta_label2,beta_vector)
    simbetas_data$beta_label1 <- factor(simbetas_data$beta_label1,levels=conditions)
    svg(paste0(sim_output,12+j+(i-1)*2,".svg"),width=4,height=3)
    p4 <- ggboxplot(simbetas_data, x="beta_label2", y="beta_vector", fill="beta_label1",
                    palette=regr_palette)
    p4 <- ggpar(p4, legend.title="",xlab=xaxis_label,ylab="Estimated Effect Size")
    dev.off()
  }
}    
