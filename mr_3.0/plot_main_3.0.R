library("ggplot2")

output = "C:/Users/David/Desktop/genetics/revision/dotplot_"

PARTIAL = F
PRESSO = T
for (k in 1:2) {
  if (k == 2) { 
    PARTIAL = T 
    part = which(pmin(covs[[1]]$LDL_p,covs[[1]]$HDL_p,covs[[1]]$TG_p) < 5e-8) 
  }
  source("mvmr_3.0.R")
  
  for (m in 1:3) {
    cov <- covs[[m]]
    ocov <- ocovs[[m]][-which(ocovs[[m]]$rsid == "rs2315065"),]
    allcov <- rbind(cov,ocov)
    
    bi <- which(colnames(cov) %in% c("LDL_beta","HDL_beta","TG_beta"))
    bai <- bi[-m]
    ci <- c(1:3)[-m]
    cov$adj_CAD <- cov$CAD_beta - (cov[,bai[1]] * ivw$coef[ci[1]]) - 
                                  (cov[,bai[2]] * ivw$coef[ci[2]])
    ocov$adj_CAD <- ocov$CAD_beta - (ocov[,bai[1]] * ivw$coef[ci[1]]) - 
                                     (ocov[,bai[2]] * ivw$coef[ci[2]])
    allcov$adj_CAD <- allcov$CAD_beta - (allcov[,bai[1]] * ivw$coef[ci[1]]) - 
                                        (allcov[,bai[2]] * ivw$coef[ci[2]])
    allcov$CAD_pred <- allcov[,bi[m]]*ivw$coef[m]
    allcov$CAD_egger <- allcov[,bi[m]]*eggers[[m]]$coef[m+1]
    ivw_res <- summary(ivw)
    allcov$CAD_lb <- allcov[,bi[m]] * 
      (ivw$coef[m]-qt(0.975,df=nrow(cov)-3)*ivw_res$coef[m,2])
    allcov$CAD_ub <- allcov[,bi[m]] * 
      (ivw$coef[m]+qt(0.975,df=nrow(cov)-3)*ivw_res$coef[m,2])
    allcov$xaxis <- rep(0,nrow(allcov))
    colnames(allcov)[bi[m]] <- "x"
    colnames(cov)[bi[m]] <- "x"
    colnames(ocov)[bi[m]] <- "x"
    
    svg(paste0(output,m+(k-1)*3,".svg"),width=5,height=2.5)
    p <- ggplot(cov, aes(x=x,y=adj_CAD)) + geom_point(size=0.2*i) +
      geom_point(ocov, mapping=aes(x=x,y=adj_CAD),size=1.5,shape=21) +
      geom_line(data=allcov,aes(x=x,y=CAD_pred, color="MVMR-IVW"),linetype=2) +
      geom_line(data=allcov,aes(x=x,y=CAD_egger,color="MVMR-EGGER"),linetype=2) +
      scale_color_manual(name="",values=c("MVMR-IVW"="black",
                                          "MVMR-EGGER"="red"),
                         guide=guide_legend(reverse=T)) +
      geom_line(data=allcov,color='black',aes(x=x,y=xaxis)) +
      geom_ribbon(data=allcov, mapping=aes(x=x,ymin=CAD_lb,
                                           ymax=CAD_ub),alpha=0.2) +
      scale_x_continuous(expand=c(0.001,0.001)) +
      xlab(paste0(c("LDL","HDL","TG")[m]," effect size")) +
      ylab("Adjusted CAD effect size") +
      theme_bw() + theme(panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.margin=margin(5,10,5,10),
                         legend.background=element_blank(),
                         legend.position=c(0.85,0.85 - (m%%2)*0.6))
  #   legend.position=c(0.85,0.25 for bottomright)
    print(p)
    dev.off()
  }
}

