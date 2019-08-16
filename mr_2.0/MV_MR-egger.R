library("ggpubr")

base_dir="C:/Users/David/Desktop/genetics/revision/snp680/"
file_in="table_680.txt"
data <- read.delim(paste0(base_dir, file_in),header=TRUE)
#pleiotropic <- c(46,377,380,429,529) # MR-PRESSO FDR
#PRESSO_SNPs rs653178 rs934287 rs1250229 rs2315065 rs12940887 

#psize = 0.2
#plot_app = "_FDR2"
#filename = "C:/Users/david/Desktop/FDR_coef.csv"

# base_dir="C:/Users/David/Desktop/genetics/revision/snp680/"
# file_in="table_179_common.txt"
# pleiotropic <- c(4,8) # MR-PRESSO GW
# psize = 0.4
# plot_app = "_GW2"

#pleiotropic <- c(31,46,350,377,380,429,458,529,588) # Cochran's Q/Bonf
#pleiotropic <- c(11,23,31,32,34,46,47,49,68,73,83,84,96,99,102,
#                  107,109,110,113,121,123,127,128,138,139,155,159,
#                  165,170,174,177,180,181,190,195,205,206,214,215,
#                  217,219,222,227,231,236,238,243,245,257,259,266,272,
#                  290,291,294,301,303,305,307,308,309,311,315,322,329,
#                  346,350,375,377,380,386,387,390,393,413,415,417,421,
#                  425,427,428,429,437,446,458,471,478,479,484,490,493,
#                  497,500,526,529,534,537,546,552,558,559,581,588,597,
#                  614,616,621,623,631,636) # Cochran's Q, uncorrected
#bmi_rs <- as.vector(unlist(read.delim("C:/Users/david/Desktop/genetics/rerevision/bmi_q5_886.txt",
#                          header=FALSE)))
#presso_rs <- c("rs653178","rs934287","rs1250229","rs2315065","rs12940887")
#ex_rs <- c(bmi_rs) # + presso_rs
#ex_ind <- as.vector(unlist(sapply(ex_rs, FUN = function(x) which(data[1] == x))))
#pleiotropic<-c(ex_ind)
#nvars = 631

lname = c("LDL","TG","HDL")
plotname = c("LDL-C","TG","HDL-C")
lipid = paste0("data$",lname,"_beta")
olipid = paste0("outlier$",lname,"_beta")
err = paste0("data$",lname,"_se")
record=data.frame(row.names=c("Beta_IVW","SE_IVW","OR_IVW",
                              "CILower_IVW","CIUpper_IVW","P_IVW","I2IVW",
                              "Beta_Egger","OR_Egger","SE_Egger",
                              "CILower_Egger","CIUpper_Egger","P_Egger",
                              "I2_Egger","Int_Egger","Int.P_Egger"))
for (i in 1:length(lipid)) {
  # orientation
  data <- read.delim(paste0(base_dir, file_in),header=TRUE)
  
  clist = c(lipid[i%%3+1],lipid[(i+1)%%3+1],"data$CAD_beta")
  for (var in clist) {
    eval(parse(text=paste0(var,"<-ifelse(",lipid[i],">0,",
                           var,",",var,"*-1)")))
  }
  eval(parse(text=paste0(lipid[i],"<-","abs(",lipid[i],")")))
  
# code for exclusion of pleiotropic variant list
  
  outlier <- c()#data[pleiotropic, ]
#  data <- data[-pleiotropic, ]
  
# code for subsetting by first-stage F-statistic
  
#  subset <- c(which((data$HDL_beta/data$HDL_se)^2 >= 100),
#              which((data$HDL_beta/data$HDL_se)^2 < 10))
#  data <- data[-subset, ]

# code for subsetting by # of variants ranked by lipid q

#  data <- data[order(pmin(data$LDL_q,data$HDL_q,data$TG_q)),]
#  data <- head(data,nvars)
  
  # ivw & egger regression
  eval(parse(text=paste0("ivw <- lm(data$CAD_beta~",lipid[i],"+",
                                 lipid[i%%3+1],"+",lipid[(i+1)%%3+1],
                                 "+0, weights = data$CAD_se^-2)")))
  eval(parse(text=paste0("egger <- lm(data$CAD_beta~",lipid[i],"+",
                                      lipid[i%%3+1],"+",lipid[(i+1)%%3+1],
                                      ",weights = data$CAD_se^-2)")))
  
  # ivw coefficients
  ivw_res <- summary(ivw)
  beta_ivw <- ivw_res$coef[1]
  se_ivw <- ivw_res$coef[1,2] / min(ivw_res$sigma,1)
  lb_ivw <- min(beta_ivw + qt(0.975,df=nrow(data)-3)*se_ivw, 
                beta_ivw - qt(0.975,df=nrow(data)-3)*se_ivw)
  ub_ivw <- max(beta_ivw + qt(0.975,df=nrow(data)-3)*se_ivw, 
                beta_ivw - qt(0.975,df=nrow(data)-3)*se_ivw)
  p_ivw <- 2*pt(abs(beta_ivw/se_ivw),df=nrow(data)-3,lower=FALSE)
  Q_ivw <- sum((ivw$resid)^2*ivw$weights)
  h2_ivw <- Q_ivw/(nrow(data)-3)
  i2h_ivw <- ((h2_ivw-1)/h2_ivw)*100
  p_Qivw <- pchisq(Q_ivw,df=nrow(data)-3,lower.tail=FALSE)
  snpwise_Q <- (ivw$resid)^2*ivw$weights
  snpwise_het <- pchisq((ivw$resid)^2*ivw$weights,df=1,lower.tail=FALSE)
  het_snps <- which(snpwise_het < 0.05/nrow(data))
  
  # egger coefficients
  egger_res <- summary(egger)
  theta1ME <- egger_res$coef[2]
  se_theta1ME.random <- egger_res$coef[2,2] / min(egger_res$sigma,1)
  lb_theta1ME = theta1ME - qt(0.975,df=nrow(data)-4)*se_theta1ME.random
  ub_theta1ME = theta1ME + qt(0.975,df=nrow(data)-4)*se_theta1ME.random
  p_theta1ME = 2*pt(abs(theta1ME/se_theta1ME.random),df=nrow(data)-4,
                  lower=FALSE)
  
  interME = egger_res$coef[1]
  se_interME.random = egger_res$coef[1,2] / min(egger_res$sigma,1)
  p_interME = 2*(1-pt(abs(interME/se_interME.random),df=nrow(data)-4))
  Q_egger <- sum((egger$resid)^2*egger$weights)
  h2_egger <- Q_egger/(nrow(data)-3)
  i2h_egger <- ((h2_egger-1)/h2_egger)*100
  p_Qegger <- pchisq(Q_egger,df=nrow(data)-4,lower.tail=FALSE)
  
  # y-adjustment for plotting
  eval(parse(text=paste0("adj_CAD <- data$CAD_beta - ",lipid[i%%3+1],
                         " * ivw_res$coef[2] - ",lipid[(i+1)%%3+1],
                         " * ivw_res$coef[3]")))
  eval(parse(text=paste0("outlier_CAD <- outlier$CAD_beta - ",
                         olipid[i%%3+1]," * ivw_res$coef[2] - ",
                         olipid[(i+1)%%3+1]," * ivw_res$coef[3]")))
  
  # plotting
  # eval(parse(text=paste0("df_plot <- data.frame(x=",lipid[i],",
  #                        y=adj_CAD)")))
  # eval(parse(text=paste0("df_plot2 <- data.frame(x=",olipid[i],",
  #                        y=outlier_CAD)")))
  # df_plot2 <- df_plot2[-4,] # hide the extreme outlier Lp(a) in plot
  # 
  # df_all <- rbind(df_plot,df_plot2)
  # df_all$y_ivw <- df_all$x*beta_ivw
  # df_all$y_egger <- df_all$x*theta1ME + interME
  # df_all$y_lb <- df_all$x*lb_ivw
  # df_all$y_ub <- df_all$x*ub_ivw
  # df_all$xaxis <- rep(0,nrow(df_all))
  # 
  # svg(paste0(base_dir,lname[i],plot_app,".svg"),width=5,height=2.5)
  # g <- ggplot(df_plot, aes(x=x,y=y)) + geom_point(size=psize) +
  #   geom_point(df_plot2, mapping=aes(x=x,y=y),size=1.5,shape=21) +
  #   geom_line(data=df_all,aes(x=x,y=y_ivw, color="MVMR-IVW"),linetype=2) +
  #   geom_line(data=df_all,aes(x=x,y=y_egger,color="MVMR-EGGER"),linetype=2) +
  #   scale_color_manual(name="",values=c("MVMR-IVW"="black",
  #                                       "MVMR-EGGER"="red"),
  #                      guide=guide_legend(reverse=T)) +
  #   geom_line(data=df_all,color='black',aes(x=x,y=xaxis)) +
  #   geom_ribbon(data=df_all, mapping=aes(x=x,ymin=y_lb,
  #                                        ymax=y_ub),alpha=0.2) +
  #   scale_x_continuous(expand=c(0.001,0.001)) +
  #   xlab(paste0(plotname[i]," effect size")) +
  #   ylab("Adjusted CAD effect size") +
  #   theme_bw() + theme(panel.grid.major=element_blank(),
  #                      panel.grid.minor=element_blank(),
  #                      plot.margin=margin(5,10,5,10),
  #                      legend.background=element_blank(),
  #                      legend.position=c(0.85,0.85))
  # # legend.position=c(0.85,0.25 for bottomright)
  # print(g)
  # dev.off()
  
  allres <- matrix(c(beta_ivw,se_ivw,exp(beta_ivw),exp(lb_ivw),exp(ub_ivw),
              p_ivw,i2h_ivw,theta1ME,exp(theta1ME),se_theta1ME.random,
              exp(lb_theta1ME),exp(ub_theta1ME),p_theta1ME,i2h_egger,
              interME,p_interME))
  colnames(allres) <- lname[i]
  record <- cbind(record,allres)
  
  # # diagnostics
  # eval(parse(text=paste0("F <- ",lipid[i],"^2/",err[i],"^2")))
  # Q_0.10 <- snpwise_Q[F < 10]
  # Q_10.100 <- snpwise_Q[(F >= 10) & (F < 100)]
  # Q_100.Inf <- snpwise_Q[F >= 100]
  # b1 <- c(Q_0.10,Q_10.100,Q_100.Inf)
  # b1lev <- c(rep(1,length(Q_0.10)),
  #            rep(2,length(Q_10.100)),
  #            rep(3,length(Q_100.Inf)))
  # hetP_0.10 <- snpwise_het[F < 10]
  # hetP_10.100 <- snpwise_het[(F >= 10) & (F < 100)]
  # hetP_100.Inf <- snpwise_het[F >= 100]
  # b2 <- c(hetP_0.10,hetP_10.100,hetP_100.Inf)
  # b2lev <- c(rep("F<10",length(hetP_0.10)),
  #            rep("10<F<100",length(hetP_10.100)),
  #            rep("F>100",length(hetP_100.Inf)))
  # nhet_0.10 <- sum(ifelse(hetP_0.10<0.05,1,0))
  # nhet_10.100 <- sum(ifelse(hetP_10.100<0.05,1,0))
  # nhet_100.Inf <- sum(ifelse(hetP_100.Inf<0.05,1,0))
  # p_het <- c(nhet_0.10/length(hetP_0.10),
  #            nhet_10.100/length(hetP_10.100),
  #            nhet_100.Inf/length(hetP_100.Inf))
  # dfbox <- data.frame(x=b2lev,y=b1)
  # dfbar <- data.frame(x=c("F<10","10<F<100","F>100"),y=p_het*100)
  # 
  # svg(paste0(base_dir,lname[i],"_barplot.svg"),width=2.166,height=2.166)
  # dfbar$x <- factor(dfbar$x,levels=c("F<10","10<F<100","F>100"))
  # p1 <- ggbarplot(dfbar,x="x",y="y") + ylim(0,max(dfbar$y)*1.25) +
  #   ylab("% of Heter. Variants") + xlab(plotname[i]) +
  #   theme(text=element_text(size=9))
  # print(p1)
  # dev.off()
  # 
  # svg(paste0(base_dir,lname[i],"_boxplot.svg"),width=2.166,height=2.166)
  # dfbox$x <- factor(dfbox$x,levels=c("F<10","10<F<100","F>100"))
  # comp <- list(c("F<10","F>100"),c("10<F<100","F>100"))
  # p2 <- ggboxplot(dfbox,x="x",y="y",size=0.5,outlier.size=0.25) +
  #   stat_compare_means(comparisons=comp,size=3) +
  #   ylim(0,max(dfbox$y)*1.25) + ylab("Cochran's Q") +
  #   xlab(plotname[i]) + theme(text=element_text(size=9))
  # print(p2)
  # dev.off()
  
}
print(record)
#write.csv(record,file=filename)