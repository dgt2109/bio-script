library("ggpubr")

# byvariant output: list of matrices of results
# col order: ivw beta, ivw se, I2 ivw, egger beta, egger se, I2 egger, egger int
source("byvariant_3.0.R")

byvar_output = "C:/Users/david/Desktop/genetics/revision/plot/byvar_"
x <- 5:615
for (i in 2:3) {
  if (i == 2) { bymin=-0.15; bymax=0.05; symax=0.11; pymin=-4 }
  if (i == 3) { bymin=0; bymax=0.3; symax=0.2; pymin=-6 }
  # plot beta trend
  betas_data <- data.frame(cbind(x),ivw_betas=byvar_results[[i]][,1],
                           egger_betas=byvar_results[[i]][,4])
  svg(paste0(byvar_output,(i-2)*4+1,".svg"),5,2.5)
  p1 <- ggplot(betas_data,aes(x=x,y=ivw_betas,colour="MVMR-IVW"))+geom_point(size=0.5) +
    geom_point(aes(x=x,y=egger_betas,colour="MVMR-EGGER"), size=0.5) +
    scale_colour_manual(name="",values=c("MVMR-IVW"="black","MVMR-EGGER"="red"),
                        guide=guide_legend(reverse=T)) +
    scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
    ylim(bymin,bymax) + xlab("Number of Variants") +
    guides(colour = guide_legend(reverse=T,override.aes=list(size=1.5))) +
    theme_bw() + theme(panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.margin=margin(5,10,5,10),
                       legend.text=element_text(size=11),
                       legend.justification=c(0,0),
                       legend.position=c(0.6,0.55))
  if (i == 2) { p1 <- p1 + ylab(expression(paste("HDL-C ",beta," Coefficient"))) }
  if (i == 3) { p1 <- p1 + ylab(expression(paste("TG ",beta," Coefficient"))) }
  print(p1)
  dev.off()
  
  # plot se trend
  ses_data <- data.frame(cbind(x),ivw_ses=byvar_results[[i]][,2],
                         egger_ses=byvar_results[[i]][,5])
  svg(paste0(byvar_output,(i-2)*4+2,".svg"),5,2.5)
  p2 <- ggplot(ses_data,aes(x=x,y=ivw_ses,colour="MVMR-IVW"))+geom_point(size=0.5) +
    geom_point(aes(x=x,y=egger_ses,colour="MVMR-EGGER"),size=0.5) +
    scale_colour_manual(name="",values=c("MVMR-IVW"="black","MVMR-EGGER"="red"),
                        guide=guide_legend(reverse=T)) +
    scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
    ylim(0,symax) + xlab("Number of Variants") +
    guides(colour = guide_legend(reverse=T,override.aes=list(size=1.5,
                                                             linetype=rep("blank",2)))) +
    geom_line(aes(x=x,y=abs(tail(betas_data$egger_betas,1)/1.96)),linetype=2) +
    theme_bw() + theme(panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.margin=margin(5,10,5,10),
                       legend.text=element_text(size=11),
                       legend.justification=c(0,0),
                       legend.position=c(0.6,0.55))
  if (i == 2) { p2 <- p2 + ylab(expression(paste("HDL-C ",beta," Coefficient SE"))) }
  if (i == 3) { p2 <- p2 + ylab(expression(paste("TG ",beta," Coefficient SE"))) }
  print(p2)
  dev.off()
  
  # plot I2 trend
  I2s_data <- data.frame(cbind(x),ivw_hets=byvar_results[[i]][,3],
                         egger_hets=byvar_results[[i]][,6])
  svg(paste0(byvar_output,(i-2)*4+3,".svg"),5,2.5)
  p3 <- ggplot(I2s_data,aes(x=x,y=egger_hets,colour="MVMR-EGGER"))+geom_point(size=0.1) +
    geom_point(aes(x=x,y=ivw_hets,colour="MVMR-IVW"),size=0.1) +
    scale_colour_manual(name="",values=c("MVMR-IVW"="black","MVMR-EGGER"="red")) +
    scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
    ylim(0,100) + xlab("Number of Variants") +
    ylab(expression(paste("Heterogeneity Statistic ",I[H]^2," (%)"))) +
    #  guides(colour = guide_legend(override.aes=list(size=1.5))) +
    guides(colour = guide_legend(reverse=T,override.aes=list(size=1.5))) +
    theme_bw() + theme(panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.margin=margin(5,10,5,10),
                       legend.text=element_text(size=11),
                       legend.justification=c(0,0),
                       legend.position=c(0.6,0.05))
  print(p3)
  dev.off()
    
  # plot egger int trend
  int_data <- data.frame(cbind(x),egger_inters=byvar_results[[i]][,7])
  svg(paste0(byvar_output,(i-2)*4+4,".svg"),5,2.5)
  p4 <- ggplot(int_data,aes(x=x,y=egger_inters)) +
    geom_point(color="red", size=0.5) +
    scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
    xlab("Number of Variants") +
    ylab(paste0("MR-EGGER Intercept,",c("LDL-C","HDL-C","TG")[i])) +
    geom_line(aes(x=x,y=0)) +
    theme_bw() + theme(panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.margin=margin(5,10,5,10))
  print(p4)
  dev.off()
}
