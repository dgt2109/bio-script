dof = 4
nmin = 1+dof
nmax = 631

ivw_all <- vector(mode="list",length=nmax-dof)
egger_all <- vector(mode="list",length=nmax-dof)

rnames <- c("ivw_hets","egger_hets","ivw_betas","ivw_ses",
            "egger_betas","egger_ses","ivw_ps","egger_ps",
            "egger_inters","egger_inter_ps")
for (i in 1:length(rnames)) {
  eval(parse(text=paste0(rnames[i],
             "<- vector(mode=\"double\",length=nmax-dof)")))
}

for (k in nmin:nmax) {
  nvars <- k
  source("MV_MR-egger.R")
  ivw_all[[k-dof]] <- ivw_res
  egger_all[[k-dof]] <- egger_res
  ivw_hets[k-dof] <- i2h_ivw
  egger_hets[k-dof] <- i2h_egger
  if(k%%10==0) { print(k) }
}

for (m in 1:(nmax-dof)) {
  ivw <- ivw_all[[m]]
  egger <- egger_all[[m]]
  ivw_betas[m] <- ivw$coef[1]
  egger_betas[m] <- egger$coef[2]

  ivw_ses[m]  <- ivw$coef[1,2]
  egger_ses[m] <- egger$coef[2,2]
  
  ivw_ps[m] <- ivw$coef[1,4]
  egger_ps[m] <- egger$coef[2,4]
  
  egger_inters[m] <- egger$coef[1]
  egger_inter_ps[m] <- egger$coef[1,4]
}

library("ggpubr")
base_dir="C:/Users/david/Desktop/genetics/revision/snp680/byvariant/TG/"
#bymin=-0.15; bymax=0.05; symax=0.11; pymin=-4 #params HDL
bymin=0; bymax=0.3; symax=0.2; pymin=-6 # params TG

#make 6 plots
N_variants = 1:(nmax-dof)
df_beta = data.frame(cbind(N_variants,ivw_betas,egger_betas))
svg(paste0(base_dir,"beta_byvariant.svg"),width=5,height=2.5)
ggplot(df_beta,aes(x=N_variants,y=ivw_betas,colour="MVMR-IVW")) +
  geom_point(size=0.5) +
  geom_point(aes(x=N_variants,y=egger_betas,colour="MVMR-EGGER"),
             size=0.5) +
  scale_colour_manual(name="",values=c("MVMR-IVW"="black",
                                       "MVMR-EGGER"="red"),
                      guide=guide_legend(reverse=T)) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
  ylim(bymin,bymax) + xlab("Number of Variants") +
  ylab(expression(paste("TG ",beta," Coefficient"))) +
  guides(colour = guide_legend(reverse=T,override.aes=list(size=1.5))) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.margin=margin(5,10,5,10),
                     legend.text=element_text(size=11),
                     legend.justification=c(0,0),
                     legend.position=c(0.6,0.55))
dev.off()

df_se = data.frame(cbind(N_variants,ivw_ses,egger_ses))
svg(paste0(base_dir,"se_byvariant.svg"),width=5,height=2.5)
ggplot(df_se,aes(x=N_variants,y=ivw_ses,colour="MVMR-IVW")) +
  geom_point(size=0.5) +
  geom_point(aes(x=N_variants,y=egger_ses,colour="MVMR-EGGER"),
             size=0.5) +
  scale_colour_manual(name="",values=c("MVMR-IVW"="black",
                                       "MVMR-EGGER"="red"),
                      guide=guide_legend(reverse=T)) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
  ylim(0,symax) + xlab("Number of Variants") +
  ylab(expression(paste("TG ",beta," Coefficient SE"))) +
  guides(colour = guide_legend(reverse=T,override.aes=list(size=1.5,
                               linetype=rep("blank",2)))) +
  geom_line(aes(x=N_variants,y=abs(tail(egger_betas,1)/1.96))
            ,linetype=2) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.margin=margin(5,10,5,10),
                     legend.text=element_text(size=11),
                     legend.justification=c(0,0),
                     legend.position=c(0.6,0.55))
dev.off()

df_het = data.frame(cbind(N_variants,ivw_hets, egger_hets))
svg(paste0(base_dir,"het_byvariant.svg"),width=5,height=2.5)
ggplot(df_het,aes(x=N_variants,y=egger_hets,colour="MVMR-EGGER")) +
  geom_point(size=0.25) +
  geom_point(aes(x=N_variants,y=ivw_hets,colour="MVMR-IVW"),
             size=0.25) +
  scale_colour_manual(name="",values=c("MVMR-IVW"="black",
                                       "MVMR-EGGER"="red")) +
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
dev.off()
# 
# df_ps = data.frame(cbind(N_variants,ivw_ps, egger_ps))
# svg(paste0(base_dir,"p_byvariant.svg"),width=5,height=2.5)
# ggplot(df_ps,aes(x=N_variants,y=log10(ivw_ps),colour="MR-IVW")) +
#   geom_point(size=0.5) +
#   geom_point(aes(x=N_variants,y=log10(egger_ps),colour="MR-EGGER"),
#              size=0.5) +
#   scale_colour_manual(name="",values=c("MVMR-IVW"="black",
#                                        "MVMR-EGGER"="red")) +
#   scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
#   ylim(pymin,0) + xlab("Number of Variants") +
#   ylab(expression(paste("log TG ",beta," Coefficient P"))) +
#   guides(colour = guide_legend(override.aes=list(size=2))) +
#   geom_line(aes(x=N_variants,y=log10(0.05)),linetype=2) +
#   theme_bw() + theme(panel.grid.major=element_blank(),
#                      panel.grid.minor=element_blank(),
#                      plot.margin=margin(5,10,5,10),
#                      legend.text=element_text(size=11),
#                      legend.justification=c(0,0),
#                      legend.position=c(0.02,0.02))
# dev.off()
# 
df_inter = data.frame(cbind(N_variants,egger_inters, egger_inter_ps))
svg(paste0(base_dir,"inter_byvariant.svg"),width=5,height=2.5)
ggplot(df_inter,aes(x=N_variants,y=egger_inters)) +
  geom_point(color="red", size=0.5) +
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
  xlab("Number of Variants") +
  ylab(paste0("MR-EGGER Intercept,TG")) +
  geom_line(aes(x=N_variants,y=0)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.margin=margin(5,10,5,10))
dev.off()
# 
# svg(paste0(base_dir,"inter_p_byvariant.svg"),width=5,height=2.5)
# ggplot(df_inter,aes(x=N_variants,y=log10(egger_inter_ps))) +
#   geom_point(color="red", size=0.5) +
#   #  scale_colour_manual(name="",values=c("MVMR-IVW"="black","MVMR-EGGER"="red")) +
#   scale_x_continuous(expand=c(0.01,0.01),breaks=c(1:6*100)) +
#   xlab("Number of Variants") +
#   ylab("MR-EGGER Inter. P, TG") +
#   geom_line(aes(x=N_variants,y=log10(0.05)),linetype=2) +
#   theme_bw() + theme(panel.grid.major=element_blank(),
#                      panel.grid.minor=element_blank(),
#                      plot.margin=margin(5,10,5,10),
#                      legend.text=element_text(size=11),
#                      legend.justification=c(0,0),
#                      legend.position=c(0.65,0.55))
# dev.off()
