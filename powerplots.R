
library("ggpubr")

base_dir="C:/Users/david/Desktop/"

power_table <- read.csv(paste0(base_dir,"tp_power.csv"),header=T)
pt_ivw <- power_table[1:3,][,c(2:4)]
pt_egger <- power_table[1:3,][,c(2,5:6)]
pplot_ivw <- data.frame(x=unlist(pt_ivw),
                        y=substr(names(unlist(pt_ivw)),2,2),
                        z=paste0(substr(names(unlist(pt_ivw)),3,3),"0"))
pplot_egger <- data.frame(x=unlist(pt_egger),
                          y=substr(names(unlist(pt_egger)),2,2),
                          z=paste0(substr(names(unlist(pt_egger)),3,3),"0"))
tiff(paste0(base_dir,"ivw_power.tif"),width=6.5,height=5,units="in",res=800)
p1<-ggplot(pplot_ivw,aes(factor(z),x,fill=factor(y))) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(name="",labels=c("Observational","MR-IVW, Bonf.",
      "MR-IVW, FDR"),values=c("#1f78b4","#ff7f00","#6a3d9a")) +
  ylab("Statistical Power") + xlab("Measurement Error") + ylim(0,1) +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p1)
dev.off()
tiff(paste0(base_dir,"egger_power.tif"),width=6.5,height=5,units="in",res=800)
p2<-ggplot(pplot_egger,aes(factor(z),x,fill=factor(y))) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(name="",labels=c("Observational","MR-EGGER, Bonf.",
      "MR-EGGER, FDR"),values=c("#1f78b4","#33a02c","#e31a1c")) +
  ylab("Statistical Power") + xlab("Measurement Error") + ylim(0,1) +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p2)
dev.off()

error_table <- read.csv(paste0(base_dir,"fp_power.csv"),header=T)
et_ivw <- error_table[1:3,][,c(2:4)]
et_egger <- error_table[1:3,][,c(2,5:6)]
eplot_ivw <- data.frame(x=unlist(et_ivw),
                        y=paste0(substr(names(unlist(et_ivw)),2,2),"0"),
                        z=substr(names(unlist(et_ivw)),3,3))
eplot_egger <- data.frame(x=unlist(et_egger),
                          y=paste0(substr(names(unlist(et_egger)),2,2),"0"),
                          z=substr(names(unlist(et_egger)),3,3))
tiff(paste0(base_dir,"ivw_error.tif"),width=6.5,height=5,units="in",res=800)
p3<-ggplot(eplot_ivw,aes(factor(z),x,fill=factor(y))) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(name="",labels=c("Observational","MR-IVW, Bonf.",
      "MR-IVW, FDR"),values=c("#1f78b4","#ff7f00","#6a3d9a")) +
  ylab("Type I Error") + xlab("Measurement Error") + ylim(0,1) +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p3)
dev.off()

tiff(paste0(base_dir,"egger_error.tif"),width=6.5,height=5,units="in",res=800)
p4<-ggplot(eplot_egger,aes(factor(z),x,fill=factor(y))) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_fill_manual(name="",labels=c("Observational","MR-EGGER, Bonf.",
      "MR-EGGER, FDR"),values=c("#1f78b4","#33a02c","#e31a1c")) +
  ylab("Type I Error") + xlab("Measurement Error") + ylim(0,1) +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p4)
dev.off()

beta_table <- read.csv(paste0(base_dir,"tp_betas.csv"),header=T)
beta_table <- beta_table[,-1]

bt_ivw <- beta_table[,c(1:3, 6:8, 11:13)]
bt_egger <- beta_table[,c(1,4:5,6,9:10,11,14:15)]

dfbeta_ivw <- data.frame(w=rep(1:9,each=100),
                         x=rep(1:3*10,each=300),
                         y=unlist(bt_ivw),
                         z=rep(rep(1:3,each=100),3))

tiff(paste0(base_dir,"ivw_sim_beta.tif"),width=6.5,height=5,units="in",res=800)
p5<-ggplot(dfbeta_ivw,aes(group=w,x=x,y=y,fill=factor(z))) + 
  geom_boxplot(alpha=0.7) + ylim(-1.5,3.5) +
  scale_fill_manual(name="",labels=c("Observational","MR-IVW, Bonf.", 
      "MR-IVW, FDR"),values=c("#1f78b4","#ff7f00","#6a3d9a")) +
  scale_x_continuous(breaks=c(10,20,30)) +
  xlab("Measurement Error") + ylab("Estimated Beta (100 trials)") +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p5)
dev.off()

dfbeta_egger <- data.frame(w=rep(1:9,each=100),
                           x=rep(1:3*10,each=300),
                           y=unlist(bt_egger),
                           z=rep(rep(1:3,each=100),3))

tiff(paste0(base_dir,"egger_sim_beta.tif"),width=6.5,height=5,units="in",res=800)
p6<-ggplot(dfbeta_egger,aes(group=w,x=x,y=y,fill=factor(z))) + 
  geom_boxplot(alpha=0.7) + ylim(-1.5,3.5) +
  scale_fill_manual(name="",labels=c("Observational","MR-EGGER, Bonf.", 
      "MR-EGGER, FDR"),values=c("#1f78b4","#33a02c","#e31a1c")) +
  scale_x_continuous(breaks=c(10,20,30)) +
  xlab("Measurement Error") + ylab("Estimated Beta (100 trials)") +
  theme_bw(base_size=17) + theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), legend.position="bottom")
print(p6)
dev.off()
