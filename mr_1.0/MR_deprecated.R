library(MendelianRandomization)
###########################################################################
# data import
###########################################################################

LDL=read.delim("C:/Users/David/Desktop/humgen/jointGwasMc_LDL_185.txt",header=FALSE)
HDL=read.delim("C:/Users/David/Desktop/humgen/jointGwasMc_HDL_185.txt",header=FALSE)
TG=read.delim("C:/Users/David/Desktop/humgen/jointGwasMc_TG_185.txt",header=FALSE)
CARDIO=read.delim("C:/Users/David/Desktop/humgen/CARDIoGRAM_GWAS_RESULTS_185.txt",header=FALSE)
UKBB=read.delim("C:/Users/David/Desktop/humgen/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517_185.txt",header=FALSE)
SUPPL=read.delim("C:/Users/David/Desktop/humgen/supplement_data.txt",header=TRUE)

colnames(LDL) <- c("loc","loc2","rsid","a1","a2","beta","se","N","p","freq")
LDL[LDL$freq > 0.5, ]$beta <- -LDL[LDL$freq > 0.5, ]$beta
LDL[LDL$freq > 0.5, ]$freq <- 1-LDL[LDL$freq > 0.5, ]$freq

colnames(HDL) <- c("loc","loc2","rsid","a1","a2","beta","se","N","p","freq")
HDL[HDL$freq > 0.5, ]$beta <- -HDL[HDL$freq > 0.5, ]$beta
HDL[HDL$freq > 0.5, ]$freq <- 1-HDL[HDL$freq > 0.5, ]$freq  

colnames(TG) <- c("loc","loc2","rsid","a1","a2","beta","se","N","p","freq")
TG[TG$freq > 0.5, ]$beta <- -TG[TG$freq > 0.5, ]$beta
TG[TG$freq > 0.5, ]$freq <- 1-TG[TG$freq > 0.5, ]$freq

colnames(CARDIO)=c("rsid","loc","a1","a2", "freq", "p", "hetp", 
                    "CARDIO.beta", "CARDIO.se", "N_case", "N_control", "model")
CARDIO[CARDIO$freq > 0.5, ]$CARDIO.beta = -CARDIO[CARDIO$freq > 0.5, ]$CARDIO.beta
CARDIO[CARDIO$freq > 0.5, ]$freq = 1-CARDIO[CARDIO$freq > 0.5, ]$freq

colnames(UKBB)=c("name","rsid","chr","loc","a1","a2", "freq", "UKBB.beta", 
                 "UKBB.se", "p","N", "exome", "info")
UKBB[UKBB$freq > 0.5, ]$UKBB.beta = -UKBB[UKBB$freq > 0.5, ]$UKBB.beta
UKBB[UKBB$freq > 0.5, ]$freq = 1-UKBB[UKBB$freq > 0.5, ]$freq

# get se for supplemental data from t-distribution
SUPPL_LDL.t=qt(SUPPL$LDL_p,df=182)
SUPPL_HDL.t=qt(SUPPL$HDL_p,df=182)
SUPPL_TG.t=qt(SUPPL$TG_p,df=182)
SUPPL_CAD.t=qt(SUPPL$CAD_p,df=182)
SUPPL_LDL.se=abs(SUPPL$LDL_beta/SUPPL_LDL.t)
SUPPL_HDL.se=abs(SUPPL$HDL_beta/SUPPL_HDL.t)
SUPPL_TG.se=abs(SUPPL$TG_beta/SUPPL_TG.t)
SUPPL_CAD.se=abs(SUPPL$CAD_beta/SUPPL_CAD.t)
SUPPL$LDL.se=SUPPL_LDL.se
SUPPL$HDL.se=SUPPL_HDL.se
SUPPL$TG.se=SUPPL_TG.se
SUPPL$CAD.se=SUPPL_CAD.se

lipid_data = data.frame(LDL$rsid,LDL$beta,LDL$se,
                   HDL$beta,HDL$se,
                   TG$beta, TG$se)

UKBB_data = merge(x = lipid_data, y = UKBB[, c("rsid","UKBB.beta","UKBB.se")], 
                 by.x="LDL.rsid", by.y="rsid")
colnames(UKBB_data)[1] = "rsid"
CARDIO_data = merge(x = lipid_data, y = CARDIO[, c("rsid","CARDIO.beta","CARDIO.se")], 
                 by.x="LDL.rsid", by.y="rsid")
colnames(CARDIO_data)[1] = "rsid"

#testing effect of rs5880
UKBB_data = UKBB_data[!UKBB_data$rsid=="rs5880", ]

###########################################################################
# compute regression using 3 methods
###########################################################################

# Original data from supplement
SUPPL_ols=lm(SUPPL$CAD_beta ~ SUPPL$LDL_beta + SUPPL$HDL_beta + SUPPL$TG_beta +0)

SUPPL_err=1/(SUPPL$CAD.se^2)

#SUPPL_wls=lm(SUPPL$CAD_beta ~ SUPPL$LDL_beta + SUPPL$HDL_beta + SUPPL$TG_beta +0,
#             weights = SUPPL_err)


# CARDIoGRAM alone
CARDIO_ols=lm(CARDIO_data$CARDIO.beta ~ CARDIO_data$LDL.beta + CARDIO_data$HDL.beta + 
              CARDIO_data$TG.beta +0)
CARDIO_err=1/(CARDIO_data$CARDIO.se^2)# + 1/(CARDIO_data$LDL.se^2) +
          #      1/(CARDIO_data$HDL.se^2) + 1/(CARDIO_data$TG.se^2)
CARDIO_wls=lm(CARDIO_data$CARDIO.beta ~ CARDIO_data$LDL.beta + CARDIO_data$HDL.beta + 
              CARDIO_data$TG.beta +0,weights = CARDIO_err)

# CARDIoGRAMC4D_UKBB
UKBB_ols=lm(UKBB_data$UKBB.beta ~ UKBB_data$LDL.beta + UKBB_data$HDL.beta + 
              UKBB_data$TG.beta +0)
UKBB_err=1/(UKBB_data$UKBB.se^2)# + 1/(UKBB_data$LDL.se^2) +
#  1/(UKBB_data$HDL.se^2) + 1/(UKBB_data$TG.se^2)
UKBB_wls=lm(UKBB_data$UKBB.beta ~ UKBB_data$LDL.beta + UKBB_data$HDL.beta + 
              UKBB_data$TG.beta +0,weights = UKBB_err)

###########################################################################
# output 3 datasets
###########################################################################

# Original data from supplement
capture.output(summary(SUPPL_ols),
               file="C:/Users/David/Desktop/humgen/output/SUPPL_ols.txt")
capture.output(summary(SUPPL_wls),
               file="C:/Users/David/Desktop/humgen/output/SUPPL_wls.txt")
imagefile="C:/Users/David/Desktop/humgen/output/SUPPL.png"
png(imagefile,width=9,height=4,units="in",res=600)
par(mfrow=c(1,3))
plot(SUPPL$LDL_beta, SUPPL$CAD_beta,main="LDL",
     xlab="LDL effect", ylab="CAD effect",
     xlim=c(-0.5,0.2),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.5,0.2,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,SUPPL_ols$coefficients[1])
plot(SUPPL$HDL_beta, SUPPL$CAD_beta,main="HDL",
     xlab="HDL effect", ylab="CAD effect",
     xlim=c(-0.4,0.5),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.4,0.5,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,SUPPL_ols$coefficients[2])
plot(SUPPL$TG_beta, SUPPL$CAD_beta,main="TG",
     xlab="TG effect", ylab="CAD effect",
     xlim=c(-0.2,0.3),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.2,0.3,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,SUPPL_ols$coefficients[3])
dev.off()

imagefile="C:/Users/David/Desktop/humgen/output/SUPPL_funnel.png"
png(imagefile,width=5,height=5,units="in",res=600)
#imagefile="C:/Users/David/Desktop/humgen/output/SUPPL_volcano.png"
plot(SUPPL$HDL_beta, SUPPL$HDL.se)
dev.off()
#png(imagefile,width=5,height=5,units="in",res=600)
#plot(SUPPL$HDL_beta, -log(SUPPL$HDL_p),main="HDL SNPs volcano",
#     xlab="HDL effect", ylab="-log p",
#     xlim=c(-0.5,0.5),
#     axes=F)
#axis(1,pos=0,seq(-0.5,0.5,0.1))
#axis(2,pos=0)
#dev.off()

# CARDIoGRAM alone
capture.output(summary(CARDIO_ols),
               file="C:/Users/David/Desktop/humgen/output/CARDIO_ols.txt")
capture.output(summary(CARDIO_wls),
               file="C:/Users/David/Desktop/humgen/output/CARDIO_wls.txt")

imagefile="C:/Users/David/Desktop/humgen/output/CARDIOGRAM_2011.png"
png(imagefile,width=9,height=4,units="in",res=600)
par(mfrow=c(1,3))
plot(CARDIO_data$LDL.beta, CARDIO_data$CARDIO.beta,main="LDL",
     xlab="LDL effect", ylab="CAD effect",
     xlim=c(-0.5,0.2),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.5,0.2,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,CARDIO_MRMVObject$Estimate[1])
plot(CARDIO_data$HDL.beta, CARDIO_data$CARDIO.beta,main="HDL",
     xlab="HDL effect", ylab="CAD effect",
     xlim=c(-0.4,0.2),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.4,0.2,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,CARDIO_MRMVObject$Estimate[2])
plot(CARDIO_data$TG.beta, CARDIO_data$CARDIO.beta,main="TG",
     xlab="TG effect", ylab="CAD effect",
     xlim=c(-0.2,0.3),ylim=c(-0.2,0.2),
     axes=F)
axis(1,pos=0,seq(-0.2,0.3,0.1))
axis(2,pos=0,seq(-0.2,0.2,0.1))
abline(0,CARDIO_MRMVObject$Estimate[3])
dev.off()

# CARDIoGRAMC4D_UKBB
capture.output(summary(UKBB_ols),
               file="C:/Users/David/Desktop/humgen/output/UKBB_ols.txt")
capture.output(summary(UKBB_wls),
               file="C:/Users/David/Desktop/humgen/output/UKBB_wls.txt")

imagefile="C:/Users/David/Desktop/humgen/output/CARDIOGRAMC4D_UKBB_2017.png"
png(imagefile,width=9,height=4,units="in",res=600)
par(mfrow=c(1,3))
plot(UKBB_data$LDL.beta, UKBB_data$UKBB.beta,main="LDL",
     xlab="LDL effect", ylab="CAD effect",
     xlim=c(-0.5,0.2),ylim=c(-0.2,0.1),
     axes=F)
axis(1,pos=0,seq(-0.5,0.2,0.1))
axis(2,pos=0,seq(-0.2,0.1,0.1))
abline(0,UKBB_MRMVObject$Estimate[1])
plot(UKBB_data$HDL.beta, UKBB_data$UKBB.beta,main="HDL",
     xlab="HDL effect", ylab="CAD effect",
     xlim=c(-0.4,0.2),ylim=c(-0.2,0.1),
     axes=F)
axis(1,pos=0,seq(-0.4,0.2,0.1))
axis(2,pos=0,seq(-0.2,0.1,0.1))
abline(0,UKBB_MRMVObject$Estimate[2])
plot(UKBB_data$TG.beta, UKBB_data$UKBB.beta,main="TG",
     xlab="TG effect", ylab="CAD effect",
     xlim=c(-0.2,0.3),ylim=c(-0.2,0.1),
     axes=F)
axis(1,pos=0,seq(-0.2,0.3,0.1))
axis(2,pos=0,seq(-0.2,0.1,0.1))
abline(0,UKBB_MRMVObject$Estimate[3])
dev.off()
