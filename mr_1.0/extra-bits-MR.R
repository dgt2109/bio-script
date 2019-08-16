


CARDIO_MRMVInputObject = mr_mvinput(bx=cbind(CARDIO_data$LDL.beta, CARDIO_data$HDL.beta, CARDIO_data$TG.beta),
                                    bxse=cbind(CARDIO_data$LDL.se, CARDIO_data$HDL.se, CARDIO_data$TG.se),
                                    by=CARDIO_data$CARDIO.beta,
                                    byse=CARDIO_data$CARDIO.se)
CARDIO_MRMVObject = mr_mvivw(CARDIO_MRMVInputObject,
                             model="default",
                             correl="FALSE",
                             distribution="normal",
                             alpha=0.05)
MRMVObject = mr_mvivw(CARDIO_MRMVInputObject)


UKBB_MRMVInputObject = mr_mvinput(bx=cbind(UKBB_data$LDL.beta, UKBB_data$HDL.beta, UKBB_data$TG.beta),
                                  bxse=cbind(UKBB_data$LDL.se, UKBB_data$HDL.se, UKBB_data$TG.se),
                                  by=UKBB_data$UKBB.beta,
                                  byse=UKBB_data$UKBB.se)
UKBB_MRMVObject = mr_mvivw(UKBB_MRMVInputObject,
                           model="default",
                           correl="FALSE",
                           distribution="normal",
                           alpha=0.05)
MRMVObject = mr_mvivw(UKBB_MRMVInputObject)


capture.output(print(CARDIO_MRMVObject),
               file="C:/Users/David/Desktop/humgen/output/CARDIO_mvivw.txt")
capture.output(print(UKBB_MRMVObject),
               file="C:/Users/David/Desktop/humgen/output/UKBB_mvivw.txt")

#library("MRPRESSO")
#SummaryStats=all_data[ , c("rsid","CARDIO.beta","CARDIO.se","HDL.beta","HDL.se")]
#mr_presso(BetaOutcome="CARDIO.beta", BetaExposure="HDL.beta", 
#          SdOutcome="CARDIO.se", SdExposure="HDL.se",
#          OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=SummaryStats,
#          NbDistribution=100000, SignifThreshold=0.05)

# cooksd by VAF
#cooks_dplot=paste(base_dir,table_id,"_cooksd.png",sep="")
#png(cooks_dplot,width=6.5,height=2.5,units="in",res=600)
#par(mfrow=c(2,2),mai=rep(0.2,4),oma=c(3,3,0,1))
#data=data[order(data$VAF), ]
#cooksd=cooks.distance(naive.model)
#plot(data$VAF,cooksd,pch=20,main="Cook's D of OLS")
#for (i in 1:length(data$VAF)) {
#  ly = seq(0,as.integer(cooksd[i]*10000))/10000
#  lines(rep(data$VAF[i],length(ly)),ly,lty=3)
#}
#lines(rep(0.05,601),seq(0,600)/1000,col=34)
#lines(seq(0,50)/100,rep(4/length(data$VAF),51),lty=2)
#title(xlab="VAF",ylab="Cook's D",outer=TRUE,line=1,cex.lab=1.2)
#dev.off()
