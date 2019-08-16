library(simex)

# files
base_dir="C:/Users/David/Desktop/genetics/data/"
table_id="185/table_179"
data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
plot_file=paste(base_dir,table_id,".png",sep="")
ols_summary=paste(base_dir,table_id,"_ols.txt",sep="")
simex_summary=paste(base_dir,table_id,"_simex.txt",sep="")

data=data[data$Common=="TRUE", ]

# regressions
W=cbind(data$LDL_beta,data$HDL_beta,data$TG_beta+0)
y=data$CAD_beta
mydata=as.data.frame(cbind(W,y))
naive.model=lm(y ~ V1+V2+V3, data=mydata, x=TRUE)
error=cbind(data$LDL_se,data$HDL_se,data$TG_se,data$CAD_se)
simex.model=simex(model=naive.model,
                  SIMEXvariable=c("V1","V2","V3","y"),
                  measurement.error=error,
                  lambda=c(0.5,1.0,1.5,2.0), 
                  B=10000,
                  fitting.method="quad", 
                  asymptotic=F)

# plotting
png(plot_file,width=13,height=4,units="in",res=600)
par(mfrow=c(1,3),mai=rep(0.2,4),oma=c(3,3,0,1))
plot(data$LDL_beta,data$CAD_beta,main="LDL-C",pch=20,cex.axis=1.6,cex.main=1.6)
abline(naive.model$coefficients[1],naive.model$coefficients[2])
test_LDL=data.frame(cbind(seq(-120,120)/200,rep(0,241),rep(0,241)))
colnames(test_LDL)=c("V1","V2","V3")
prd_LDL = predict(naive.model,newdata=test_LDL,interval=c("confidence"),
                    level=0.95,type="response")
lines(test_LDL$V1,prd_LDL[,2],lty=3)
lines(test_LDL$V1,prd_LDL[,3],lty=3)
  
plot(data$HDL_beta,data$CAD_beta,main="HDL-C",pch=20,cex.axis=1.6,cex.main=1.6)
abline(naive.model$coefficients[1],naive.model$coefficients[3])
test_HDL=data.frame(cbind(rep(0,201),seq(-100,100)/200,rep(0,201)))
colnames(test_HDL)=c("V1","V2","V3")
prd_HDL = predict(naive.model,newdata=test_HDL,interval=c("confidence"),
                    level=0.95,type="response")
lines(test_HDL$V2,prd_HDL[,2],lty=3)
lines(test_HDL$V2,prd_HDL[,3],lty=3)
  
plot(data$TG_beta,data$CAD_beta,main="TG",pch=20,cex.axis=1.6,cex.main=1.6)
abline(naive.model$coefficients[1],naive.model$coefficients[4])
test_TG=data.frame(cbind(rep(0,201),rep(0,201),seq(-100,100)/200))
colnames(test_TG)=c("V1","V2","V3")
prd_TG = predict(naive.model,newdata=test_TG,interval=c("confidence"),
                  level=0.95,type="response")
lines(test_TG$V3,prd_TG[,2],lty=3)
lines(test_TG$V3,prd_TG[,3],lty=3)

title(xlab="Lipid beta",ylab="CAD beta",outer=TRUE,line=1,cex.lab=1.6)
  
dev.off()

# print to text file
capture.output(summary(naive.model),file=ols_summary)
capture.output(summary(simex.model),file=simex_summary)


