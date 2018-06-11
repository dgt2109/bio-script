
library(simex)

# files
base_dir="C:/Users/David/Desktop/genetics/data/simple/"
table_id="HDL_simple"
plot_title=strsplit(table_id,"_")[[1]][1]
if(plot_title != "TG") {
  plot_title=paste(plot_title,"-C",sep="")
}
xlabel=paste(plot_title," beta",sep="")
data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
plot_file=paste(base_dir,table_id,".png",sep="")
ols_summary=paste(base_dir,table_id,"_ols.txt",sep="")
simex_summary=paste(base_dir,table_id,"_simex.txt",sep="")

data=data[data$Common=="TRUE", ]

# regressions
w=data$lipid_beta
y=data$CAD_beta

naive.model=lm(y ~ w, x=TRUE)

error=cbind(data$lipid_se,data$CAD_se)
simex.model=simex(model=naive.model,
                  SIMEXvariable=c("w","y"),
                  measurement.error=error,
                  lambda=c(0.5,1.0,1.5,2.0), 
                  B=100,
                  fitting.method="quad", 
                  asymptotic=F)

# plotting
png(plot_file,width=4,height=4.5,units="in",res=600)
#par(mfrow=c(1,1),mai=rep(0.2,4),oma=c(3,3,0,1))
plot(data$lipid_beta,data$CAD_beta,main=plot_title,xlab=xlabel,ylab="CAD beta",
     pch=20,cex.axis=1.0,cex.main=1.0)
abline(naive.model$coefficients[1],naive.model$coefficients[2])
test_points=data.frame(w=seq(-100,100)/200)
prd_CI = predict(naive.model,newdata=test_points,interval=c("confidence"),level=0.95,type="response")
lines(test_points$w,prd_CI[,2],lty=3)
lines(test_points$w,prd_CI[,3],lty=3)
dev.off()

capture.output(summary(naive.model),file=ols_summary)
capture.output(summary(simex.model),file=simex_summary)

