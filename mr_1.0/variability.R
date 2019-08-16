
#files
base_dir="C:/Users/David/Desktop/genetics/data/"
table_id="185/table_179"
data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
plot_file=paste(base_dir,table_id,"_variability-2.png",sep="")

for (i in 1:length(data$VAF)) {
  if (data$VAF[i]>0.5) {
    data$VAF[i]=1-data$VAF[i]
  }
}
invrootN = 1/sqrt(data$VAF)

#plotting
png(plot_file,width=13,height=5,units="in",res=600)
par(mfrow=c(2,2),mai=rep(0.2,4),oma=c(3,3,0,1))

plot(data$VAF,data$LDL_se,pch=20,main="LDL-C SE")
for (i in 1:length(data$VAF)) {
  ly = seq(0,as.integer(data$LDL_se[i]*10000))/10000
  lines(rep(data$VAF[i],length(ly)),ly,lty=3)
}
lines(rep(0.05,101),seq(0,100)/1000,col=34)

res1=lm(data$LDL_se ~ invrootN)
a = seq(0,500)/1000
b = res1$coefficients[1] + res1$coefficients[2] * 1/sqrt(a)
lines(a,b,col=28,lty=5)
  
plot(data$VAF,data$HDL_se,pch=20,main="HDL-C SE")
for (i in 1:length(data$VAF)) {
  ly = seq(0,as.integer(data$HDL_se[i]*10000))/10000
  lines(rep(data$VAF[i],length(ly)),ly,lty=3)
}
lines(rep(0.05,101),seq(0,100)/1000,col=34)

res2=lm(data$HDL_se ~ invrootN)
a = seq(0,500)/1000
b = res2$coefficients[1] + res2$coefficients[2] * 1/sqrt(a)
lines(a,b,col=28,lty=5)

plot(data$VAF,data$TG_se,pch=20, main="TG SE")
for (i in 1:length(data$VAF)) {
  ly = seq(0,as.integer(data$TG_se[i]*10000))/10000
  lines(rep(data$VAF[i],length(ly)),ly,lty=3)
}
lines(rep(0.05,101),seq(0,100)/1000,col=34)

res3=lm(data$TG_se ~ invrootN)
a = seq(0,500)/1000
b = res3$coefficients[1] + res3$coefficients[2] * 1/sqrt(a)
lines(a,b,col=28,lty=5)

plot(data$VAF,data$CAD_se,pch=20, main="CAD SE")
for (i in 1:length(data$VAF)) {
  ly = seq(0,as.integer(data$CAD_se[i]*10000))/10000
  lines(rep(data$VAF[i],length(ly)),ly,lty=3)
}
lines(rep(0.05,101),seq(0,100)/1000,col=34)

res4=lm(data$CAD_se ~ invrootN)
a = seq(0,500)/1000
b = res4$coefficients[1] + res4$coefficients[2] * 1/sqrt(a)
lines(a,b,col=28,lty=5)

title(main="Lipid beta",xlab="VAF",ylab="Std. Error",outer=TRUE,
      line=1,cex.lab=1.6)
dev.off()



