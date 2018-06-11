

#files
base_dir="C:/Users/David/Desktop/genetics/data/update/"
table_id="update_table"
data=read.delim(paste(base_dir,table_id,".txt",sep=""),header=TRUE)
plot_file=paste(base_dir,table_id,"-2.png",sep="")


png(plot_file,width=13,height=7,units="in",res=600)
par(mfrow=c(1,2))#,mai=rep(0.6,4),oma=c(3,3,0,1))

plot(data$CAD2011_beta,data$CAD2017_beta,main="CAD beta by variant in each dataset",
     pch=20,
     xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),
     xlab="CAD 2011 beta", ylab="CAD 2017 beta")
res=lm(data$CAD2017_beta ~ data$CAD2011_beta)
abline(res$coefficients[1],res$coefficients[2])
abline(0,1,lty=2)

data=data[order(-data$CAD2011_se), ]
barplot(data$CAD2017_se-data$CAD2011_se,
        main="Change in CAD SE by variant",
        xlab="Variant (ranked by CAD 2011 SE)",
        ylab="Change in CAD SE from 2011 to 2017 CAD data")

dev.off()