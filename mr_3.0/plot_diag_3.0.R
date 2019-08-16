library("ggpubr")
library("rmeta")

# diag output: Q_byF, pct_het, F_res
het_output = "C:/Users/David/Desktop/genetics/revision/plot/hetplot_"
forest_output = "C:/Users/David/Desktop/genetics/revision/plot/forest.svg"

PRESSO = T
source("diag_3.0.R")

for(i in 1:3) {
   colnum <- 1:3 + 3*(i-1)
   Q_lipid <- Q_byF[colnum]
   het_lipid <- pct_het[colnum]
   
   F_level <- c(
     rep("F<10",length(Q_lipid[[1]])),
     rep("10<F<100",length(Q_lipid[[2]])),
     rep("F>100",length(Q_lipid[[3]]))
   )
   Q_data <- data.frame(x=F_level,y=unlist(Q_lipid))
   Q_data$x <- factor(Q_data$x,levels=c("F<10","10<F<100","F>100"))
#   comp <-list(c("F<10","F>100"),c("10<F<100","F>100"))
   svg(paste0(het_output,i,".svg"),width=2.166,height=2.166)
   p <- ggboxplot(Q_data,x="x",y="y",size=0.5,outlier.size=0.25) +
#     stat_compare_means(comparisons=comp) +
     ylim(0,max(Q_data$y)*1.25) + ylab("Cochran's Q") +
     xlab(c("LDL-C","HDL-C","TG")[i]) + theme(text=element_text(size=9))
   print(p)
   dev.off()
   
   pct_data <- data.frame(x=c("F<10","10<F<100","F>100"),y=het_lipid)
   pct_data$x <- factor(pct_data$x,levels=c("F<10","10<F<100","F>100"))
   svg(paste0(het_output,i+3,".svg"),width=2.166,height=2.166)
   p2 <- ggbarplot(pct_data,x="x",y="y") + ylim(0,max(pct_data$y)*1.25) +
     ylab("% of Heter. Variants") + xlab(c("LDL-C","HDL-C","TG")[i]) +
     theme(text=element_text(size=9))
   print(p2)
   dev.off()
}

# forest plot of F_res
PRESSO = T
PARTIAL = F
source("mvmr_3.0.R")
res_allF <- ivw
F_res <- append(append(append(F_res,list(res_allF),9),list(res_allF),6),list(res_allF),3)

beta_byF <- matrix(nrow=3,ncol=15)
n_variants <- vector()
for(j in 1:12) {
  k <- ifelse(j < 5, 1, ifelse(j < 9, 2, 3))
  beta_byF[1,j+k] <- F_res[[j]]$coef[k]
  beta_byF[2,j+k] <- F_res[[j]]$coef[k] - qt(0.975,F_res[[j]]$df.residual)*summary(F_res[[j]])$coef[k,2]
  beta_byF[3,j+k] <- F_res[[j]]$coef[k] + qt(0.975,F_res[[j]]$df.residual)*summary(F_res[[j]])$coef[k,2]
  n_variants[j] <- F_res[[j]]$df.residual+3
}
OR_byF <- exp(beta_byF)

forest_labels <- cbind(
  c("",rep("LDL-C",4),"",rep("HDL-C",4),"",rep("TG",4),""),
  c("F-statistic",rep(c("F<10","10<F<100","F>100","All",""),3)),
  c("  # Instr.",n_variants[1:4],"",n_variants[5:8],"",n_variants[9:12],""),
  c("OR",round(OR_byF[1,2:5],digits=2),"",round(OR_byF[1,7:10],digits=2),
    "",round(OR_byF[1,12:15],digits=2),"")
)

svg(forest_output,width=5,height=4)
p3 <- forestplot(forest_labels,OR_byF[1,],OR_byF[2,],OR_byF[3,],
           is.summary=c(TRUE,rep(c(FALSE,FALSE,FALSE,TRUE,TRUE),3)),zero=1)
print(p3)
dev.off()