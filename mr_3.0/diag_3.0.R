# store Q, heterogeneity, and regressions for variants binned by F-stat
# Q list: 1-3 LDL, 4-6 HDL, 7-9 TG by increasing F bin
Q_byF <- list()
pct_het <- vector()
F_res <- list()
# assume se is in the column after beta in table
for(k in 1:3) {
  PARTIAL = F
  source("mvmr_3.0.R")
  Q_allF <- Q_snp
  
  c <- covs[[k]]
  bi <- which(colnames(c) %in% c("LDL_beta","HDL_beta","TG_beta"))
  bai <- bi[-k]
  cov_pred <- lm(c[,bi[k]] ~ c[,bai[1]] + c[,bai[2]] +0)
  adj_cov <- c[,bi[k]] - cov_pred$coef[1]*c[,bai[1]] - cov_pred$coef[2]*c[,bai[2]]
  
  F_stat <- adj_cov^2/c[,bi[k]+1]^2
  F_0.10 <- which(F_stat < 10)
  F_10.100 <- which(F_stat >= 10 & F_stat < 100)
  F_100.Inf <- which(F_stat >= 100)
  F_ranges <- list(F_0.10,F_10.100,F_100.Inf)
  
  PARTIAL = T
  for(m in 1:3) {
    Q_F <- Q_allF[F_ranges[[m]]]
    phet <- pchisq(Q_F,df=1,lower.tail=FALSE)
    pct_het[(k-1)*3+m] <- length(which(phet < 0.05))/length(phet) * 100
    Q_byF[[(k-1)*3+m]] <- Q_F
    
    part <- F_ranges[[m]]
    source("mvmr_3.0.R")
    F_res[[(k-1)*3+m]] <- ivw
  }
}
