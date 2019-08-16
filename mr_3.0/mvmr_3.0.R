input="C:/Users/David/Desktop/genetics/revision/table_620.txt"

data <- read.delim(input,header=TRUE,stringsAsFactors=FALSE)
data <- data[order(pmin(data$LDL_q,data$HDL_q,data$TG_q)),]

lname <- c("LDL_beta","HDL_beta","TG_beta")
# orientation
covs <- list()
ocovs <- list()
errors <- list()
for(i in 1:3) {
  cov <- data
  beta_index <- which(colnames(cov) %in% c("LDL_beta","HDL_beta","TG_beta","CAD_beta"))
  lb <- beta_index[i]
  beta_index <- beta_index[-i]
  for(j in beta_index) { cov[,j] <- ifelse(cov[,lb]>0,cov[,j],cov[,j]*-1) }
  cov[,lb] <- abs(cov[,lb])
  
  if(PARTIAL == T) { cov <- cov[part,] }
  oind <- which(cov$rsid %in% c("rs2315065","rs934287","rs1250229","rs653178","rs12940887"))
  if(PRESSO == T) { 
    ocovs[[i]] <- cov[oind,]; 
    if(length(oind) > 0) { cov <- cov[-oind,] } 
  }
  covs[[i]] <- cov
}

# regression & heterogeneity analysis
eggers <- list()
Q_snps <- list()
I = matrix(nrow=2,ncol=3)
for(i in 1:3) {
  W <- covs[[i]][,c("LDL_beta","HDL_beta","TG_beta","CAD_beta","CAD_se")]
  ivw <- lm(CAD_beta ~ LDL_beta+HDL_beta+TG_beta +0, data=W, weights=CAD_se^-2)
  egger <- lm(CAD_beta ~ LDL_beta+HDL_beta+TG_beta, data=W, weights=CAD_se^-2)
  
  Q_snp <- ivw$resid^2*ivw$weights
  Q_1 <- sum(Q_snp)
  Q_2 <- sum(egger$resid^2*egger$weights)
  h2_1 <- Q_1/(nrow(cov)-3)
  h2_2 <- Q_2/(nrow(cov)-4)
  I[,i] <- c((h2_1-1)/h2_1,(h2_2-1)/h2_2)*100
  
  eggers[[i]] <- egger
}
