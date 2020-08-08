library(rmeta)

data <- read.table("C:/Genetics/SNP_Table/table_three-trait.txt", header = TRUE,
                   stringsAsFactors = FALSE)
param <- read.table("C:/Genetics/Results/three-trait_MRPRESSO-results.txt", header=TRUE, stringsAsFactors=FALSE)

CAD_dataset <- read.table("C:/Genetics/Datasets/CAD_data.txt",header=TRUE,stringsAsFactors=FALSE)
chrpos <- CAD_dataset[, c("rsid","Chr","Pos")]
data <- merge(data, chrpos, by="rsid")

loci <- data.frame(Gene=c("LIPG","LIPG N396S","LCAT","LPL","MLXIPL","CETP","LIPC","FADS","APOC3","ANGPTL4"),
                   Chr=c(18,18,16,8,7,16,15,11,11,19),
                   Pos=c(47087069,47087069,67978656,19796582,73038903,56995835,58702953,61583675,116700624,8429011),
                   Analysis=c("HDL","HDL","HDL","TG","TG","HDL","HDL","TG","TG","TG"))

SNPs <- list()
res_all <- list()
I2_all <- vector()
for (i in 1:nrow(loci)) {
  locus <- data[data$Chr == loci[i,2] & data$Pos > loci[i,3] - 5e5 & data$Pos < loci[i,3] + 5e5, ]
  if (i == 1) { locus <- locus[-c(3), ] }
  if (i == 2) { 
    locus <- rbind(locus[3, ],locus[3, ],locus[3, ])
    OR_Voight <- log(c(0.99, 0.88, 1.11))
    SE_Voight <- mean(OR_Voight[3]-OR_Voight[1], OR_Voight[1]-OR_Voight[2])/1.96
    locus[1,c("CAD_Beta","CAD_SE","CAD_P")] <- c(OR_Voight[1]*-1, SE_Voight, 0.41) # add Voight CAD data
    locus[2,c("CAD_Beta","CAD_SE","CAD_P")] <- c(0.101958,0.0508976,0.0451558)# add Nikpay CAD data
  }
  if(loci[i,4] == "HDL") { 
    res <- lm(locus$CAD_Beta ~ locus$HDL_Beta +0, weights = locus$CAD_SE^-2)
    cad_resid <- locus$CAD_Beta - log(param[1,2])*locus$LDL_Beta - log(param[3,2])*locus$TG_Beta
    res_adj <- lm(cad_resid ~ locus$HDL_Beta +0, weights = locus$CAD_SE^-2) 
  }
  if (loci[i,4] == "TG") {
    res <- lm(locus$CAD_Beta ~ locus$TG_Beta +0, weights = locus$CAD_SE^-2)
    cad_resid <- locus$CAD_Beta - log(param[1,2])*locus$LDL_Beta - log(param[2,2])*locus$HDL_Beta
    res_adj <- lm(cad_resid ~ locus$TG_Beta +0, weights = locus$CAD_SE^-2)    
  }
  locus$CAD_adj <- cad_resid
  h2 <- sum(res$resid^2*res$weights)/(nrow(locus)-1)
  I2H <- (h2-1)/h2*100
  I2_all[i] <- I2H
  SNPs[[i]] <- locus
  res_all[[i]] <- list(res,res_adj)
}

# single SNP MR
OR_all <- list()
for (i in 1:nrow(loci)) {
  SNP_list <- SNPs[[i]]
  OR <- data.frame(matrix(nrow=6,ncol=nrow(SNP_list)+1))
  rownames(OR) <- c("OR", "Upper", "Lower", "Adj OR", "Adj Upper", "Adj Lower")
  colnames(OR) <- c(SNP_list$rsid, paste0(nrow(SNP_list)," ", loci$Gene[i], " SNPs"))
  if (i == 8 ) { colnames(OR)[ncol(OR)] <- "1 FADS1 SNP" }
  for (j in 1:nrow(SNP_list)) {
    SNP <- SNP_list[j, ]
    if (loci[i,4] == "HDL") {
      beta_iv <- SNP$CAD_Beta / SNP$HDL_Beta
      beta_adj <- SNP$CAD_adj / SNP$HDL_Beta
      se_iv <- sqrt( SNP$CAD_SE^2/SNP$HDL_Beta^2 + SNP$CAD_Beta^2/SNP$HDL_Beta^4*SNP$HDL_SE^2 )
    }
    if (loci[i,4] == "TG") {
      beta_iv <- SNP$CAD_Beta / SNP$TG_Beta
      beta_adj <- SNP$CAD_adj / SNP$TG_Beta
      se_iv <- sqrt( SNP$CAD_SE^2/SNP$TG_Beta^2 + SNP$CAD_Beta^2/SNP$TG_Beta^4*SNP$TG_SE^2 )
    }
    # use CAD data to extract proper # of df for each SNP
    
    ############################
    OR[1, j] <- exp(beta_iv) # mean
    OR[2, j] <- exp(beta_iv + 1.96*se_iv) # upper bound
    OR[3, j] <- exp(beta_iv - 1.96*se_iv) # lower bound
    OR[4, j] <- exp(beta_adj) # mean, adj
    OR[5, j] <- exp(beta_adj + 1.96*se_iv) # upper bound, adj
    OR[6, j] <- exp(beta_adj - 1.96*se_iv) # lower bound, adj
  }
  OR[1,nrow(SNP_list)+1] <- exp(res_all[[i]][[1]]$coef)
  OR[2,nrow(SNP_list)+1] <- exp(res_all[[i]][[1]]$coef + qt(0.975,df=res_all[[i]][[1]]$df.residual) * summary(res_all[[i]][[1]])$coef[2])
  OR[3,nrow(SNP_list)+1] <- exp(res_all[[i]][[1]]$coef - qt(0.975,df=res_all[[i]][[1]]$df.residual) * summary(res_all[[i]][[1]])$coef[2])
  OR[4,nrow(SNP_list)+1] <- exp(res_all[[i]][[2]]$coef)
  OR[5,nrow(SNP_list)+1] <- exp(res_all[[i]][[2]]$coef + qt(0.975,df=res_all[[i]][[2]]$df.residual) * summary(res_all[[i]][[2]])$coef[2])
  OR[6,nrow(SNP_list)+1] <- exp(res_all[[i]][[2]]$coef - qt(0.975,df=res_all[[i]][[2]]$df.residual) * summary(res_all[[i]][[2]])$coef[2])

  OR_all[[i]] <- OR  
}

fig_width = c(4,8,4,4,4,4,4,5,4,4)
fig_height = c(2.25,2.25,2.25,3.25,2.25,5.5,3,2,4,4)
axislabel = paste0("OR CAD per s.d. ", loci$Analysis)
# forest plot
for (i in 1:nrow(loci)) {
  print(i)
  OR <- data.matrix(OR_all[[i]])
  tiff(paste0("C:/Genetics/Figures/",loci[i,1],".tif"),width=fig_width[i],height=fig_height[i],units="in",res=500)
  if (i == 2) {
    OR <- OR[, 1:ncol(OR)-1]
    labels <- cbind(c("",paste0(rep("LIPG N396S"),c(", Voight et al. meta-analysis", ", CARDIoGRAMplusC4D", ", CARDIoGRAMplusC4D+UKBB CAD"))),
                    c("n Cases","20,913","60,801","122,733"),
                    c("n Controls","95,407","123,504","424,528"))
    forestplot(labels, c(NA,log(OR[1,])), c(NA,log(OR[3,])), c(NA,log(OR[2,])), is.summary=c(TRUE,rep(FALSE,ncol(OR))),
               clip=c(-1.0,1.6), xlog=T, xlab=axislabel[i])  
  } else {
    if (ncol(OR) > 2) { # reorder columns
      OR <- cbind(OR[,1:ncol(OR)-1][,order(OR[,1:ncol(OR)-1][1,])],OR[,ncol(OR)])
      colnames(OR)[ncol(OR)] <- colnames(OR_all[[i]])[ncol(OR)]
    }
    labels <- cbind(colnames(OR),rep("",ncol(OR)))
    forestplot(labels, log(OR[1,]), log(OR[3,]), log(OR[2,]), is.summary=c(rep(FALSE,ncol(OR)-1),TRUE),
               clip=c(-0.7,1.6), xlog=T, xlab=axislabel[i])
  }
  dev.off()
#  forestplot(labels, OR[4,], OR[6,], OR[5,], is.summary=c(rep(FALSE,ncol(OR)-1),TRUE),clip=c(0,3))
}
