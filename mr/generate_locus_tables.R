
# inherit list of arrays SNPs 1 x 11
# inherit list of trait-of-interest loci$Analysis
# read SNP table 343x6 for crossref data

data <- read.table("C:/Genetics/SNP_Table/343x6.txt",header=T,stringsAsFactors=F)
param <- read.table("C:/Genetics/Results/sixtrait_MRPRESSO-results.txt",header=T,stringsAsFactors=F,
                    row.names=NULL)
rownames(param) <- make.unique(param[,1])
param <- param[,-1]
betas <- log(param[,1])[1:6]

crossref <- list()
for (i in 1:length(SNPs)) {
  crossref[[i]] <- merge(SNPs[[i]][,c(1,14:15)], data)
  crossref[[i]]$CAD_adj <- crossref[[i]]$CAD_Beta - betas[1]*crossref[[i]]$LDL_Beta - 
                                                    betas[4]*crossref[[i]]$BMI_Beta - 
                                                    betas[6]*crossref[[i]]$T2D_Beta
  crossref[[i]]$CAD_adj <- crossref[[i]]$CAD_adj - ifelse(is.na(crossref[[i]]$SBP_Beta),0,
                                                    betas[5]*crossref[[i]]$SBP_Beta)
  beta_indexes <- 1:8*3+1
  if(loci$Analysis[i] == "HDL") {
    crossref[[i]]$CAD_adj <- crossref[[i]]$CAD_adj - betas[3]*crossref[[i]]$TG_Beta
    crossref[[i]]$OR <- exp(crossref[[i]]$CAD_Beta / crossref[[i]]$HDL_Beta)
    crossref[[i]]$OR_adj <- exp(crossref[[i]]$CAD_adj / crossref[[i]]$HDL_Beta)
    
    ref_col_index = beta_indexes[2]
    for (j in beta_indexes) {
      if (j == ref_col_index) { next }
      crossref[[i]][, j] <- ifelse(crossref[[i]][, ref_col_index]>0,
                                   crossref[[i]][, j],
                                   crossref[[i]][, j]*-1)
    }
    crossref[[i]][, ref_col_index] <- abs(crossref[[i]][, ref_col_index])
  }
  else if (loci$Analysis[i] == "TG") {
    crossref[[i]]$CAD_adj <- crossref[[i]]$CAD_adj - betas[2]*crossref[[i]]$HDL_Beta
    crossref[[i]]$OR <- exp(crossref[[i]]$CAD_Beta / crossref[[i]]$TG_Beta)
    crossref[[i]]$OR_adj <- exp(crossref[[i]]$CAD_adj / crossref[[i]]$TG_Beta)
    
    ref_col_index = beta_indexes[3]
    for (j in beta_indexes) {
      if (j == ref_col_index) { next }
      crossref[[i]][, j] <- ifelse(crossref[[i]][, ref_col_index]>0,
                                   crossref[[i]][, j],
                                   crossref[[i]][, j]*-1)
    }
    crossref[[i]][, ref_col_index] <- abs(crossref[[i]][, ref_col_index])
  }
  crossref[[i]] <- crossref[[i]][,-c(2:3,1:7*3+2)]
  crossref[[i]] <- crossref[[i]][order(crossref[[i]]$OR),]
  crossref[[i]]$Locus <- loci$Gene[i]
  
  write.table(crossref[[i]], "C:/Genetics/locus_tables-3.txt", append=T, sep='\t', row.names=F)
}
