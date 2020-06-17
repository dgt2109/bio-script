
data <- read.table("C:/Genetics/SNP_Table/table_three-trait.txt", header = TRUE,
                   stringsAsFactors = FALSE)

outlier_data <- read.table("C:/Genetics/MRPRESSO/three-trait_MRPRESSO.txt", header = TRUE,
                           stringsAsFactors = FALSE)
rsid_MRPRESSO <- outlier_data[outlier_data$Pvalue < 0.05, 1]
data <- data[!(data$rsid %in% rsid_MRPRESSO), ]

# basic IVW MR
ivw <- lm(data$CAD_Beta ~ data$LDL_Beta + data$HDL_Beta + data$TG_Beta +0,
          weights = data$CAD_SE^-2)

# reorient for Egger
traits <- c("LDL","HDL","TG","CAD")
beta_indexes <- (1:length(traits)*3)-1
eggers <- list()
for (ref_col_index in beta_indexes[1:3]) {
  for (i in beta_indexes) {
    if (i == ref_col_index) { next }
    data[, i] <- ifelse(data[, ref_col_index]>0,
                        data[, i],
                        data[, i]*-1)
  }
  data[, ref_col_index] <- abs(data[, ref_col_index])
  
  eggers[[(ref_col_index+1)/3]] <- lm(data$CAD_Beta ~ data$LDL_Beta + data$HDL_Beta + 
                                        data$TG_Beta, weights = data$CAD_SE^-2)
}
names(eggers) <- colnames(data)[beta_indexes][1:3]

# compile results in output table
lipids <- c("LDL","HDL","TG")
bound_const <- qt(0.975,df=nrow(data)-length(lipids))
h2 <- sum(ivw$resid^2*ivw$weights)/(nrow(data)-3)
I2H <- (h2-1)/h2*100
ivw_table <- cbind(exp(ivw$coef),
                   exp(ivw$coef - bound_const*summary(ivw)$coef[,2]),
                   exp(ivw$coef + bound_const*summary(ivw)$coef[,2]),
                   summary(ivw)$coef[,4],
                   I2H,
                   nrow(data),
                   NA,NA)
colnames(ivw_table) <- c("OR", "LB", "UB", "P", "I2H", "N.Var","Int","Int.P")
rownames(ivw_table) <- lipids

egger_table <- vector()
bound_const <- qt(0.975,df=nrow(data)-length(lipids)-1)
for (i in 1:length(lipids)) {
  egger <- eggers[[i]]
  h2 <- sum(egger$resid^2*egger$weights)/(nrow(data)-4)
  I2H <- (h2-1)/h2*100
  eggers_row <- c(exp(egger$coef[i+1]),
                  exp(egger$coef[i+1] - bound_const*summary(egger)$coef[i+1,2]),
                  exp(egger$coef[i+1] + bound_const*summary(egger)$coef[i+1,2]),
                  summary(egger)$coef[i+1,4],
                  I2H,
                  nrow(data),
                  egger$coef[1],
                  summary(egger)$coef[1,4])
  egger_table <- rbind(egger_table, eggers_row)
}
colnames(egger_table) <- c("OR", "LB", "UB", "P", "I2H", "N.Var", "Int", "Int.P")
rownames(egger_table) <- lipids

results_table <- rbind(ivw_table, egger_table)
write.table(results_table, "C:/Genetics/Results/three-trait_MRPRESSO-results.txt",
            sep='\t', row.names=TRUE, col.names=NA)

# ivw <- lm(data$CAD_Beta ~ data$LDL_Beta + data$HDL_Beta + data$TG_Beta +0,
#           weights = data$CAD_SE^-2)
# Q_snp <- ivw$resid^2*ivw$weights
# p_CQ <- pchisq(Q_snp, df=1, lower.tail=FALSE)
# rsid_CochransQ <- data[p_CQ < 0.05/nrow(data), 1]
# data <- data[!(data$rsid %in% rsid_CochransQ), ]
