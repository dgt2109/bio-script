
library("MRPRESSO")

data <- read.delim("C:/Genetics/SNP_Table/table_three-trait.txt", 
                   header=TRUE, stringsAsFactors = FALSE)
  
SummaryStats <- data.frame(cbind(data$LDL_Beta,
                                 data$LDL_SE,
                                 data$HDL_Beta,
                                 data$HDL_SE,
                                 data$TG_Beta,
                                 data$TG_SE,
                                 data$CAD_Beta,
                                 data$CAD_SE))
  
colnames(SummaryStats) <- c("E1_effect",
                            "E1_se",
                            "E2_effect",
                            "E2_se",
                            "E3_effect",
                            "E3_se",
                            "Y_effect",
                            "Y_se")
  
res <- mr_presso(BetaOutcome = "Y_effect", 
                 BetaExposure = c("E1_effect", "E2_effect","E3_effect"), 
                 SdOutcome = "Y_se", 
                 SdExposure = c("E1_se", "E2_se","E3_se"), 
                 OUTLIERtest = TRUE, 
                 DISTORTIONtest = FALSE, 
                 data = SummaryStats, 
                 NbDistribution = 3.5e5, 
                 SignifThreshold = 0.05)
  
outlier_test <- res$`MR-PRESSO results`$`Outlier Test`
presso_output <- cbind(data$rsid, outlier_test)
  
write.table(presso_output, "C:/Genetics/MRPRESSO/three-trait_MRPRESSO.txt", 
            sep='\t', row.names=TRUE, col.names=NA)

# P_achieved = N.Var / NbDist 
# T = N.Var * NbDist * (0.00023 sec)
