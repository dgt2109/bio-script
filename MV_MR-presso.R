library("MRPRESSO")

base_dir="C:/Users/David/Desktop/genetics/revision/snp680/"
file_in="table_680_wcadj3.txt"
#base_dir="C:/Users/David/Desktop/genetics/data/GW/"
#file_in="table_179.txt"

data <- read.delim(paste(base_dir, file_in, sep=""),header=TRUE)

SummaryStats <- data.frame(cbind(data$LDL_beta,
                                 data$LDL_se,
                                 data$HDL_beta,
                                 data$HDL_se,
                                 data$TG_beta,
                                 data$TG_se,
                                 data$CAD_beta,
                                 data$CAD_se))

colnames(SummaryStats) <- c("E1_effect",
                            "E1_se",
                            "E2_effect",
                            "E2_se",
                            "E3_effect",
                            "E3_se",
                            "Y_effect",
                            "Y_se")
# sample code
res <- mr_presso(BetaOutcome = "Y_effect", 
          BetaExposure = c("E1_effect", "E2_effect","E3_effect"), 
          SdOutcome = "Y_se", 
          SdExposure = c("E1_se", "E2_se","E3_se"), 
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = SummaryStats, 
          NbDistribution = 128000, 
          SignifThreshold = 0.05)

ot <- res$`MR-PRESSO results`$`Outlier Test`
ex <- ot[ot$Pvalue < 0.05, ]

write.csv(ot,paste0(base_dir,"MV_outliertest_wcadj3_128000.csv"))
write.csv(ex,paste0(base_dir,"MV_outlier_wcadj3_128000.csv"))

# 1000 * 2^4 = 16000 simulations for p <0.05 ~3 hr
# 1000 * 2^7 = 128000 simulations for p <0.0005 ~12 hr
