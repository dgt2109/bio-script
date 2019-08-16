library(MendelianRandomization)
library(MRPRESSO)

base_dir="C:/Users/David/Desktop/genetics/revision/snp680/"
file_in="table_680.txt"
data <- read.delim(paste(base_dir, file_in, sep=""),header=TRUE)
#data <- data[-429,]

W <- cbind(data$LDL_beta, data$TG_beta)
err <- cbind(data$LDL_se, data$TG_se)
x <- data$HDL_beta
xerr <- data$HDL_se
y <- data$CAD_beta
yerr <- data$CAD_se
  
MRMVInputObject_1 <- mr_mvinput(bx = W,
                                bxse = err,
                                by = y,
                                byse = yerr)
MRMVInputObject_2 <- mr_mvinput(bx = W, 
                                bxse = err, 
                                by = x, 
                                byse = xerr)

MRMVObject_1<- mr_mvivw(MRMVInputObject_1)
MRMVObject_2<- mr_mvivw(MRMVInputObject_2)
adjusted_y <- y - MRMVObject_1$Estimate[1]*W[,1] - MRMVObject_1$Estimate[2]*W[,2]
adjusted_x <- x - MRMVObject_2$Estimate[1]*W[,1] - MRMVObject_2$Estimate[2]*W[,2]

adjusted_y[adjusted_x < 0] <- -adjusted_y[adjusted_x < 0]
adjusted_x[adjusted_x < 0] <- -adjusted_x[adjusted_x < 0]

MRInputObject <- mr_input(bx = adjusted_x, bxse = xerr, by = adjusted_y,
                          byse = yerr, exposure="HDL-C", outcome="CAD")

IVWObject <- mr_ivw(MRInputObject, model = "default", robust = FALSE,
                    penalized = FALSE, correl = FALSE, weights = "simple",
                    psi = 0, distribution = "normal", alpha = 0.05)

EggerObject <- mr_egger(MRInputObject, robust = FALSE, penalized = FALSE,
                        correl = FALSE, distribution = "normal", alpha = 0.05)

SummaryStats <- data.frame(cbind(adjusted_x,xerr,adjusted_y,yerr))
colnames(SummaryStats) <- c("HDL","HDL_err","CAD","CAD_err")
res_presso <- mr_presso(BetaOutcome = "CAD",
                        BetaExposure = "HDL",
                        SdOutcome = "CAD_err",
                        SdExposure = "HDL_err",
                        OUTLIERtest = TRUE,
                        DISTORTIONtest = TRUE,
                        data = SummaryStats,
                        NbDistribution = 100000, 
                        SignifThreshold = 0.05)
ot <- res_presso$"Outlier Test"
ex <- ot[ot$Pvalue < 0.05, ]

