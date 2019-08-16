base_dir="C:/Users/david/Desktop/genetics/rerevision/pleiotropy_886/"

data <- read.delim(paste0(base_dir,"table_886p.txt"))

rownames(data) <- data[,1]
data <- data <- data[,-1]
betas <- data[,c(1,5,9,13,17,21,25)]

betas$LDL_beta <- ifelse(betas$HDL_beta > 0, betas$LDL_beta, betas$LDL_beta * -1)
betas$TG_beta <- ifelse(betas$HDL_beta > 0, betas$TG_beta, betas$TG_beta * -1)
betas$CAD_beta <- ifelse(betas$HDL_beta > 0, betas$CAD_beta, betas$CAD_beta * -1)
betas$BMI_beta <- ifelse(betas$HDL_beta > 0, betas$BMI_beta, betas$BMI_beta * -1)
betas$T2D_beta <- ifelse(betas$HDL_beta > 0, betas$T2D_beta, betas$T2D_beta * -1)
betas$SBP_beta <- ifelse(betas$HDL_beta > 0, betas$SBP_beta, betas$SBP_beta * -1)

betas$HDL_beta <- abs(betas$HDL_beta)
weights <- data$CAD_se^-2

ivw <- lm(betas$CAD_beta ~ betas$HDL_beta + betas$LDL_beta + betas$TG_beta +
          betas$T2D_beta + betas$SBP_beta +0)
egger <- lm(betas$CAD_beta ~ betas$HDL_beta + betas$LDL_beta + betas$TG_beta +
            betas$T2D_beta + betas$SBP_beta )