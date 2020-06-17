
library("TwoSampleMR")

# set exposures in parent call

pool <- data.frame()
for (i in 1:length(exposures)) {
  dat <- read.table(paste0("C:/Genetics/Variants/Significant/",exposures[i],"_SNPs.txt"), 
                    header = TRUE, stringsAsFactors = FALSE)
  print(paste0("Adding ", nrow(dat), " ", exposures[i], " trait SNPs..."))
  pool <- rbind(pool,dat)
}

pool <- pool[!duplicated(pool$SNP), ]
instruments <- clump_data(pool, clump_kb = 10000, clump_r2 = 0.01)
colnames(instruments) <- c("rsid","pval.selection","job")
